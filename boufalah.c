#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


#define NEUTREADDICTION 0
#define NEUTREMAX -2147483647
#define MAXIMUM 100
#define SUM 1
#define MAX 2
#define PREFIX 1
#define SUFIX 2

#define max(a,b) (((a)>(b))?(a):(b))

struct tablo {
  int * tab;
  int size;
};


int operator(int value1, int value2, int op){
  switch(op){
    case 1 : // sum
      return value1 + value2;
    case 2 : // max
      return max(value1 ,value2);
    default :
      exit(-1);
  }
}

void printArray(struct tablo * tmp) {
  printf("---- Array of size %i ---- \n", tmp->size);
  int size = tmp->size;
  int i;
  for (i = 0; i < size; ++i) {
    printf("%i ", tmp->tab[i]);
  }
  printf("\n");
}

struct tablo * allocateTablo(int size) {
  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->size = size;
  tmp->tab = calloc(size,sizeof(int));
  return tmp;
}



void montee(struct tablo * source, struct tablo * destination, int op) {
 for (int i = log2(source->size)-1; i > -1; i--) {
   int fin = pow(2,i+1);
   #pragma omp parallel for
    for (int j = fin/2; j < fin ; j++) {
      destination->tab[j] = operator(destination->tab[2*j],destination->tab[2*j+1],op);
    }
  }
}


void descente(struct tablo * a, struct tablo * b, int op) {
  // implementation de la phase de descente
  switch(op){
    case 1 : // sum
      b->tab[1] = NEUTREADDICTION;
      break;
    case 2 : // MAX
      b->tab[1] = NEUTREMAX;
      break;
  }
  for (int i = 1; i <= log2(b->size/2); i++) {
    int fin =pow(2,i+1);
    #pragma omp parallel for
    for (int j = fin/2; j < fin ; j++) {
      if (j%2) {
          //impair
          b->tab[j] = operator(b->tab[j/2],a->tab[j-1],op);
      }else{
        //pair
          b->tab[j] = b->tab[j/2];
      }
    }
  }
}

void final(struct tablo * a, struct tablo *b, int op) {
  //phase finale
  #pragma omp parallel for
  for (int i =  b->size/2; i < b->size; i++) {
    b->tab[i] = operator(b->tab[i],a->tab[i],op);
  }
}


void initPrefix(struct tablo* source,struct tablo* a){
  int firstCase = a->size - source->size;
  #pragma omp parallel for
  for (int i = 0; i < source->size; i++) {
    a->tab[i+firstCase] = source->tab[i];
  }
}


void initSufix(struct tablo* source,struct tablo* a){
  int firstCase = a->size - source->size;
  #pragma omp parallel for
  for (int i = 0; i < source->size; i++) {
    a->tab[i+firstCase] = source->tab[source->size-i-1];
  }
}

void initScan(struct tablo* source,struct tablo* a, int type){
  switch(type){
    case 1 : // prefix
      initPrefix(source, a);
      break;
    case 2 : // Sufix
      initSufix(source, a);
      break;
  }
}

void resultSufix(struct tablo* b,struct tablo* result){
  int deb = result->size*2;
  #pragma omp parallel for
  for (int i =0 ; i < result->size; i++) {
    result->tab[i] = b->tab[deb-i-1];
  }
}

void resultPrefix(struct tablo* b,struct tablo* result){
  int deb = result->size;
  #pragma omp parallel for
  for (int i = deb; i < b->size; i++) {
    result->tab[i-deb] = b->tab[i];
  }
}



void resultScan(struct tablo* source,struct tablo* b, int type){
  switch(type){
    case 1 : // prefix
      resultPrefix(source, b);
      break;
    case 2 : // Sufix
      resultSufix(source, b);
      break;
  }
}

void scan(struct tablo* source,struct tablo* result,int operator ,int type){
	// type = [prefix 1|sufix |] // operator [1 sum | 2 max]  
  struct tablo * a = allocateTablo(source->size*2);
  struct tablo * b = allocateTablo(source->size*2);


  initScan(source, a, type); // type = prefix or sufix

  montee(source, a,operator);
  descente(a, b,operator);
  final(a,b,operator);

  resultScan(b,result,type);

  free(a->tab);
  free(a);
  free(b->tab);
  free(b);
}



void submax(struct tablo * source, struct tablo * pMax , struct tablo * sSum,struct tablo * sMax , struct tablo * pSum, struct tablo * m){
  #pragma omp parallel for
  for (int i = 0; i < source->size; i++) {
    m->tab[i] = pMax->tab[i] - sSum->tab[i] + source->tab[i] + sMax->tab[i] - pSum->tab[i];
  }
}

void resultMaxSubArray(struct tablo * m, struct tablo * source ){
  // fin first max
  int deb = 0 ;
  int fin = 0 ;
  int endofsub = 0;
  for (int i = 1; i < m->size; i++) {
    if (m->tab[deb] < m->tab[i] ) {
      deb = i;
      fin = i;
      endofsub = 0;

    }else{
	    if (m->tab[deb] == m->tab[i] && endofsub ==0) {
	      fin = i;
	    }else{
	    	endofsub = 1;
	    }
    }
  }
  printf("%d ",m->tab[deb]);
  for (int i = deb; i <= fin; i++) {
    printf("%d ",source->tab[i]);
  }
  printf("\n");
}

void maxSubArray(struct tablo* source){
  struct tablo * pSum = allocateTablo(source->size);
  struct tablo * sSum = allocateTablo(source->size);
  #pragma omp parallel
  {
    scan(source,pSum,SUM,PREFIX); // prefixSum
    scan(source,sSum,SUM,SUFIX); // sufixSum
  }

  //printf("------------------ sufixSum ------------------\n");
  //printArray(sSum);

  //printf("------------------ prefixSum ------------------ \n");
  //printArray(pSum);

  struct tablo * sMax = allocateTablo(source->size);
  struct tablo * pMax = allocateTablo(source->size);
  #pragma omp parallel
  {
  scan(sSum, pMax,MAX,PREFIX);
  scan(pSum, sMax,MAX,SUFIX);
  }
  //printf("------------------ sufixMax ------------------ \n");
  //printArray(sMax);

  //printf("------------------ prefixMax ------------------ \n");
  //printArray(pMax);

  struct tablo * m = allocateTablo(source->size);

  submax(source,pMax ,sSum,sMax ,pSum,m); // dernier étape
  //printf("------------------ Result ------------------ \n");
  //printArray(m);

  resultMaxSubArray(m, source);

  free(sSum->tab);
  free(sSum);
  free(pSum->tab);
  free(pSum);
  free(pMax->tab);
  free(pMax);
  free(sMax->tab);
  free(sMax);
  free(m->tab);
  free(m);
}

void generateArray(struct tablo * s) {
  //construction d'un tableau pour tester
  s->size=8; //1048576;
  s->tab=malloc(s->size*sizeof(int));
  s->tab[0]=-1;
  s->tab[1]=2;
  s->tab[2]=-3;
  s->tab[3]=6;
  s->tab[4]=-1;
  s->tab[5]=8;
  s->tab[6]=-10;
  s->tab[7]=-5;



  /*s->size=16;
  s->tab=malloc(s->size*sizeof(int));
  s->tab[0]=3;
  s->tab[1]=2;
  s->tab[2]=-7;
  s->tab[3]=11;
  s->tab[4]=10;
  s->tab[5]=-6;
  s->tab[6]=4;
  s->tab[7]=9;
  s->tab[8]=-6;
  s->tab[9]=1;
  s->tab[10]=-2;
  s->tab[11]=-3;
  s->tab[12]=4;
  s->tab[13]=-3;
  s->tab[14]=0;
  s->tab[15]=2;*/
}

void generateArrayRandom(struct tablo * s, int size) {
  double taille = log2(size);
  int m  =(int) taille;
  if(taille > (double)m){
    size = pow(2,m+1);
  }
  s->size = size;
  srand(time(NULL));
  s->size = size;
  s->tab=malloc(s->size*sizeof(int));
  for (int i = 0; i < s->size; i++) {
    s->tab[i] =  (rand()%MAX);
  }
}

void checkArray(struct tablo* s){
  double size = log2(s->size);
  int m  =(int) size;
  if(size > (double)m){
    // add element neutre
    int newsize = pow(2,m+1);
    s->tab =realloc( s->tab , newsize * sizeof(int) );
    for (int i = s->size; i < newsize; i++) {
      s->tab[i] = 0;
    }
    s->size =newsize;
  }
}

void readArray(struct tablo * source, char* file){
  int number;
  FILE* in_file = fopen(file, "r");
  if (! in_file ){
   exit(-1);
  }
  fseek(in_file, 0L, SEEK_END);
  source->size = ftell(in_file);
  rewind(in_file);
  int i= 0;
  source->tab = malloc(source->size* sizeof(int));
  while ( fscanf(in_file, "%d", & number ) == 1 ) {
      source->tab[i]=number;
      i++;
  }
  source->size = i;
  fclose(in_file);
}


int main(int argc, char **argv) {

  struct tablo source;
  if(argc > 1){
      readArray(&source, argv[1]);
  }else{
     generateArray(&source);
  }
  //checkArray(&source); // verifie la taille du tablo 2^N
  //printArray(&source);

  maxSubArray(&source);
// desalocation de memoire
  free(source.tab);
}
