#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Maxm 1024

typedef unsigned long int RANDOM; 

int maxT(RANDOM a, RANDOM b,RANDOM m);
int findT(RANDOM seed,RANDOM a, RANDOM b, RANDOM m);
int primes(unsigned long int * p, RANDOM m);
int coprimes(unsigned long int * p,int d, RANDOM a, RANDOM b);
int Atilde(unsigned long int * p, int d, RANDOM atilde,RANDOM m);


int main(int argc, char ** argv ){

  RANDOM a = 27;//(RANDOM) strtol(argc[1],NULL,10);
  RANDOM b = 7;//(RANDOM)strtol(argc[2],NULL,10);
  RANDOM m = 739;//(RANDOM)strtol(argc[3],NULL,10);
  RANDOM seed = 2;//(RANDOM)strtol(argc[4],NULL,10);
  int T,check;

  if(argc != 5){
    fprintf(stderr, " Inserire il moltiplicatore, l'incremento il modulo e il seme\n");
    exit (1);
  }

  T = findT(seed, a, b, m);
  check = maxT(a,b,m);
  if(check && (T == m)){
    fprintf(stderr, "I parametri del generatore soddisfano le condizioni per un massimo periodo\n il periodo calcolato e' %d\n", T);
  } else if((check == 0) && (T != m)){
    fprintf(stderr,"I parametri del generatore non soddisfano le condizioni per un massimo periodo\n il periodo calcolato e' %d\n", T);
  } else {
      fprintf(stderr, "Il programma e' sbagliato\n");
  }
  return 0;
}


int maxT(RANDOM a, RANDOM b,RANDOM m){
 unsigned long int *p;//array dei numeri primi tra 2 e m
 int d, check = 0, i;
 RANDOM atilde = a-1;

 p = (unsigned long int *) malloc(m * sizeof (unsigned long int));//voglio un array di m caselle.. posso?
 if(p == NULL){
   fprintf(stderr, "non e' stato possibile trovare i numeri primi minori di m \n");//no...
   exit(1);
 }
 d = primes(p, m);//dà il numero di primi tra 2 ed m e li mette nell'array p
 p = (unsigned long int *) realloc(p,d* sizeof(unsigned long int));//ridimensiono p.. posso?

 if(p == NULL){
   fprintf(stderr, "non e' stato possibile ridimensionare la malloc\n");
   exit(1);//no...
   }
 check = coprimes(p,d,b,m);// m e b sono coprimi??
 check += Atilde(p,d,atilde,m);//atilde è divisibile per tutti i fattori primi di m?

 if((m%4) == 0){
   if((atilde%4) == 0){
     check += 1;
   }
   if(check == 3){
     return 1;
   }else{
     return 0;
   }
 }else{ 
   if (check == 2){
     return 1;
   } else {
     return 0;
   }
 }
}

int primes(unsigned long int * p, RANDOM m){
  int d = 0;
unsigned long  int i, factor, maxfactor;
  for(i = 2;i <= m; i++){// 1 non e' primo
    factor = 2; // factor sono i possibili divisori di i
    maxfactor = (int)sqrt(i)+1;

    while((factor < maxfactor) && (i%factor)){
      factor++;
    }
    if(factor == maxfactor){//se nessuno dei possibili divisori e' un divisore
      *(p+d) = i;
      d++;
   }
  }
  return d;
}
int coprimes(unsigned long int * p, int d, RANDOM a, RANDOM b){
  int i;
  for(i = 0; i < d; i++){
    if(((a%(*(p+i))) == 0) && ((b%(*(p+i))) == 0)){
      return 0;//ho trovato un div comune quindi non coprimi
    }
  }
  return 1;// non ho trovato div comuni, quindi coprimi
}

int Atilde(unsigned long int * p, int d, RANDOM atilde,RANDOM m){
  int i,div_com = 0, div_m = 0;
  for(i=0; i<d; i++){
    if((m%(*(p+i))) == 0){
      div_m++;
    }
    if(((atilde%(*(p+i))) == 0) && ((m%(*(p+i))) == 0)){
      div_com ++;
    }
  }

  if(div_com == div_m) {
    return 1;
  } else {
    return 0;
  }
}
int findT(RANDOM seed,RANDOM a, RANDOM b, RANDOM m){
  RANDOM *tau;//numeri generati
  RANDOM r,In=seed;
  int i=1,stop=0, j;
  tau= calloc(m, sizeof(RANDOM));
  while(stop==0){
    r=(a * In+b)%m;//genera un random
    In=r; //cambia il seed
    fprintf(stdout, " %u %d \n", r, i);//stampalo
    if((*(tau+r)) == 0){
      *(tau+r)=i;//aggiorna l'array se non è mai stato generato
      i++;
    }else{
       stop=1;//altrimenti fermati
    }
  }
  return i - *(tau + r); 
}
