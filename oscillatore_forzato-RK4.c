#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define beta 0.5
#define k1 1.
#define k2 2./3.

struct motion {
  double a;
  double w;
  double t;
};

struct K {
  double x;
  double v;
};

//Prototipi
struct motion *cauchy(double a0, double w0);
struct motion *pickmethod(struct motion*old, double dt, double fzero, int n);
struct motion *eulerocromer(struct motion* old, double dt, double fzero);
struct motion *RK4(struct motion* old, double dt, double fzero);
double force(struct motion* phasespace, double fzero);

int main(int argc, char** argv) {
double dt = strtod(argv[1], NULL);
double f0,Tmax = 5*6280*dt, a0 = (M_PI*0.5), w0 = 0.;
struct motion awt;
 int i,j, n = 0;
char nomefile[100];
FILE* file;

if( argc != 2 ) {
  fprintf(stderr, "Metti il dt!\n");
  exit(1);
}
while( n < 2 ) {
  for( j = 1; j >= 0; j-- ) {
    for( i = 0; i < 5; i++ ) {
      if( j ) {
        awt =* cauchy(a0,w0); 
      } else {
        awt = *cauchy((a0-0.001),(w0+0.001));
      }

      sprintf(nomefile, "metodo%i.f0_%i.%i.dat",n , i,j);
      file = fopen(nomefile, "w");
      fprintf(file, "%.14e %.14e %.14e\n", awt.a,awt.w,awt.t);
      if( i == 0 ) {
        f0 = 0.9;
      } else if ( i == 1 ) {
        f0 = 1.07;
      } else if ( i==2 ) {
        f0 = 1.15;
      } else if( i == 3 ) {
        f0 = 1.47;
      } else {
        f0 = 1.5;
      }
      while( awt.t <= Tmax ) {
        awt = *pickmethod(&awt, dt, f0, n);
        fprintf(file, "%.14e %.14e %.14e\n",awt.a, awt.w, awt.t);
      }
     }
  }
  n++;
 }
  fclose(file);
  return 0;
}

struct motion *pickmethod( struct motion*old, double dt, double fzero, int n ) {
  if( n ) {
    return RK4(old,dt,fzero);
  } else {
    return eulerocromer(old,dt, fzero);
  }
}
struct motion *eulerocromer(struct motion* old, double dt, double fzero){
  static struct motion new;
  new.w = old->w+force(old,fzero)*dt;
  new.a = old->a+new.w*dt;
  new.t = old->t+dt;
  return &new;
}
struct motion* cauchy(double a0, double w0){
  static struct motion start;
  start.a=a0;
  start.w=w0;
  start.t=0.;
  return &start;
}
struct motion *RK4(struct motion* old, double dt, double fzero){
  struct K K1,K2,K3,K4;
  static struct motion new;
  new.a=old->a;
  new.w=old->w;
  new.t=old->t+dt*0.5;
//yn=(xn, vn)
//gn=(vn,fn(x,v,t))
//fn=-k*sin(an)-b*vn-f0*cos(w*tn)
  //K1=g(tn,yn)*dt
  K1.x = old->w*dt;//è una posizione
  K1.v = force(old, fzero)*dt; // è una velocità
  //K2=g(tn*dt/2,yn+k1/2)*dt
  new.a = old->a+0.5*K1.x;
  new.w = old->w+0.5*K1.v;
  K2.x = new.w*dt; // posizione
  K2.v = force(&new, fzero)*dt; //velocità
  //K3=g(tn+dt/2,yn+k2/2)*dt
  new.a = old->a+0.5*K2.x;
  new.w = old->w+0.5*K2.v;
  K3.x = new.w*dt; //posizione
  K3.v = force(&new, fzero)*dt;//velocità
  //K4=g(tn+dt, yn+k3)*dt
  new.a = old->a+K3.x;
  new.w = old->w+K3.v;
  new.t = old->t+dt;
  K4.x = new.w*dt;//posizione
  K4.v = force(&new, fzero)*dt;// velocità
  //ora calcoliamo il passo successivo
  new.a = old->a + 1./6.*(K1.x+2*K2.x+2*K3.x+K4.x);
  new.w = old->w + 1./6.*(K1.v+2*K2.v+2*K3.v+K4.v);
  return &new;
}
double force(struct motion* phasespace,double fzero){
  return -k1*sin(phasespace->a)-beta*phasespace->w+fzero*cos(k2*phasespace->t);
}
