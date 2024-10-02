#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define gamma0 0.5
#define w_sq 1.
#define x0 1.
#define v0 0.

#define Tmax 200. //periodo dopo il quale consideriamo il moto critico
#define Amin 0.00000000000002
struct motion {
  double x;
  double v;
  double t;
};
struct point {
  double x;
  double y;
};

//Prototipi
struct motion *cauchy(void);
struct motion *verlet(struct motion* old, double dt, double gamma);
struct motion * countzeros(struct motion* old, double dt, double gamma, int* counter);
double max_fabs(double a, double b);
double force(double gamma, struct motion* phasespace);
double fit(struct point* p0, struct point* p1, double halfA);

int main(int argc, char** argv){
  int zeros = 0, i = 0, disclaimer;
  double dt = strtod(argv[1], NULL);
  double gamma = gamma0, A;
  double sx = A, dx = A;//intervallo
  double tau, oldtau = 0., maxmin, time = 0.;
  struct point in, fin;//estremi dell'intervallo su cui fare il fit
  struct motion xvt, before;
  FILE*fp;
  if(argc != 2){
    fprintf(stderr, "Metti il dt!\n");
  }

  //trovo gamma critico
  
  while(gamma == gamma0 || zeros == 2){
    zeros = 0;
    xvt = *cauchy();
    while(xvt.t < Tmax && zeros < 2){
      xvt = *countzeros(&xvt, dt, gamma, &zeros); //verlet incrementa il contatore di zeri
    }
    if(zeros == 2){
      gamma = 1.1*gamma;
    }
  }
  fprintf(stderr, "Gamma critico è %.14lf\n", gamma);
  fprintf(stdout,"\n\n");
  //studiamo il moto smorzato
  
  xvt = *cauchy();
  A = x0;
  fp = fopen("tau.dat", "w");
  dx = A;
  while(A > Amin){
    disclaimer = 1;
     while(A*0.5 < fabs(dx)){//Cerca un max e un min successivi che contengono A/2
      sx = dx;//...che comincia dove finiva quello "sbagliato"..
      //..e finisce nel massimo successivo che troviamo
      in.x = time;
      in.y = sx;
      xvt = *verlet(&xvt, dt, gamma0);//mi allontano dall'ultimo max che ho trovato
      while(disclaimer){
	before = xvt;
	xvt = *verlet(&xvt, dt, gamma0);
	if((xvt.v*before.v) < 0){//se sono in un massimo o in un minimo
	  disclaimer = 0;//eureka
	  maxmin = max_fabs(before.x,xvt.x);//il massimo/minimo è
	  if(maxmin == before.x){ //e lo troviamo al tempo
	    time = before.t;
	  } else {
	    time = xvt.t;
	  }
	  dx = maxmin;
	}
      }
    }
      fin.x = time;
      fin.y = dx;
      tau = fit(&in, &fin, A*0.5);
      A *= 0.5;//adesso dimezzo ancora
      fprintf(fp, "A %.14lf SX %.14lf T0 %.14lf DX %.14lf T1 %.14lf I %i TAU_I %.14lf DIFF %.14lf\n", A, in.y,in.x, fin.y,fin.x, i,tau, tau-oldtau);
      oldtau = tau;
      i++;
  }
  fclose(fp);
  return 0;
}
struct motion* cauchy(void){
  static struct motion start;
  start.x = x0;
  start.v = v0;
  start.t = 0;
  return &start;
}
struct motion *verlet(struct motion* old, double dt, double gamma){
  static struct motion new;
  new.x = old->x+old->v*dt+0.5*force(gamma, old)*dt*dt;
  new.v = (old->v+0.5*(force(gamma, old)-w_sq*new.x)*dt)/(1.+gamma*dt*0.5);
  new.t = old->t+dt;
  fprintf(stdout, "%.14lf %.14lf %.14lf\n", new.x, new.v, new.t);
  return &new;
}
struct motion* countzeros(struct motion* old, double dt, double gamma, int* counter){
  static struct motion new;
  new=*verlet(old,dt, gamma);
  if(old->x*new.x <= 0.){
    *counter += 1;
  } 
  return &new;
}
double force(double gamma, struct motion* phasespace){
  return -w_sq*phasespace->x-gamma*phasespace->v;
}
double fit(struct point* p0, struct point* p1, double halfA){
  double t;
  double m;
  m = (p1->y-p0->y)/(p1->x-p0->x);
  t = p1->x+(halfA-p1->y)/m;
  return t;
}
double max_fabs(double a, double b){
  double max_f;
  max_f = fabs(a)>fabs(b) ? a:b;
  return max_f;
}