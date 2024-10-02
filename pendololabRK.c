#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 1 //siamo in un nuovo sistema di unità di misura dove g è 1, il pendolo ha inoltre  massa 1 e lunghezza 1
#define w0 0.
#define da 0.017 //decremento di amax, in rad è circa 1 grado
#define TArm 2.*M_PI/k
#define smallenough 0.001
#define EULER 0
#define EULER_CROMER 1
#define MIDPOINT 2
#define VERLET 3
#define RK4 4
struct motion {
  double a;
  double w;
  double t;
};

struct K {
  double x;
  double v;
};

struct point {
  double x;
  double y;
};

// Prototipi
//1) Forze in gioco
double force(double a); //F=-k*sin(a), k è definita sopra
double forceapprox(double a); //Piccoli angoli, F=-k*a
//2) Integrazione numerica
struct motion *pickmethod (struct motion* old, double dt,int algorithm, int f);
struct motion *euler(struct motion* old, double dt, int f);
struct motion *euler_cromer(struct motion* old, double dt, int f);
struct motion *midpoint(struct motion* old, double dt, int force);
struct motion *verlet(struct motion* old, double dt, int f);
struct motion *RungeKutta4(struct motion* old, double dt, int f);
//3) Energia
double energy(struct motion*awt);
//4) Condizioni inizial
struct motion *cauchy(double amax);// la velocità iniziale è nulla
struct motion *findamax(struct motion* old, double dt, double* tempi,int* pcounter);
//4) Interpolazione
double fit(struct point* p0, struct point* p1);
int main(int argc, char** argv){
  int counter, algorithm, disclaimer = 0, i, j;
  double a0 = strtod(argv[1], NULL);
  double dt0 = strtod(argv[2], NULL);
  double E, E0, T = 0., amax = a0, tempi[2], dt = dt0, Tmax = 6280*dt0;// L'array tempi contiene il tempo a cui a si azzera per la prima volta e quello dove a si azzera per la 21ma, se prendo dt0=0.01, in Tmax dovrei avere 10 oscillazioni
  char nomefile[100], file_ene[100];
  FILE* fp;
  FILE*ene;

  struct motion awt;
//Cerco amax
  while(fabs((TArm-T)/TArm) > smallenough){
    if(disclaimer){
      amax *= 0.9;
    } else {
      disclaimer=1;
    }
    awt = *cauchy(amax);
    counter=0;
    while(counter <= 21){
      awt = *findamax(&awt, dt, tempi, &counter);
    }
    T = (*(tempi+1)-*tempi)*0.1;
    printf("%.14lf %lf\n", T, amax);
  }
  E0 = energy(cauchy(amax));
  for(i = 0;i < 5; i++){
    if(i != 0){
      dt *= 0.5;
    }
//Ora integro numericamente
    for(algorithm = 0; algorithm < 5; algorithm++){
      sprintf(nomefile, "metodo%i.%lf.dat", algorithm,dt);
      fp = fopen(nomefile, "w");
      // sprintf(file_ene, "energie%i.%i.dat", algorithm, i);
      sprintf(file_ene, "energie%i.dat", algorithm);
      ene = fopen(file_ene, "a+");
      awt = *cauchy(amax);
      fprintf(fp, "%.14e %.14e %.14e %.14e\n", awt.a,awt.w,awt.t, E0);
      while(awt.t < Tmax){
        awt = *pickmethod(&awt, dt, algorithm, 1);
        E = energy(&awt);
        fprintf(fp, "%.14e %.14e %.14e %.14e\n", awt.a,awt.w,awt.t, E);
      }
      fprintf(ene,"%.14e %.14e %.14e\n", dt, fabs((E-E0)/E0), awt.t);
      fprintf(fp,"\n\n#Approssimazione con oscillatore armonico\n\n");
      awt = *cauchy(amax);
      fprintf(fp, "%.14e %.14e %.14e %.14e\n", awt.a,awt.w,awt.t, E0);
      while(awt.t < 10.*Tmax){
        awt = *pickmethod(&awt,dt,algorithm,0);
        E = energy(&awt);
        fprintf(fp, "%.14e %.14e %.14e %.14e\n", awt.a,awt.w,awt.t, E);
      }
      fclose(fp);
      fclose(ene);
    }
  }
  return 0;
}
double force(double a){
  return -k*sin(a);
}
double forceapprox(double a){
  return -k*a;
}
struct motion* pickmethod (struct motion* old, double dt,int algorithm, int force){
  if (algorithm == EULER){
    return euler(old, dt, force);
  } else if (algorithm == EULER_CROMER){
      return euler_cromer(old,dt, force);
    } else if (algorithm == MIDPOINT){
        return midpoint(old,dt, force);
    } else if (algorithm == VERLET){
          return verlet(old, dt, force);
   } else {
     return RungeKutta4(old, dt, force);
   }
}
struct motion * euler(struct motion* old, double dt, int f){
  static struct motion new;
  new.a = old->a+old->w*dt;
  if(f){
    new.w = old->w+force(old->a)*dt;
  } else {
    new.w = old->w+forceapprox(old->a)*dt;
  }
  new.t = old->t+dt;
  return &new;
}
struct motion * euler_cromer(struct motion* old, double dt, int f){
  static struct motion new;
  if(f){
    new.w = old->w+force(old->a)*dt;
  } else {
    new.w = old->w+forceapprox(old->a)*dt;
  }
  new.a = old->a+new.w*dt;
  new.t = old->t+dt;
  return &new;
}
struct motion * midpoint(struct motion* old, double dt, int f){
  static struct motion new;
  if(f){
    new.w = old->w+force(old->a)*dt;
  } else {
    new.w = old->w+forceapprox(old->a)*dt;
  }
  new.a = old->a+0.5*(new.w+old->w)*dt;
  new.t = old->t+dt;
  return &new;
}
struct motion *verlet(struct motion* old, double dt, int f){
  static struct motion new;
  if(f){
    new.a = old->a+old->w*dt+0.5*force(old->a)*dt*dt;
    new.w = old->w+0.5*(force(old->a)+force(new.a))*dt;
  } else {
    new.a = old->a+old->w*dt+0.5*forceapprox(old->a)*dt*dt;
    new.w = old->w+0.5*(forceapprox(old->a)+forceapprox(new.a))*dt;
  }
  new.t = old->t+dt;
  return &new;
}
struct motion *RungeKutta4(struct motion* old, double dt, int f){
  struct K K1,K2,K3,K4;
  static struct motion new;
//yn=(xn, vn)
//gn=(vn,fn(x,v,t))
//fn=-k*sin(an)-b*vn-f0*cos(w*tn)
  //K1=g(tn,yn)*dt
  K1.x = old->w*dt;//è una posizione
  if(f) {
    K1.v = force(old->a)*dt;
  } else {
    K1.v = forceapprox(old->a)*dt;
  }
  //K2=g(tn*dt/2,yn+k1/2)*dt
  new.a=old->a+0.5*K1.x;
  new.w=old->w+0.5*K1.v;
  K2.x = new.w*dt; // posizione
  if(f){
  K2.v = force(new.a)*dt;
  } else {
    K2.v = forceapprox(new.a)*dt; //velocità
  }  
//K3=g(tn+dt/2,yn+k2/2)*dt
  new.a=old->a+0.5*K2.x;
  new.w=old->w+0.5*K2.v;
  K3.x = new.w*dt; //posizione
  if(f){
    K3.v = force(new.a)*dt;
  } else {
  K3.v = forceapprox(new.a)*dt;//velocità
  }
  //K4=g(tn+dt, yn+k3)*dt
  new.a = old->a+K3.x;
  new.w = old->w+K3.v;
  K4.x = new.w*dt;//posizione
  if(f){
    K4.v = force(new.a)*dt;
  } else {
    K4.v = forceapprox(new.a)*dt;// velocità
  }
  //ora calcoliamo il passo successivo
  new.a = old->a + 1./6.*(K1.x+2*K2.x+2*K3.x+K4.x);
  new.w = old->w + 1./6.*(K1.v+2*K2.v+2*K3.v+K4.v);
  new.t = old->t+dt;
  return &new;
}
struct motion * findamax(struct motion* old, double dt, double* tempi, int* pcounter){
  static struct motion new;
  new = *verlet(old, dt, 1);
  if((old->a*new.a) < 0){
    *pcounter += 1;
  }
  if(*pcounter == 1 || *pcounter == 21){
    int i;
    struct point in, fin;
    if(*pcounter == 1){
      i = 0;
    } else {
      i = 1;
    }
    in.x = old->t;
    in.y = old->a;
    fin.x = new.t;
    fin.y = new.a;
    *(tempi+i) = fit(&in, &fin);
 }
 return &new;
}
double energy(struct motion*awt){
  return 0.5*k*k*awt->w*awt->w+1-cos(awt->a);
}
struct motion* cauchy(double amax){
  static struct motion cauchy;
  cauchy.a = amax;
  cauchy.w = w0;
  cauchy.t = 0.;
  return &cauchy;
}
double fit(struct point* p0, struct point* p1){
  double t;
  double m;
  m = (p1->y-p0->y)/(p1->x-p0->x);
  t = p0->x-p0->y/m;
  return t;
}