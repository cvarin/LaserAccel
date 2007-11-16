#ifndef _MainHeader_H_
#define _MainHeader_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "nrutil.h"  /*En-têtes des routines du NR*/

#define co 2.99792458e8 /* Vitesse de la lumiere dans le vide (Codata)*/
#define Pi 3.1415926535897932

#define particule 0  /*Définitions utilisées par initialize.c*/
#define faisceau 1
#define integrateur 2
#define phase_initiale 3

#define run 0  /*Définitions utilisées par notes.c*/
#define scan 1

void odeint(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double,
	double [], double *, double *, void (*)(double, double [], double [])));

//Pour Runge Kutta (Cash-Karp)
void rkck(double y[], double dydx[], int n, double x, double h,
	double yout[], double yerr[], void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x,
	double htry, double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
	
//Pour Bulirsh-Stoer
void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
    double eps, double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []));
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
    int nstep, double yout[], void (*derivs)(double, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
    int nv);

void derivs(double x,double y[],double dydx[]);
void initialize(int n);
void notes();

#endif /* _MainHeader_H_ */
