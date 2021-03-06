#ifndef _MainHeader_H_
#define _MainHeader_H_

#include <cstdlib>
#include <stdio.h>
#include <math.h>

#include "nrutil.h"

#define Pi           3.1415926535897932
#define co           2.99792458e8    // Vitesse de la lumi�re dans le vide [m/s]
#define q            -1.602176e-19   // Charge de l'�lectron [C]
#define m            9.109381e-031   // Masse de l'�lectron [kg]
#define m_mev        0.510998902     // Masse de l'�lectron [MeV]

/******************** D�finitions utilis�es par initialize.c ******************/
#define faisceau 1
#define integrateur 2
/******************************************************************************/

void odeint(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double,
	double [], double *, double *, void (*)(double, double [], double [])));

void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
    double eps, double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []));
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
    int nstep, double yout[], void (*derivs)(double, double[], double[]));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
    int nv);

double gasdev(long *idum);
double ran1(long *idum);

void derivs(double x,double y[],double dydx[]);
void initialize(int n);
void notes(void);

#endif /* _MainHeader_H_ */
