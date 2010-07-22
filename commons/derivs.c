/******************************************************************************
  Équations du mouvement 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "constants.h"
#include "cust_math.h"

extern double dz,nz,tprime,impulsion,waist,gouy;
extern double phase,Ez,mag;
extern double zf,z_rayleigh,zpo,T,omega,A;
extern double k;
extern int nrhs;

complex double Ez_complex;
complex double Z,kR;
extern double z_conf;

/******************************************************************************/
void TM01_paraxial(double x,double y[],double dydx[])
{
    nrhs++;
    
    dz = y[1]-zf;
    nz = dz/z_rayleigh;
    tprime = x - dz/co + zpo/co;
    impulsion = exp(-(tprime*tprime)/(T*T));
    waist = 1/sqrt(1+(nz*nz));
    gouy = atan(dz/z_rayleigh);
    
    Ez = 2*A*(waist*waist)*impulsion*sin(omega*tprime+2*gouy-phase);
    
    mag = sqrt(1-((y[2]*y[2])/(co*co)));
    dydx[1] = y[2];
    dydx[2] = q/m*mag*(1-((y[2]*y[2])/(co*co)))*Ez;
}

/******************************************************************************/
void TM01(double x,double y[],double dydx[])
// Ref: A. April, Opt. Lett. 331563 (2008).
//      C. J. R. Sheppard and S. Saghafi, Opt. Lett 24, 1543 (1999).
{
     nrhs++;
    
     dz = y[1]-zf;
     nz = dz/z_rayleigh;
     tprime = x - dz/co + zpo/co;
     impulsion = exp(-(tprime*tprime)/(T*T));
     waist = 1/sqrt(1+(nz*nz));
     gouy = atan(dz/z_rayleigh);
     
     Z = dz + z_conf*I;
     kR = k*csqrt(Z*Z);
     
     // Manque l'amplitude du champ électrique
     // E0 = 1j*sqrt(P/2.0/(pi*a/(2.0*k))*(4.0*eta0))*exp(-k*a)
     Ez_complex = 2.0/3.0*(spherical_j0(kR) + spherical_j2(kR))
                   *cexp(-omega*x*I);
     Ez = creal(Ez_complex);
     
     mag = sqrt(1-((y[2]*y[2])/(co*co)));
     dydx[1] = y[2];
     dydx[2] = q/m*mag*(1-((y[2]*y[2])/(co*co)))*Ez;
}
/****************** End of file ***********************************************/