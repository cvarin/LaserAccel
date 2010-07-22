/******************************************************************************
  Équations du mouvement 
*******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "constants.h"

extern double dz,nz,tprime,impulsion,waist,gouy,coeff;
extern double phase,Ez,mag;
extern double zf,z_rayleigh,zpo,T,omega,A;
extern int nrhs;

void TM01_paraxial(double x,double y[],double dydx[])
{
    nrhs++;
    
    dz = y[1]-zf;
    nz = dz/z_rayleigh;
    tprime = x - dz/co + zpo/co;
    impulsion = exp(-(tprime*tprime)/(T*T));
    waist = 1/sqrt(1+(nz*nz));
    gouy = atan(dz/z_rayleigh);
    coeff = 1/(T*omega);
    
    Ez = 2*A*(waist*waist)*impulsion*sin(omega*tprime+2*gouy-phase);
    
    mag = sqrt(1-((y[2]*y[2])/(co*co)));
    dydx[1] = y[2];
    dydx[2] = q/m*mag*(1-((y[2]*y[2])/(co*co)))*Ez;
}

