/*
  Name: initialize.c
  Author: Charles Varin
  Date: 08-12-05 10:43
  Description: Initialisation des équations du mouvement 
  et des paramètres de simulation.
*/

#include "MainHeader.h"

extern double Wo,vo,zini;
extern double P,Imax,lambda,zf,wo,dT,zpo;
extern double ka,omega,z_rayleigh,Eo,A,T;
extern double eps,h1,hmin,x1,x2;
extern double phaseo;
extern double dz,nz,tprime,impulsion,waist,gouy,coeff;
extern double phase,Ez,mag;

extern int kmax,kount;
extern double *xp,**yp,dxsav;

extern int nrhs;


void derivs(double x,double y[],double dydx[])
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

void initialize(int n) /*Conditions initiales et paramètres d'intégration*/
{
 FILE *entree;

 switch(n){
    case 1:
        if((entree = fopen("./input/faisceau.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'faisceau.arg\'\n\n");
            exit(1);}
        else{
            fscanf(entree, "Puissance (Watts) : %lf", &P);
            fscanf(entree, "\nLongueur d'onde (m) : %lf", &lambda);
            fscanf(entree, "\n\nPosition du foyer (m) : %lf", &zf);
            fscanf(entree, "\n\n\nDimension du faisceau au foyer (m) : %lf", &wo);
            fscanf(entree, "\n\n\n\nLargeur de l'impulsion (fs) : %lf", &T);
            fscanf(entree, "\n\n\n\n\nPhase (x Pi rads) : %lf", &phaseo);
            fclose(entree);
            }
            ka = 2*Pi/lambda;
            omega = ka*co;
            z_rayleigh = ka*(wo*wo)/2;
            Imax = 2*P/(Pi*exp(1)*wo*wo);
            Eo = sqrt(2*120*Pi*Imax);
            A = 0.371*lambda/wo*Eo;
            T *= 1.0e-15;      // Conversion de femtoseconde à seconde
            dT = omega*T/2/Pi;
        break;

    case 2:
        if((entree = fopen("./input/integrateur.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'integrateur.arg\'\n\n");
            exit(1);}
        else{
            fscanf(entree, "Precision (eps) : %lf", &eps);
            fscanf(entree, "\nPas de depart (h1) : %lf", &h1);
            fscanf(entree, "\n\nPas minimal permis (hmin) : %lf", &hmin);
            fscanf(entree, "\n\n\nTemps initial (x1, secondes) : %lf", &x1);
            fscanf(entree, "\n\n\n\nTemps final (x2, secondes) : %lf", &x2);
            fclose(entree);
            }
        break;
    
    default:
            printf("\nMauvais paramètres d\'initialisation...\n\n");
            exit(1);
    }
}
