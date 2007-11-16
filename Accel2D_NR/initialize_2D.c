/*
  Name: initialize_2D.c
  Author: Charles Varin
  Description: Initialisation des équations du mouvement 
  et des paramètres de simulation.
*/

#include <cstdlib>

#include "MainHeader.h"

extern int NPTS;
extern double rscale,zscale,norm;

extern double P,Imax,lambda,zf,wo,dT;
extern double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
extern double dz,nz,tprime,impulsion,waist,gouy;
extern double zpo,phase,phaseo,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;

extern double eps,h1,hmin,x1,x2;

extern int nrhs;

/********************* Équations du mouvement *********************************/
void derivs(double x,double y[],double dydx[])
{
  nrhs++;
    
  dz = y[1]-zf;
  nz = dz/z_rayleigh;                                     /*Position en z normalisée*/
  tprime = x - dz/co +zpo/co;                    /*Temps retardé*/ 
  impulsion = exp(-(tprime*tprime)/(T*T));        /*Impulsion gaussienne*/
  waist = 1/(1+(nz*nz));      /*Envelope de diffraction du faisceau gaussien*/
  gouy = atan(dz/z_rayleigh);                             /*Déphasage de Gouy*/
  space_env = exp(-y[3]*y[3]*waist/(wo*wo));/*Enveloppe gaussienne transversale*/
  courbure = dz+z_rayleigh*z_rayleigh/(dz+1.0e-30);   /*Paramètre de courbure du front d'onde*/
//    courbure = y[1]+z_rayleigh*z_rayleigh/y[1];   /*Paramètre de courbure du front d'onde*/
  Psi_c = ka*y[3]*y[3]/2/courbure; /*Déphasage associé à la courbure du front d'onde*/
  Psi_0 = omega*tprime+2*gouy-Psi_c-phase; /*Phase de la porteuse*/
    
  Er = Ero*waist*impulsion*y[3]*space_env*cos(Psi_0);
  Ez = Ezo*waist*space_env*impulsion*(
         (1-y[3]*y[3]*waist/(wo*wo))*sin(Psi_0)-
          ka*y[3]*y[3]/2/courbure*cos(Psi_0));
    
  mag = sqrt(1-(y[2]*y[2])/(co*co)-(y[4]*y[4])/(co*co));
    
  /*************** Dérivée de z ***********************************************/ 
  dydx[1] = y[2];
  /*************** Dérivée de la vitesse en z *********************************/ 
  dydx[2] = q/m*mag*((1-(y[2]*y[2])/(co*co))*Ez + y[4]/co*(1-y[2]/co)*Er);
  /* ************* Dérivée de r ***********************************************/ 
  dydx[3] = y[4];
  /*************** Dérivée de la vitesse en r *********************************/                
  dydx[4] = q/m*mag*((1-y[2]/co-(y[4]*y[4])/(co*co))*Er - y[2]*y[4]/(co*co)*Ez); 
}

/******** Lecture des conditions initiales et des paramètres d'intégration ****/
void initialize(int n) 
{
 FILE *entree;

 switch(n){

    case 1:
        if((entree = fopen("faisceau.arg", "r")) == NULL){
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
            Ero = Eo*exp(0.5)*sqrt(2)/wo;
            Ezo = Eo*exp(0.5)*2*sqrt(2)/ka/wo;
            T *= 1.0e-15;      // Conversion de femtoseconde à seconde
            dT = omega*T/2/Pi;
        break;

    case 2:if((entree = fopen("integrateur.arg", "r")) == NULL){
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
        printf("Options demandée non-disponible dans 'initialize.c'.\n");
        printf("\nAppuyer sur [ENTER] pour sortir du programme.\n");
        getchar();
        exit(1);
    }
}
/******************************************************************************/
