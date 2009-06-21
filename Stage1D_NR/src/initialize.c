/*
  Name: initialize.c
  Author: Charles Varin
  Date: 08-12-05 10:43
  Description: Initialisation des équations du mouvement 
  et des paramètres de simulation.
*/

#include "MainHeader.h"

extern double Wo,q,m,m_mev,vo,zini,zinter;
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
    double Ez1,Ez2;
    double waist2;
    double tprime2;
    double impulsion2;
    const double fac = 2.0;
    const double d = 0.0;
    
    nrhs++;
    
    dz = y[1]-zf;
    nz = dz/z_rayleigh;
    gouy = atan(dz/z_rayleigh);
    coeff = 1/(T*omega);
    
    // Pulse 1
    tprime = x - dz/co + zpo/co;
    impulsion = exp(-(tprime*tprime)/(T*T));
    waist = 1/sqrt(1+(nz*nz));
    
    // Pulse 2
    tprime2 = x - dz/co + zpo/co + d*zpo/co;
    impulsion2 = exp(-(tprime2*tprime2)/(1000.0*T*T));
    waist2 = 1/sqrt(1+(nz*nz)/(fac*fac));
    
    Ez1 = 2.0*A*(waist*waist)*impulsion*sin(omega*tprime+2*gouy-phase);
//     Ez2 = 2.0*A*(waist2*waist2)*impulsion2*sin(omega*tprime2+2*gouy+phase);
     Ez2 = 0.0;
    Ez = (Ez1+Ez2)/2.0;
    
    mag = sqrt(1-((y[2]*y[2])/(co*co)));
    dydx[1] = y[2];
    dydx[2] = q/m*mag*(1-((y[2]*y[2])/(co*co)))*Ez;
    
}

void initialize(int n) /*Conditions initiales et paramètres d'intégration*/
{
     FILE *entree;
     int resultat;
     switch(n)
     {
          case 0:
               if((entree = fopen("./input/particule.arg", "r")) == NULL)
               {
                    printf("\nImpossible d\'ouvrir le fichier \'particule.arg\'\n\n");
                    exit(1);
               }
               else
               {
                    resultat = fscanf(entree, "Energie (MeV): %lf\n", &Wo);
                    resultat = fscanf(entree, "Charge (coulombs) : %lf\n", &q);
                    resultat = fscanf(entree, "Masse (kg) : %lf\n", &m);
                    resultat = fscanf(entree, "Masse (MeV) : %lf\n", &m_mev);
                    fclose(entree);
               } 
               vo = co*sqrt(1-((m_mev*m_mev)/(Wo*Wo)));
          break;

          case 1:
               if((entree = fopen("./input/faisceau.arg", "r")) == NULL)
               {
                    printf("\nImpossible d\'ouvrir le fichier \'faisceau.arg\'\n\n");
                    exit(1);
               }
               else
               {
                    resultat = fscanf(entree, "Puissance (Watts) : %lf\n", &P);
                    resultat = fscanf(entree, "Longueur d'onde (m) : %lf\n", &lambda);
                    resultat = fscanf(entree, "Position du foyer (m) : %lf\n", &zf);
                    resultat = fscanf(entree, "Dimension du faisceau au foyer (m) : %lf\n", &wo);
                    resultat = fscanf(entree, "Largeur de l'impulsion (multiple de la periode) : %lf\n", &dT);
                    resultat = fscanf(entree, "Phase (x Pi rads) : %lf", &phaseo);
                    fclose(entree);
               }
               ka = 2*Pi/lambda;
               omega = ka*co;
               z_rayleigh = ka*(wo*wo)/2;
               Imax = 2*P/(Pi*exp(1)*wo*wo);
               Eo = sqrt(2*120*Pi*Imax);
               A = 0.371*lambda/wo*Eo;
               T = dT*2*Pi/omega;
          break;

          case 2:
               if((entree = fopen("./input/integrateur.arg", "r")) == NULL)
               {
                    printf("\nImpossible d\'ouvrir le fichier \'integrateur.arg\'\n\n");
                    exit(1);
               }
               else
               {
                    resultat = fscanf(entree, "Precision (eps) : %lf\n", &eps);
                    resultat = fscanf(entree, "Pas de depart (h1) : %lf\n", &h1);
                    resultat = fscanf(entree, "Pas minimal permis (hmin) : %lf", &hmin);
                    fclose(entree);
               }
          break;
    
          default:
               printf("\nMauvais paramètres d\'initialisation...\n\n");
               exit(1);
          break;
    }
}
