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


void readfile(const char *inputfile)
{
     FILE *entree = fopen(inputfile, "r");
     const int nparams = 16;
     int resultat[16];
      
     if(entree == NULL)
     {
          printf("\nImpossible d\'ouvrir le fichier \'%s\'\n\n",inputfile);
          exit(1);
     }
     else
     {    
          resultat[0] = !fscanf(entree, "# Paramètres du faisceau\n");
          resultat[1] = fscanf(entree, "Puissance (Watts) : %lf\n", &P);
          resultat[2] = fscanf(entree, "Longueur d'onde (m) : %lf\n", &lambda);
          resultat[3] = fscanf(entree, "Position du foyer (m) : %lf\n", &zf);
          resultat[4] = fscanf(entree, "Dimension du faisceau au foyer (m) : %lf\n", &wo);
          resultat[5] = fscanf(entree, "Largeur de l'impulsion (multiple de la periode) : %lf\n", &dT);
          resultat[6] = fscanf(entree, "Phase (x Pi rads) : %lf\n", &phaseo);
          
          resultat[7]  = !fscanf(entree, "# Paramètres de la particule\n");
          resultat[8]  = fscanf(entree, "Energie (MeV): %lf\n", &Wo);
          resultat[9]  = fscanf(entree, "Charge (coulombs) : %lf\n", &q);
          resultat[10] = fscanf(entree, "Masse (kg) : %lf\n", &m);
          resultat[11] = fscanf(entree, "Masse (MeV) : %lf\n", &m_mev);
          
          resultat[12] = !fscanf(entree, "# Paramètres de l\'intégrateur\n");
          resultat[13] = fscanf(entree, "Precision (eps) : %lf\n", &eps);
          resultat[14] = fscanf(entree, "Pas de depart (h1) : %lf\n", &h1);
          resultat[15] = fscanf(entree, "Pas minimal permis (hmin) : %lf\n", &hmin);
          
          for(int i=nparams;i--;) assert(resultat[i]);
          
          fclose(entree);
     } 
     // Paramètres du faisceau
     ka = 2*Pi/lambda;
     omega = ka*co;
     z_rayleigh = ka*(wo*wo)/2;
     Imax = 2*P/(Pi*exp(1)*wo*wo);
     Eo = sqrt(2*120*Pi*Imax);
     A = 0.371*lambda/wo*Eo;
     T = dT*2*Pi/omega;
     
     // Paramètre de la particule
     vo = co*sqrt(1-((m_mev*m_mev)/(Wo*Wo)));
     
     // Vérification d'usage
     if (Wo < 0.511)
     { 
          std::cout << "Attention! ";
          std::cout << "L'énergie initiale de l'électron doit être ";
          std::cout << "égale ou supérieure à l'énergie de masse. ";
          std::cout << "(0.511 = au repos)" << std::endl;
          std::cout << "Corriger et redémarrer.";
          sleep(1);
          exit(EXIT_FAILURE);
     }
     
     if (Wo == 0.511)
     {
          std::cout << "Le code n'a pas été développé pour ";
          std::cout << "traiter la cas d'un électron au repos.";
          std::cout << std::endl;
          std::cout << "Corriger et redémarrer." << std::endl;
          sleep(1);
          exit(EXIT_FAILURE);
     }
}


