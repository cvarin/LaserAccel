/*
  Name: main.c (bunch_accel_stage)
  Author: Charles Varin
  Date: 21-07-05 13:27
  Description: Simulation des équations du mouvement en 2D d'électrons relativistes
  au centre d'une impulsion TM01. Cette routine utilise l'algorithme Runge-Kutta d'ordre 5 
  avec pas adaptatif (Cash-Karp) et contrôle d'erreur du Numerical recipes en C 
  (chapitre 16). Tout est en unité MKS.
*/

#include <iostream>  // pour maintenir 
#include <string>    // l'invite de commande 
using namespace std; // ouverte après exécution.

#include <math.h>
#include <stdio.h>

#include "bunch_accel.h"  /*En-tête des routines d'intégration*/
#include "nrutil.h"  /*En-têtes des routines du NR*/

#define N 4              /*Nombre d'equations a resoudre*/

int NPTS;
float rscale,zscale;
int kmax,kount;        // Paramètres utilisés 
double *xp,**yp,dxsav; // par odeint_double.c
int nrhs;

double Wo,q,m,m_mev,vo;                      

double zpo,zini,zinter;

double ro,vro_norm,vro;                     
double Imax,lambda,zf,wo,dT,pulse_delay;
double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
double eps,h1,hmin,x1,x2;
double phaseo,phasei,phasef,inti,intf;
unsigned int npt,npti;
double dz,nz,tprime,impulsion,waist,gouy;
double phase,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;

double roo,zoo,gam,angle_degres;
double zfinal_onaxis,rfinal_onaxis,angle_onaxis,Wfinal_onaxis;
double position_finale_enveloppe,position_finale_porteuse;

int sleep_time = 1;

main()
{

int i,nbad,nok;
double *ystart;

int recouv;
int option;
float position;

long idum=(-1);

FILE *fichier_onaxis;
FILE *log;
FILE *fichier1;
FILE *fichier2;

// Attention, les paramètres de simulations sont à 3 microns!!!!!!
// 15 ps pour wo = 3 microns
// 40 ps pour wo = 10 microns
NPTS = 2000;
rscale = 0.0125; //0.1, d=400 nm ; 0.025, d=100 nm ; 0.0125, d=50 nm
zscale = 0.0125;

kmax=200000;
dxsav=(x2-x1)/1e2;

while(option!=2)
 {
  cout << "\33[2J"; //Réinitialisation de l'affichage
  cout << "\nOptions possibles:\n"
       <<"1. Simuler\n"
       <<"2. Quitter\n"
       <<"Choix (Appuyer sur ENTER pour démarrer):";     
  cin  >> option;
  cout << endl;

  switch(option)
       {            
          case 1:    
               initialize(particule);       /*Initialisation des paramètres de la particule*/
               initialize(faisceau);        /*Initialisation des paramètres du faisceau*/
               initialize(integrateur);     /*Initialisation des paramètres de l'intégrateur*/
               
               cout << "Point de rencontre (en unités de zR)";
               cout << "\n(valeur précédente:" << position << "):";
               cin >> position;
               cout << endl;
               
               recouv = 8;
               zinter = position*z_rayleigh;               
               zini=-(recouv*vo*T/(1-vo/co)-zinter);   //Position initiale relative
               zpo=-(recouv*co*T/(1-vo/co)-zinter);    //Position initiale du centre de l'impulsion
               
		       //Temps de la simulation
               x1 = 0;
               x2 = 2*(-zpo)/co;

               /*Calcul du mouvement de la particule de référence accélérée sur l'axe*/
               nrhs=0;
    
               ystart=dvector(1,N);
               xp=dvector(1,kmax);
               yp=dmatrix(1,N,1,kmax);
    
               phase = phaseo*Pi;
               ystart[1]=zf+zini;             /*Condition initiale sur la position en z*/
               ystart[2]=vo;             /*Condition initiale sur la vitesse en z*/
               ystart[3]=0;              /*Condition initiale sur la position en r*/
               ystart[4]=0;              /*Condition initiale sur la vitesse en r*/
    
               odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
               fprintf(log,"\n%s\n","Particule de référence (sur l'axe)");
               fprintf(log,"\n%s %13s %3d\n","successful steps:"," ",nok);
               fprintf(log,"%s %20s %3d\n","bad steps (corrected):"," ",nbad);
               fprintf(log,"%s %9s %3d\n","function evaluations:"," ",nrhs);
               fprintf(log,"%s %3d\n\n","stored intermediate values:    ",kount);
    
               gam = 1/sqrt(1-yp[2][kount]*yp[2][kount]/(co*co)-yp[4][kount]*yp[4][kount]/(co*co));
               zfinal_onaxis = yp[1][kount];
               rfinal_onaxis = yp[3][kount];
               angle_onaxis = atan(yp[4][kount]/yp[2][kount]);
               Wfinal_onaxis = gam*m_mev;
               position_finale_enveloppe = exp(-((xp[kount] - (yp[1][kount]-zf)/co - pulse_delay*T)*(xp[kount] - (yp[1][kount]-zf)/co - pulse_delay*T))/(T*T));
               position_finale_porteuse = asin(sin((2*Pi/lambda)*co*xp[kount]-(2*Pi/lambda)*(yp[1][kount]-zf)+2*atan((yp[1][kount]-zf)/z_rayleigh)-phase))/Pi;
                        
               //Écriture du fichier des résultats
               fichier_onaxis = fopen ("on_axis.dat","w");
               fprintf(fichier_onaxis, "%.12f \t %.12f \t %.12f \t %.12f \t %.12f \t %.12f \n",
               zfinal_onaxis,rfinal_onaxis*1e6,angle_onaxis,Wfinal_onaxis-Wo,position_finale_enveloppe,position_finale_porteuse);
               fclose(fichier_onaxis);
                                             
               free_dmatrix(yp,1,N,1,kmax);
               free_dvector(xp,1,kmax);
               free_dvector(ystart,1,N);
               
               /* Calcul du mouvement du gaz d'électrons */
                 
               nrhs=0; 
                  
               log = fopen("./log.txt", "w");
               fichier1 = fopen ("initial.dat","w");
               fichier2 = fopen ("final.dat","w");

               for (int l=1;l<=NPTS;l++) 
                   {
                       ystart=dvector(1,N);
                       xp=dvector(1,kmax);
                       yp=dmatrix(1,N,1,kmax);
 
                       //Sortie à l'écran de la progression du calcul
                       cerr << "\r" << l << "/" << NPTS;
                       
                       roo = rscale*1.0e-6*gasdev(&idum);
                       zoo = zscale*1.0e-6*gasdev(&idum)+zf+zini;
 
                       phase = phaseo*Pi;
                       ystart[1]=zoo;            /*Condition initiale sur la position en z*/
                       ystart[2]=vo;               /*Condition initiale sur la vitesse en z*/
                       ystart[3]=roo;              /*Condition initiale sur la position en r*/
                       ystart[4]=vro;              /*Condition initiale sur la vitesse en r*/

                       odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
 
                       //Sortie du journal d'exécution
                       fprintf(log,"\n%s %3d\n","Point #",l);
                       fprintf(log,"\n%s %13s %3d\n","successful steps:"," ",nok);
                       fprintf(log,"%s %8s %3d\n","bad steps (corrected):"," ",nbad);
                       fprintf(log,"%s %9s %3d\n","function evaluations:"," ",nrhs);
                       fprintf(log,"%s %3d\n\n","stored intermediate values:    ",kount);

                       gam = 1/sqrt(1-yp[2][kount]*yp[2][kount]/(co*co)-yp[4][kount]*yp[4][kount]/(co*co));
                       angle_degres = atan(yp[4][kount]/yp[2][kount])*360/2/Pi;

                       fprintf(fichier1, "%.12f \t %.12f \t %.12f \n",(zoo-zf)*1e6,roo*1e6,Wo-m_mev);
                       fprintf(fichier2, "%.12f \t %.12f \t %.12f \t %.12f\n",(yp[1][kount]-zfinal_onaxis)*1e6,yp[3][kount]*1e6,angle_degres,(gam*m_mev)-Wo);
                          
                       notes();
                          
                       free_dmatrix(yp,1,N,1,kmax);
                       free_dvector(xp,1,kmax);
                       free_dvector(ystart,1,N);
                   }//Fin du for
 
               fclose(log);
               fclose (fichier2);
               fclose (fichier1);

               cout << "\n\nCalcul terminé" << endl;
               sleep(sleep_time);
                     
               break;//Fin du cas 1
                         
          default: 
               break; //Fin de "default"
       //Fin du switch
       }
         
  //Fin du while
  } 
  
  return (0);
//Fin du main
}
