/******************************************************************************\
|** initialize_2D.c                                                          **|
|**                                                                          **|
|** Recueil de fonctions pour la simulation avec le faisceau                 **|
|** gaussien TM01 2-D pulsé.                                                 **|
|**                                                                          **|
|** derivs : Definition des equations du mouvement                           **|
|** (equation de Lorentz) sur l'axe d'un faisceau gaussien TM01 pulsé.       **|
|**                                                                          **|
|** initialize : Initialisation des paramètres d'integration.                **|
|**                                                                          **|
|** Charles Varin                                                            **|
|**                                                                          **|
|** 20 juillet 2005                                                          **|
\******************************************************************************/

#include <stdlib.h>
#include "bunch_accel.h"

extern double Wo,q,m,m_mev,vo;
extern double zpo,zini,zinter;
extern double ro,vro_norm,vro;
extern double Imax,lambda,zf,wo,dT;
extern double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
extern double eps,h1,hmin,x1,x2;
extern double phaseo,phasei,phasef,inti,intf;
extern unsigned int npt,npti;
extern double dz,nz,tprime,impulsion,waist,gouy;
extern double phase,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;
extern int nrhs;

void derivs(double x,double y[],double dydx[])
{
    nrhs++;
    
    dz = y[1]-zf;
    nz = dz/z_rayleigh;                               /*Position en z normalisée*/
    tprime = x - dz/co + zpo/co;      /*Temps retardé*/ 
    impulsion = exp(-(tprime*tprime)/(T*T));/*Impulsion gaussienne*/
    waist = 1/(1+(nz*nz));             /*Envelope de diffraction du faisceau gaussien*/
    gouy = atan(dz/z_rayleigh);                   /*Déphasage de Gouy*/
    space_env = exp(-y[3]*y[3]*waist/(wo*wo)); /*Enveloppe gaussienne transversale*/
    courbure = y[1]+z_rayleigh*z_rayleigh/y[1];  /*Paramètre de courbure du front d'onde*/
    Psi_c = ka*y[3]*y[3]/2/courbure;  /*Déphasage associé à la courbure du front d'onde*/
    Psi_0 = omega*tprime+2*gouy-Psi_c-phase; /*Phase de la porteuse*/
    
    Er = Ero*waist*impulsion*y[3]*space_env*cos(Psi_0);
    Ez = Ezo*waist*space_env*impulsion*(
         (1-y[3]*y[3]*waist/(wo*wo))*sin(Psi_0)-
          ka*y[3]*y[3]/2/courbure*cos(Psi_0));

    mag = sqrt(1-(y[2]*y[2])/(co*co)-(y[4]*y[4])/(co*co)); /*Inverse du facteur gamma*/
    dydx[1] = y[2];                                 /*Dérivée de z*/ 
    dydx[2] = q/m*mag*((1-(y[2]*y[2])/(co*co))*Ez + y[4]/co*(1-y[2]/co)*Er); /*Dérivée de la vitesse en z*/ 
    dydx[3] = y[4];                                 /*Dérivée de r*/ 
    dydx[4] = q/m*mag*((1-y[2]/co-(y[4]*y[4])/(co*co))*Er - y[2]*y[4]/(co*co)*Ez); /*Dérivée de la vitesse en r*/
}

void initialize(int n) /*Conditions initiales et paramètres d'intégration*/
{
 FILE *entree;

 switch(n){
    case 0:if((entree = fopen("particule.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'particule.arg\'\n\n");
            exit(1);}
        else{
            fscanf(entree, "Energie (MeV): %lf", &Wo);
            fscanf(entree, "\nCharge (coulombs) : %lf", &q);
            fscanf(entree, "\n\nMasse (kg) : %lf", &m);
            fscanf(entree, "\n\n\nMasse (MeV) : %lf", &m_mev);
            fclose(entree);
            }
            vo = co*sqrt(1-((m_mev*m_mev)/(Wo*Wo)));            
            vro = vro_norm*co;
        break;

    case 1:if((entree = fopen("faisceau.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'faisceau.arg\'\n\n");
            exit(1);}
        else{
            fscanf(entree, "Intensite maximale (W/m^2) : %lf", &Imax);
            fscanf(entree, "\nLongueur d'onde (m) : %lf", &lambda);
            fscanf(entree, "\n\nPosition du foyer (m) : %lf", &zf);
            fscanf(entree, "\n\n\nDimension du faisceau au foyer (m) : %lf", &wo);
            fscanf(entree, "\n\n\n\nLargeur de l'impulsion (multiple de la periode) : %lf", &dT);
            fscanf(entree, "\n\n\n\n\nPhase (x Pi rads) : %lf", &phaseo);
            fclose(entree);
            }
            ka = 2*Pi/lambda;
            omega = ka*co;
            z_rayleigh = ka*(wo*wo)/2;
            Eo = sqrt(2*120*Pi*Imax);
            Ero = Eo*exp(0.5)*sqrt(2)/wo;
            Ezo = Eo*exp(0.5)*2*sqrt(2)/ka/wo; //0.3710927*lambda/wo*Eo;
            T = dT*2*Pi/omega;
        break;

    case 2:if((entree = fopen("integrateur.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'integrateur.arg\'\n\n");
            exit(1);}
        else{
            fscanf(entree, "Precision (eps) : %lf", &eps);
            fscanf(entree, "\nPas de depart (h1) : %lf", &h1);
            fscanf(entree, "\n\nPas minimal permis (hmin) : %lf", &hmin);
            fclose(entree);
            }
        break;

    default:printf("\nMauvais paramètres d\'initialisation...\n\n");
        exit(1);
    }
}
