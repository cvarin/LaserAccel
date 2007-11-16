/******************************************************************************\
|** notes.c                                                                  **|
|**                                                                          **|
|** Écriture des fichiers des paramètres utilisés pour les simulations avec  **|
|** le faisceau gaussien TM01 1-D pulsé.                                     **|
|**                                                                          **|
|**                                                                          **|
|** Charles Varin                                                            **|
|**                                                                          **|
|** 28 novembre 2003                                                         **|
\******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "bunch_accel.h"
#include <time.h>

//char *asctime(const struct tm *timeptr);

extern double Wo,q,m,m_mev,vo;
extern double ro,vro_norm,vro;
extern double Imax,lambda,zf,wo,dT,pulse_delay;
extern double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
extern double eps,h1,hmin,x1,x2;
extern double phaseo,phasei,phasef,inti,intf;
extern unsigned int npt,npti;
extern double dz,nz,tprime,impulsion,waist,gouy;
extern double phase,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;

extern int NPTS;
extern float rscale,zscale;
extern int nrhs;
extern int kmax,kount;
extern double *xp,**yp,dxsav;

void notes(void){

    FILE *sortie;
        
    sortie = fopen("notes.txt", "w"); 
    fprintf(sortie, "Nombre de points : %d\n",NPTS);
    fprintf(sortie, "Échelle du faisceau en r : %f\n",rscale);
    fprintf(sortie, "Échelle du faisceau en z : %f\n",zscale);

    fprintf(sortie, "\nphase.arg\n");
    fprintf(sortie, "Phase de l'impulsion au foyer : %f Pi rads\n", phaseo);
    
    fprintf(sortie, "\n\nfaisceau.arg\n\n");
    fprintf(sortie, "Intensite maximale : %e W/cm^2\n", Imax/(100*100));
    fprintf(sortie, "Longueur d'onde : %f microns\n", lambda*1e6);
    fprintf(sortie, "Fréquence angulaire : %f rad/fs\n", omega/1e15);
    fprintf(sortie, "Position du foyer : %f m\n", zf);
    fprintf(sortie, "Dimension du faisceau au foyer : %f microns\n", wo*1e6);
    fprintf(sortie, "Largeur de l'impulsion (multiple de 2*Pi/omega) : %f\n", dT);
    fprintf(sortie, "Largeur de l'impulsion : %f fs\n", T*1e15);
    fprintf(sortie, "Delai de l'impulsion : %f*T\n", pulse_delay);
    fprintf(sortie, "Distance de Rayleigh : %f microns\n", z_rayleigh*1e6);
    
    fprintf(sortie, "\n\nparticule.arg\n\n");
    fprintf(sortie, "Energie cinétique initiale : %f MeV\n", Wo-m_mev);
    
    fprintf(sortie, "\n\nintegrateur.arg\n\n");
    fprintf(sortie, "Precision (eps) : %e\n", eps);
    fprintf(sortie, "Pas de depart (h1) : %e\n", h1);
    fprintf(sortie, "Pas minimal permis (hmin) : %e\n", hmin);
    fprintf(sortie, "Temps initial (x1, secondes) : %e\n", x1);
    fprintf(sortie, "Temps final (x2, secondes) : %e\n", x2);  
    
    fclose(sortie);  
}
