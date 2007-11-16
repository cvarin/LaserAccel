
#include <cstdlib>

#include "MainHeader.h"

extern int NPTS;
extern double rscale,zscale,norm;

extern double P,Imax,lambda,zf,wo,dT;
extern double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
extern double dz,nz,tprime,impulsion,waist,gouy;
extern double zpo,phase,phaseo,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;

extern double eps,h1,hmin,x1,x2;

extern int kmax,kount;
extern double *xp,**yp,dxsav;
extern int nrhs;

void notes(void){

    FILE *sortie;
        
    sortie = fopen("notes.txt", "w"); 
    
    fprintf(sortie, "Nombre d'électrons : %d\n",NPTS);
    fprintf(sortie, "Échelle du nuage en r : %f\n",rscale);
    fprintf(sortie, "Échelle du nuage en z : %f\n",zscale);
    
    fprintf(sortie, "\nParamètres du faisceau\n\n");
    fprintf(sortie, "Puissance : %g W\n", P);
    fprintf(sortie, "Intensité : %g W/cm^2\n", Imax/(100*100));
    fprintf(sortie, "Longueur d'onde : %g microns\n", lambda*1e6);
    fprintf(sortie, "Dimension du faisceau au foyer : %g microns\n", wo*1e6);
    fprintf(sortie, "Durée de l'impulsion (fs) : %g fs\n", T*1e15);
    fprintf(sortie, "Phase de l'impulsion au foyer : %g Pi rads\n", phaseo);            
    fprintf(sortie, "Position du foyer : %g m\n", zf);
    fprintf(sortie, "Fréquence angulaire : %g rad/fs\n", omega/1e15);    
    fprintf(sortie, "Largeur de l'impulsion (multiple de 2*Pi/omega) : %g\n", dT);    
    fprintf(sortie, "Largeur de l'impulsion (multiple de zR) : %g\n", T*co/z_rayleigh);
    fprintf(sortie, "Position initiale de l'impulsion (unités de zR) : %g\n", zpo/z_rayleigh);
    fprintf(sortie, "Distance de Rayleigh (zR) : %g microns\n", z_rayleigh*1e6);
    
    fprintf(sortie, "\n\nParamètres de l'integrateur.arg\n\n");
    fprintf(sortie, "Precision (eps) : %e\n", eps);
    fprintf(sortie, "Pas de depart (h1) : %e\n", h1);
    fprintf(sortie, "Pas minimal permis (hmin) : %e\n", hmin);
    fprintf(sortie, "Temps initial (x1, secondes) : %e\n", x1);
    fprintf(sortie, "Temps final (x2, secondes) : %e\n", x2);  
    
    fclose(sortie);  
}
