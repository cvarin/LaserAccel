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

void notes(){

    FILE *sortie;
    
    sortie = fopen("./notes.txt", "w");

    fprintf(sortie, "\nParamètres du faisceau\n\n");
    fprintf(sortie, "Puissance : %g W\n", P);
    fprintf(sortie, "Intensité : %g W/cm^2\n", Imax/(100*100));
    fprintf(sortie, "Longueur d'onde : %g microns\n", lambda*1e6);
    fprintf(sortie, "Dimension du faisceau au foyer : %g microns\n", wo*1e6);
    fprintf(sortie, "Durée de l'impulsion (fs) : %g fs\n", T*1e15);
    fprintf(sortie, "Phase de l'impulsion au foyer (pour run seulement) : %g Pi rads\n", phaseo);            
    fprintf(sortie, "\nPosition du foyer : %g m\n", zf);
    fprintf(sortie, "Fréquence angulaire : %g rad/fs\n", omega/1e15);    
    fprintf(sortie, "Largeur de l'impulsion (multiple de 2*Pi/omega) : %g\n", dT);    
    fprintf(sortie, "Largeur de l'impulsion (multiple de zR) : %g\n", T*co/z_rayleigh);
    fprintf(sortie, "Position initiale de l'impulsion (unités de zR) : %g\n", zpo/z_rayleigh);
    fprintf(sortie, "Distance de Rayleigh (zR) : %g microns\n", z_rayleigh*1e6);
    
    fprintf(sortie, "\n\nParamètres de la particule\n\n");
    fprintf(sortie, "Energie initiale : %g MeV\n", Wo);
    fprintf(sortie, "Vitesse initiale (v/c): %g\n", vo/co);
    fprintf(sortie, "Position initiale (unite de zR): %g\n", zini/z_rayleigh);
    fprintf(sortie, "Masse (MeV) : %g\n", m_mev);
    fprintf(sortie, "Charge (coulombs) : %g\n", q);
    fprintf(sortie, "Masse (kg) : %g\n", m);

    fprintf(sortie, "\nFACTEUR DE RECOUVREMENT INITIAL : %g\n", -(zpo-zini)/(T*co));
    fprintf(sortie, "LIMITE INFÉRIEURE POUR L'INJECTION (unités de zR) : %g\n", vo*T/(1-vo/co)/z_rayleigh);
         
    fprintf(sortie, "\n\nParamètres de l'intégrateur\n\n");
    fprintf(sortie, "Precision (eps) : %g\n", eps);
    fprintf(sortie, "Pas de depart (h1) : %g\n", h1);
    fprintf(sortie, "Pas minimal permis (hmin) : %g\n", hmin);
    fprintf(sortie, "Temps initial (x1, secondes) : %g\n", x1);
    fprintf(sortie, "Temps final (x2, secondes) : %g\n", x2);  
            
    fclose(sortie);  
}
