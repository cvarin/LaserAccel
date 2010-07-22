/******************************************************************************
  Description: Initialisation des paramètres de simulation.
*******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "constants.h"

extern double P,Imax,lambda,zf,wo,dT;
extern double ka,omega,z_rayleigh,Eo,A,T;
extern double phaseo;
extern double eps,h1,hmin,x1,x2;

void initialize(int n) /*Conditions initiales et paramètres d'intégration*/
{
 FILE *entree;

 switch(n){
    case 1:
        if((entree = fopen("./input/faisceau.arg", "r")) == NULL){
            printf("\nImpossible d\'ouvrir le fichier \'faisceau.arg\'\n\n");
            exit(1);}
        else{
            int result;
            result = fscanf(entree, "Puissance (Watts) : %lf", &P);
            result = fscanf(entree, "\nLongueur d'onde (m) : %lf", &lambda);
            result = fscanf(entree, "\n\nPosition du foyer (m) : %lf", &zf);
            result = fscanf(entree, "\n\n\nDimension du faisceau au foyer (m) : %lf", &wo);
            result = fscanf(entree, "\n\n\n\nLargeur de l'impulsion (fs) : %lf", &T);
            result = fscanf(entree, "\n\n\n\n\nPhase (x Pi rads) : %lf", &phaseo);
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
	    int result;
            result = fscanf(entree, "Precision (eps) : %lf", &eps);
            result = fscanf(entree, "\nPas de depart (h1) : %lf", &h1);
            result = fscanf(entree, "\n\nPas minimal permis (hmin) : %lf", &hmin);
            result = fscanf(entree, "\n\n\nTemps initial (x1, secondes) : %lf", &x1);
            result = fscanf(entree, "\n\n\n\nTemps final (x2, secondes) : %lf", &x2);
            fclose(entree);
            }
        break;
    
    default:
            printf("\nMauvais paramètres d\'initialisation...\n\n");
            exit(1);
    }
}
