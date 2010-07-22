/*
  Name: main.c
  Author: Charles Varin
  Date: 07-12-05 21:56
  Description: Simulation des équations du mouvement d'un électron sur l'axe de 
  propagation d'une impulsion TM01 (1-D).
  Tout est en unité MKS.
*/

#include <iostream>  
#include <string>
#include <cstdlib>

// Les en-têtes du programme sont appellées ici
#include "MainHeader.h"

// Nombre d'équations à résoudre (une pour z + une pour vz = 2)
#define N 2              

/*************** Déclaration des variables globales ***************************/
double Wo,vo,zini;
double P,Imax,lambda,zf,wo,dT,zpo;
double ka,omega,z_rayleigh,Eo,A,T;
double eps,h1,hmin,x1,x2;
double phaseo;
double dz,nz,tprime,impulsion,waist,gouy,coeff;
double phase,Ez,mag;

int kmax,kount;
double *xp,**yp,dxsav;
int nrhs;

int sleep_time = 1;
/******************************************************************************/

int main(void)
{
 FILE *sortie;
 FILE *log;
 
 int i,nbad,nok;
 double *ystart;
 
 int option;
 float position;
 
 // Paramètres utilisés pour le balayage de la phase
 float phasei = 0.0;
 float phasef = 2.0;
 unsigned int npt = 100;
 
 // Paramètre utilisé par la fonction odeint du Numerical Recipes
 kmax=200000;
 dxsav=(x2-x1)/1e2;
 
 // Paramètres de normalisation du temps et de la position de l'électron
 float t_norm = 1e12; // 1e12 = Temps en picosecondes
 float z_norm = 1e3;  // 1e3 = Position en millimètres
 // Note : l'énergie écrite dans les fichiers est toujours en MeV
 
 while(option!=3)
 {
  /********* Initialisation de l'écran (vider l'écran) ************************/
  system("clear"); // Si ça ne fonctionne pas, essayer system("cls");
  /****************************************************************************/
  
  /*************** Options affichées à l'écran ********************************/
  std::cout << "\nOptions possibles:\n\n"
            << "1. Trajectoire\n"
            << "2. Balayage de la phase\n" 
            << "3. Quitter\n\n"
            << "Entrez votre choix et appuyez sur [ENTER] : ";
  std::cin  >> option;
  std::cout << std::endl;
  
  /****************************************************************************/
  /** Si l'option 3 (Quitter) n'est pas choisie, initialiser les variables ****/
  /****************************************************************************/
  if(option!=3)
    {
     switch(option)
      {case 1: case 2:     
       /***********************************************************************/
       /* Les paramètres de la simulations sont spécifiés dans les fichiers   */
       /* se terminant par l'extension .arg                                   */
       /***********************************************************************/
          initialize(faisceau);      
          initialize(integrateur);
    
          
          std::cout << "Position initiale de l'électron (en unités de zR)";
          std::cout << std::endl;
          std::cout << "(valeur précédente:" << position << "):";
          std::cin  >> position;
          std::cout << std::endl;
          
          zini=position*z_rayleigh;  // Position initiale de l'électron [m]
	  zpo=-0.0*co*T;               // Position initiale de l'impulsion [m] 	  //zpo=zf;        
          break;
        /**********************************************************************/          
        default: 
          break;
      }
    }
  
  /****************************************************************************/               
  /*********** Selon l'option choisie, exécution ******************************/
  /****************************************************************************/
  switch(option)
   {                       
    case 1:// Détail de la simulation pour une trajectoire
           nrhs=0;
                     
           // Allocation de la mémoire
           ystart=dvector(1,N);
           xp=dvector(1,kmax);
           yp=dmatrix(1,N,1,kmax);
                     
           // Conditions initiales de l'électron
           ystart[1]=zf+zini;//Position initiale
           ystart[2]=0.0;    //Vitesse initiale
           Wo = m_mev/sqrt(1-pow(ystart[2]/co,2));
           
           phase = phaseo*Pi; // Phase de l'impulsion à l'étranglement
          
           // Simulation
           odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);

           // Diagnostique à l'écran
           printf("%s %13s %3d\n","successful steps:"," ",nok);
           printf("%s %8s %3d\n","bad steps (corrected):"," ",nbad);
           printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
           printf("%s %3d\n\n","stored intermediate values:    ",kount);
           
           /*******************************************************************/          
           /* Écriture du fichier des résultats                               */
           /*******************************************************************/             
           sortie = fopen("./output/trajectoire.dat", "w");
           // Premier point (à t = 0).
           fprintf(sortie, "%.12f %.12f %.12f \n",
                           x1*t_norm,        // Temps
                           zini*z_norm,      // Position
                           Wo-m_mev          // Énergie
                   );    
           // Points subséquents
           for (i=1;i<=kount;i++)
               {
                fprintf(sortie, "%.12f %.12f %.12f\n",
                          xp[i]*t_norm,               
                         (yp[1][i]-zf)*z_norm,
                          m_mev/sqrt( 1-pow(yp[2][i]/co,2)) - Wo 
                        );
               }
           fclose(sortie);
           
           // Écriture du fichier de notes
           notes();

           // Désallocation de la mémoire
           free_dmatrix(yp,1,N,1,kmax);
           free_dvector(xp,1,kmax);
           free_dvector(ystart,1,N);
           sleep(sleep_time); 
           break;
    
         
    case 2:// Détail de la simulation pour un balayage de la phase
           nrhs=0;
                     
           //Ouverture du fichier d'écriture des résultats
           sortie = fopen("./output/balayage.dat", "w");
           log = fopen("./output/balayage.log", "w");
                     
           //Boucle de balayage sur la phase
           for(int st=0;st<=npt;st++)
              {
               // Allocation de la mémoire
               ystart=dvector(1,N);
               xp=dvector(1,kmax);
               yp=dmatrix(1,N,1,kmax);
    
               // Sortie à l'écran de la progression du calcul
               std::cerr << "\r" << st << "%";
                             
               // Conditions initiales
               ystart[1]=zf+zini; //Position initiale
               ystart[2]=0.0;      //Vitesse initiale
               
               Wo = m_mev/sqrt(1-pow(ystart[2]/co,2));
    
               // Phase de l'impulsion à l'étranglement
               phase = (phasei + (phasef - phasei)*st/npt)*Pi; 
    
               // Simulation
               odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
    
               //Sortie du journal d'exécution
               fprintf(log,"%d",st);
               fprintf(log,"\n%s %13s %f\n","Phase:","  ",phase/Pi);
               fprintf(log,"%s %13s %3d\n","successful steps:"," ",nok);
               fprintf(log,"%s %8s %3d\n","bad steps (corrected):","",nbad);
               fprintf(log,"%s %9s %3d\n","function evaluations:"," ",nrhs);
               fprintf(log,"%s %3d\n\n","stored intermediate values:    ",kount);
                             
               //Écriture du fichier de résultats
               fprintf(sortie, "%.12f %.12f \n", phase/Pi, 
               m_mev/sqrt( 1-pow(yp[2][kount]/co,2)) - Wo );
                             
               notes();

               free_dmatrix(yp,1,N,1,kmax);
               free_dvector(xp,1,kmax);
               free_dvector(ystart,1,N);
              } // Fin de la boucle de balayage
              
           //Fermeture des fichiers d'écriture       
           fclose(log);
           fclose(sortie);
                            
           std::cout << "\n\nBalayage terminé" << std::endl;
           sleep(sleep_time);
           break;
                            
    default: 
           break;
   }// Fin du switch d'exécution                   
 }// Fin du while. Le programme y demeure tant que l'option 3 n'est pas choisie
 return (0);
}
/************************* Fin du main ****************************************/
