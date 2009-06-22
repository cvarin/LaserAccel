     /*
  Name: main.c
  Author: Charles Varin
  Date: 07-12-05 21:56
  Description: Simulation des équations du mouvement d'un électron sur l'axe de 
  propagation d'une impulsion TM01 (1-D).
  Tout est en unité MKS.
*/

#include <cstdlib>
#include <string>

// Les en-têtes du programme sont appellées ici
#include "MainHeader.h"

// Nombre d'équations à résoudre (une pour z + une pour vz = 2)
#define N 2              

/*************** Déclaration des variables globales ***************************/
double Wo,q,m,m_mev,vo,zini,zinter;
double P,Imax,lambda,zf,wo,dT,zpo;
double ka,omega,z_rayleigh,Eo,A,T;
double eps,h1,hmin,x1,x2;
double phaseo;
double dz,nz,tprime,impulsion,waist,gouy,coeff;
double phase,Ez,mag;

int kmax,kount;
double *xp,**yp,dxsav;
int nrhs;

const int sleep_time = 1;

const char *inputfile  = "./input/file1.arg";
const char *output1    = "./output/trajectoire.dat";
const char *output2    = "./output/balayage.dat";
const char *output_log = "./output/balayage.log";
const char *notefile   = "./output/notes.txt";

/******************************************************************************/
int main(void)
{
     FILE *sortie;
     FILE *log;
     
     int i,nbad,nok;
     double *ystart;
     
     int recouv;
     int option;
     double position;
     
     // Paramètres utilisés pour le balayage de la phase
     const double phasei = 0.0;
     const double phasef = 2.0;
     const int npt = 100;
     
     // Paramètre utilisé par la fonction odeint du Numerical Recipes
     kmax=200000;
     dxsav=(x2-x1)/1e2;
     
     // Paramètres de normalisation du temps et de la position de l'électron
     const double t_norm = 1.0e9; // 1e9 = Temps en nanosecondes
     const double z_norm = 1.0e3; // 1e3 = Position en millimètres
     // Note : l'énergie écrite dans les fichiers est toujours en MeV
     while(option!=3)
     {
          int status;
          /********* Initialisation de l'écran (vider l'écran) ************************/
          status = system("clear"); // Si ça ne fonctionne pas, essayer system("cls");
          assert(!status);
          /****************************************************************************/
          
          /*************** Options affichées à l'écran ********************************/
          std::cout << "\nOptions possibles:\n\n"
                    << "1. Trajectoire\n"
                    << "2. Balayage de la phase\n" 
                    << "3. Quitter\n\n"
                    << "Entrez votre choix et appuyez sur [ENTER] : ";
          std::cin  >> option;
          std::cout << std::endl;
          
          /******** Initialisation des paramètres de la simulation ********************/
          if(option!=3)
          {
               readfile(inputfile);
               
               std::cout << "Point de rencontre (en unités de zR)";
               std::cout << std::endl;
               std::cout << "(valeur précédente:" << position << "):";
               std::cin  >> position;
               std::cout << std::endl;
                    
               recouv = 8;
               zinter = position*z_rayleigh;               
               zini=-(recouv*vo*T/(1-vo/co)-zinter);   //Position initiale relative
               zpo=-(recouv*co*T/(1-vo/co)-zinter);    //Position initiale du centre de l'impulsion
                    
               // Temps de la simulation
               x1 = 0.0;
               x2 = 2.0*(-zpo)/co;
          }
     
          /****************************************************************************/               
          /*********** Selon l'option choisie, exécution ******************************/
          switch(option)
          {            
               case 1:
                    nrhs=0;
                              
                    // Allocation de la mémoire
                    ystart=dvector(1,N);
                    xp=dvector(1,kmax);
                    yp=dmatrix(1,N,1,kmax);
                              
                    // Conditions initiales de l'électron
                    ystart[1]=zf+zini; //Position initiale
                    ystart[2]=vo;      //Vitesse initiale
                              
                    // Phase de l'impulsion à l'étranglement
                    phase = phaseo*Pi; 
          
                    // Simulation
                    odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);
          
                    // Diagnostique à l'écran
                    printf("%s %13s %3d\n","successful steps:"," ",nok);
                    printf("%s %8s %3d\n","bad steps (corrected):"," ",nbad);
                    printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
                    printf("%s %3d\n\n","stored intermediate values:    ",kount);
                    
                    /*******************************************************************/          
                    // écriture du fichier des résultats
                    sortie = fopen(output1, "w");
                    // Premier point (à t = 0).
                    fprintf(sortie, "%.12f %.12f %.12f \n",
                                   x1*t_norm,          // Temps
                                   zini*z_norm,        // Position
                                   Wo-Wo);             // énergie
                    // Points subséquents
                    for (i=1;i<=kount;i++)
                    {
                         fprintf(sortie, "%.12f %.12f %.12f\n",
                                   xp[i]*t_norm,               
                              (yp[1][i]-zf)*z_norm,
                                   m_mev/sqrt(1-((yp[2][i]*yp[2][i])/(co*co)))-Wo);
                    }
                    fclose(sortie);
                    /*******************************************************************/          
                    // écriture du fichier de notes
                    notes(notefile);
          
                    // Désallocation de la mémoire
                    free_dmatrix(yp,1,N,1,kmax);
                    free_dvector(xp,1,kmax);
                    free_dvector(ystart,1,N);
                    sleep(sleep_time); 
               break;
     
          
               case 2:
                    nrhs=0;
                              
                    //Ouverture du fichier d'écriture des résultats
                    sortie = fopen(output2, "w");
                    log = fopen(output_log, "w");
                              
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
                         ystart[2]=vo;      //Vitesse initiale
               
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
                                        
                         //écriture du fichier de rsultats
                         fprintf(sortie, "%.12f %.12f \n", phase/Pi, 
                         m_mev/sqrt(1-((yp[2][kount]*yp[2][kount])/(co*co)))-Wo );
                                        
                         notes(notefile);
               
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
