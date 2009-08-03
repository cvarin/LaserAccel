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
double Wo,q,m,m_mev,vo,position,zini,zinter;
double P,Imax,lambda,zf,wo,dT,zpo;
double ka,omega,z_rayleigh,Eo,A,T;
double Pini,Pfin;
int    Np;
double eps,h1,hmin,x1,x2;
double phaseo;
double dz,nz,tprime,impulsion,waist,gouy,coeff;
double phase,Ez,mag;

int kmax,kount;
double *xp,**yp,dxsav;
int nrhs;

const int sleep_time = 1;

// Paramètres de normalisation du temps et de la position de l'électron
const double t_norm = 1.0e9; // 1e9 = Temps en nanosecondes
const double z_norm = 1.0e3; // 1e3 = Position en millimètres

const char maxchar = 90;

const char *outputfolder = "./output/";
const char *inputfile  = "./input/file1.arg";
const char *output1    = "./output/trajectoire.dat";
const char *notefile   = "./output/notes.txt";

/******************************************************************************/     
inline void print_dashes(void){std::cout << "------------------------------------------------------\n";};
inline void print_stars(void){std::cout <<  "******************************************************\n";};
void calculate_trajectory(void);
double *Get_Optimums(double *phase, double *energy, int count);
double *phase_scan(void);
void power_scan(void);

/******************************************************************************/
int main(void)
{    
     int recouv;
     int option;
     double *OPT = (double*)calloc(4,sizeof(double));
     
     // Paramètre utilisé par la fonction odeint du Numerical Recipes
     kmax=200000;
     dxsav=(x2-x1)/1e2;
     
     // Note : l'énergie écrite dans les fichiers est toujours en MeV
     while(option!=0)
     {    
          /*************** Options affichées à l'écran ********************************/
          std::cout << std::endl;
          print_stars();
          std::cout << "\nOptions possibles:\n\n"
                    << "1. Trajectoire\n"
                    << "2. Balayage de la phase\n"
                    << "3. Balayage de l\'intensité\n" 
                    << "0. Quitter\n\n"
                    << "Entrez votre choix et appuyez sur [ENTER] : ";
          std::cin  >> option;
          std::cout << std::endl;
          
          /******** Initialisation des paramètres de la simulation ********************/
          if(option!=0)
          {
               readfile(inputfile);
                    
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
               case 1: calculate_trajectory (); break;
               case 2: print_dashes();
                       OPT = phase_scan();
                       puts("Valeurs optimales du balayage de la phase");
                       printf("Max: W = %g MeV à la phase = %g pi-rads\n",OPT[1],OPT[0]/Pi);
                       printf("Min: W = %g MeV à la phase = %g pi-rads\n",OPT[3],OPT[2]/Pi);
               break;
               case 3: power_scan(); break;
               default: break;
          }// Fin du switch d'exécution                   
     }// Fin du while. Le programme y demeure tant que l'option 3 n'est pas choisie
     free(OPT);
     return (0);
}

/******************************************************************************/
/************************* Fonctions locales **********************************/
/******************************************************************************/
void calculate_trajectory(void)
{
     nrhs=0;
     int i,nbad,nok;
     double *ystart;
     FILE *sortie = fopen(output1, "w");

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
}

/******************************************************************************/
double *Get_Optimums(double *phase, double *energy, int count)

{
     int i;
     double *opt = (double*)calloc(4,sizeof(double));
     opt[1] = energy[0];
     opt[3] = energy[0];
     for(i=count;i--;)
     {
         if(energy[i] > opt[1]) opt[0] = phase[i],opt[1] = energy[i];
         if(energy[i] < opt[3]) opt[2] = phase[i],opt[3] = energy[i];
     }
     return opt;
     free(opt);
}

/******************************************************************************/
double *phase_scan(void)
{
     nrhs=0;
     double *ystart;
     int nbad,nok;
     
     char file[maxchar];
     sprintf(file,"%sbalayage_phase_%e.dat",outputfolder,P);
     
     clock_t t1,t2;
     t1 = clock();
     
     //Ouverture du fichier d'écriture des résultats
     FILE *sortie = fopen(file, "w");

     // Paramètres utilisés pour le balayage de la phase
     // NOTE IMPORTANTE : phase est une variable globale réservée!
     const double phasei = 0.0;
     const double phasef = 2.0;
     const int npt = 100;
     
     double *phase_vec = (double*)calloc(npt+1,sizeof(double));
     double *energy_vec = (double*)calloc(npt+1,sizeof(double));

     // Allocation de la mémoire
     ystart=dvector(1,N);
     xp=dvector(1,kmax);
     yp=dmatrix(1,N,1,kmax);

     //Boucle de balayage sur la phase
     printf("Puissance = %e W/cm^2, W0 = %g MeV\n",P,Wo);
     for(int st=0;st<=npt;st++)
     {
          // Sortie à l'écran de la progression du calcul
          std::cerr << "\r" << st << "%";
                         
          // Conditions initiales
          ystart[1]=zf+zini; //Position initiale
          ystart[2]=vo;      //Vitesse initiale

          // Phase de l'impulsion à l'étranglement
          phase = phase_vec[st] = (phasei + (phasef - phasei)*st/npt)*Pi; 

          // Simulation
          odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
                         
          // écriture du fichier de résultats
          energy_vec[st] = m_mev/sqrt(1-((yp[2][kount]*yp[2][kount])/(co*co)))-Wo;
          fprintf(sortie,"%.12f %.12f\n",phase_vec[st]/Pi,energy_vec[st]);
          notes(notefile); 
     } // Fin de la boucle de balayage

     free_dmatrix(yp,1,N,1,kmax);
     free_dvector(xp,1,kmax);
     free_dvector(ystart,1,N);

     //Fermeture des fichiers d'écriture       
     fclose(sortie);
     
     // Sortie à l'écran
     t2 = clock();
     clock_t tourshorloge = t2 - t1;
     std::cout << ", temps d'exécution : "
               <<  tourshorloge
               << " tours d\'horloge ("
               << (double)tourshorloge/CLOCKS_PER_SEC
               << " s)"
               << std::endl;
     return Get_Optimums(phase_vec,energy_vec,npt+1);
     free(phase_vec);
     free(energy_vec);
}

/******************************************************************************/
void power_scan(void)
{    
     double *OPT = (double*)calloc(4,sizeof(double));
     const double pas = (Pfin - Pini)/Np;
     const double WtoTW = 1.0e-12;
     
     /************* Fichier de sortie *****************************************/
     char scanfile[maxchar];
     sprintf(scanfile,"%sbalayage_puissance.dat",outputfolder);
     FILE *sortie = fopen(scanfile,"w");
     if(sortie==NULL) printf("N'a pu créer %s\n",scanfile),abort();
     
     for(int n = 0; n <= Np; n++)
     {
          P = Pini + n*pas;
          Imax = 2.0*P/(Pi*exp(1.0)*wo*wo);
          Eo = sqrt(2.0*120.0*Pi*Imax);
          A = 0.371*lambda/wo*Eo;
          
          print_dashes();
          printf("%d/%d\n",n,Np);
          OPT = phase_scan();
          puts("Valeurs optimales du balayage de la phase");
          printf("Max: W = %g MeV à la phase = %g pi-rads\n",OPT[1],OPT[0]/Pi);
          printf("Min: W = %g MeV à la phase = %g pi-rads\n",OPT[3],OPT[2]/Pi);
          
          //               1  2  3  4  5  6    1        2     3     4      5       6
          fprintf(sortie,"%e %e %e %e %e %e\n",P*WtoTW,Imax,OPT[1],OPT[0],OPT[3],OPT[2]);
     }
     
     fclose(sortie);
     free(OPT);
}

/****************** End of file ***********************************************/
