/*
  Name: main.c
  Author: Charles Varin
  Date: 07-12-05 21:56
  Description: Simulation des �quations du mouvement d'un �lectron sur l'axe de 
  propagation d'une impulsion TM01 (1-D).
  Tout est en unit� MKS.
*/

#include <iostream>  
#include <string>
#include <cstdlib>

// Les en-t�tes du programme sont appell�es ici
#include "MainHeader.h"

// Nombre d'�quations � r�soudre (une pour z + une pour vz = 2)
#define N 2              

/*************** D�claration des variables globales ***************************/
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

int sleep_time = 1;
/******************************************************************************/

int main(void)
{
 FILE *sortie;
 FILE *log;
 
 int i,nbad,nok;
 double *ystart;
 
 int recouv;
 int option;
 float position;
 
 // Param�tres utilis�s pour le balayage de la phase
 float phasei = 0.0;
 float phasef = 2.0;
 unsigned int npt = 100;
 
 // Param�tre utilis� par la fonction odeint du Numerical Recipes
 kmax=200000;
 dxsav=(x2-x1)/1e2;
 
 // Param�tres de normalisation du temps et de la position de l'�lectron
 float t_norm = 1e9; // 1e9 = Temps en nanosecondes
 float z_norm = 1e3; // 1e3 = Position en millim�tres
 // Note : l'�nergie �crite dans les fichiers est toujours en MeV
 
 while(option!=3)
 {
  /********* Initialisation de l'�cran (vider l'�cran) ************************/
  system("clear"); // Si �a ne fonctionne pas, essayer system("cls");
  /****************************************************************************/
  
  /*************** Options affich�es � l'�cran ********************************/
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
      {
       /***********************************************************************/
       case 1:
       case 2:     
              // Initialization des param�tres de la simulations
              // � sp�cifier dans les fichiers se terminant par .arg
              initialize(particule);
              initialize(faisceau);      
              initialize(integrateur);
              
              if (Wo < 0.511)
                 { 
                  std::cout << "Attention! ";
                  std::cout << "L'�nergie initiale de l'�lectron doit �tre ";
                  std::cout << "�gale ou sup�rieure � l'�nergie de masse. ";
                  std::cout << "(0.511 = au repos)" << std::endl;
                  std::cout << "Corriger et red�marrer.";
                  sleep(3);
                  exit(EXIT_FAILURE);
                 }
              
              if (Wo == 0.511)
                 {
                  std::cout << "Le code n'a pas �t� d�velopp� pour ";
                  std::cout << "traiter la cas d'un �lectron au repos.";
                  std::cout << std::endl;
                  std::cout << "Corriger et red�marrer." << std::endl;
                  sleep(3);
                  exit(EXIT_FAILURE);
                 }
              
              std::cout << "Point de rencontre (en unit�s de zR)";
              std::cout << std::endl;
              std::cout << "(valeur pr�c�dente:" << position << "):";
              std::cin  >> position;
              std::cout << std::endl;
                   
              recouv = 8;
              zinter = position*z_rayleigh;               
              zini=-(recouv*vo*T/(1-vo/co)-zinter);   //Position initiale relative
              zpo=-(recouv*co*T/(1-vo/co)-zinter);    //Position initiale du centre de l'impulsion
                  
              // Temps de la simulation
              x1 = 0.0;
              x2 = 2.0*(-zpo)/co;
              
              break;
        /**********************************************************************/          
        default: 
              break;
      }
    }
  
  /****************************************************************************/               
  /*********** Selon l'option choisie, ex�cution ******************************/
  switch(option)
   {            
    case 1:
           nrhs=0;
                     
           // Allocation de la m�moire
           ystart=dvector(1,N);
           xp=dvector(1,kmax);
           yp=dmatrix(1,N,1,kmax);
                     
           // Conditions initiales de l'�lectron
           ystart[1]=zf+zini; //Position initiale
           ystart[2]=vo;      //Vitesse initiale
                     
           // Phase de l'impulsion � l'�tranglement
           phase = phaseo*Pi; 

           // Simulation
           odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);

           // Diagnostique � l'�cran
           printf("%s %13s %3d\n","successful steps:"," ",nok);
           printf("%s %8s %3d\n","bad steps (corrected):"," ",nbad);
           printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
           printf("%s %3d\n\n","stored intermediate values:    ",kount);
           
           /*******************************************************************/          
           // �criture du fichier des r�sultats
           sortie = fopen("./trajectoire.dat", "w");
           // Premier point (� t = 0).
           fprintf(sortie, "%.12f %.12f %.12f \n",
                           x1*t_norm,          // Temps
                           zini*z_norm,        // Position
                           Wo-Wo);             // �nergie
           // Points subs�quents
           for (i=1;i<=kount;i++)
               {
                fprintf(sortie, "%.12f %.12f %.12f\n",
                          xp[i]*t_norm,               
                         (yp[1][i]-zf)*z_norm,
                          m_mev/sqrt(1-((yp[2][i]*yp[2][i])/(co*co)))-Wo);
               }
           fclose(sortie);
           /*******************************************************************/          
           // �criture du fichier de notes
           notes();

           // D�sallocation de la m�moire
           free_dmatrix(yp,1,N,1,kmax);
           free_dvector(xp,1,kmax);
           free_dvector(ystart,1,N);
           sleep(sleep_time); 
           break;
    
         
    case 2:
           nrhs=0;
                     
           //Ouverture du fichier d'�criture des r�sultats
           sortie = fopen("./balayage.dat", "w");
           log = fopen("./balayage.log", "w");
                     
           //Boucle de balayage sur la phase
           for(int st=0;st<=npt;st++)
              {
               // Allocation de la m�moire
               ystart=dvector(1,N);
               xp=dvector(1,kmax);
               yp=dmatrix(1,N,1,kmax);
    
               // Sortie � l'�cran de la progression du calcul
               std::cerr << "\r" << st << "%";
                             
               // Conditions initiales
               ystart[1]=zf+zini; //Position initiale
               ystart[2]=vo;      //Vitesse initiale
    
               // Phase de l'impulsion � l'�tranglement
               phase = (phasei + (phasef - phasei)*st/npt)*Pi; 
    
               // Simulation
               odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
    
               //Sortie du journal d'ex�cution
               fprintf(log,"%d",st);
               fprintf(log,"\n%s %13s %f\n","Phase:","  ",phase/Pi);
               fprintf(log,"%s %13s %3d\n","successful steps:"," ",nok);
               fprintf(log,"%s %8s %3d\n","bad steps (corrected):","",nbad);
               fprintf(log,"%s %9s %3d\n","function evaluations:"," ",nrhs);
               fprintf(log,"%s %3d\n\n","stored intermediate values:    ",kount);
                             
               //�criture du fichier de r�sultats
               fprintf(sortie, "%.12f %.12f \n", phase/Pi, 
               m_mev/sqrt(1-((yp[2][kount]*yp[2][kount])/(co*co)))-Wo );
                             
               notes();

               free_dmatrix(yp,1,N,1,kmax);
               free_dvector(xp,1,kmax);
               free_dvector(ystart,1,N);
              } // Fin de la boucle de balayage
              
           //Fermeture des fichiers d'�criture       
           fclose(log);
           fclose(sortie);
                            
           std::cout << "\n\nBalayage termin�" << std::endl;
           sleep(sleep_time);
           break;
                            
    default: 
           break;
   }// Fin du switch d'ex�cution                   
 }// Fin du while. Le programme y demeure tant que l'option 3 n'est pas choisie
 return (0);
}
/************************* Fin du main ****************************************/
