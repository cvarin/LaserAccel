/*
  Nom: main.c
  Auteur: Charles Varin
  Description: Simulation de l'acc�l�ration d'un nuage d'�lectrons par une
  une impulsion TM01. 
  Tout est en unit� MKS.
*/

// Les en-t�tes du programme sont appell�es ici
#include "MainHeader.h" 

// Nombre d'�quations � r�soudre (deux pour position + deux pour vitesses)
#define N 4

/*************** D�claration des variables globales ***************************/
int    NPTS;
double rscale,zscale,norm;

double gam_onaxis,zfinal_onaxis,Wfinal_onaxis;
double roo,zoo,gam,angle_degres;
                                      
double P,Imax,lambda,zf,wo,dT;
double ka,omega,z_rayleigh,Eo,Ezo,Ero,T;
double dz,nz,tprime,impulsion,waist,gouy;
double zpo,phase,phaseo,Ez,Er,mag,Psi_0,space_env,courbure,Psi_c;

double eps,h1,hmin,x1,x2;

int kmax,kount;        
double *xp,**yp,dxsav; 
int nrhs;
/******************************************************************************/

int main(void)
{
    /*************** Nombre d'�lectrons et dimension initiale du nuage ********/
    NPTS = 100;
    rscale = 0.0125; //0.1, d=400 nm ; 0.025, d=100 nm ; 0.0125, d=50 nm
    zscale = 0.0125;
    
    /***** Normalisation pour la position dans les fichiers *******************/
    norm = 1e6;  // 1e6 = Position en microns dans les fichiers
    // Mais Attention! zfinal_onaxis est en m�tres (zf + ... )
    // Note : l'�nergie �crite dans les fichiers est toujours en MeV
    
    /************** D�clarations des variables du 'main' **********************/
    int i,nbad,nok;
    double *ystart;
    
    long idum=(-1);
    
    FILE *fichier_onaxis;
    FILE *fichier1;
    FILE *fichier2;
    
    /*************** Param�tres pour 'odeint' *********************************/
    nrhs=0;
    kmax=200000;
    dxsav=(x2-x1)/1;
    
    /******* Initialisation (lecture des param�tres dans les fichiers) ********/
    initialize(faisceau);
    initialize(integrateur);
    zpo=-4*co*T;               // Position initiale de l'impulsion [m] 
    
    /******* Calcul du mouvement d'un �lectron acc�l�r�e sur l'axe ************/
    fichier_onaxis = fopen ("on_axis.dat","w");
    ystart=dvector(1,N);
    xp=dvector(1,kmax);
    yp=dmatrix(1,N,1,kmax);
    
    /************** Conditions initiales pour l'�lectron sur l'axe ************/
    phase = phaseo*Pi;
    ystart[1]=zf;            // z_0
    ystart[2]=0.0;           // vz_0
    ystart[3]=0.0;           // r_0
    ystart[4]=0.0;           // vr_0
    
    /************** Simulation ************************************************/
    odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
    printf("\n%s\n","Trajectoire de la particule de r�f�rence");
    printf("\n%s %13s %3d\n","successful steps:"," ",nok);
    printf("%s %20s %3d\n","bad steps:"," ",nbad);
    printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
    printf("%s %3d\n\n","stored intermediate values:    ",kount);
    
    /************** Position et �nergie finales sur l'axe *********************/
    gam_onaxis = 1/sqrt(1-yp[2][kount]*yp[2][kount]/(co*co));
    zfinal_onaxis = yp[1][kount];
    Wfinal_onaxis = gam_onaxis*m_mev;
    fprintf(fichier_onaxis, "%.12f \t %.12f \n",zfinal_onaxis,Wfinal_onaxis);
    
    fclose(fichier_onaxis);
    free_dmatrix(yp,1,N,1,kmax);
    free_dvector(xp,1,kmax);
    free_dvector(ystart,1,N);
    
    /*************** Simulation du nuage d'�lectrons **************************/
    fichier1 = fopen ("initial.dat","w");
    fichier2 = fopen ("final.dat","w");
    	for (int l=1;l<=NPTS;l++) 
        {
             ystart=dvector(1,N);
             xp=dvector(1,kmax);
             yp=dmatrix(1,N,1,kmax);
             
             /****** Tirage al�atoire des positions initiales *****************/
             /****** Le nuage est centr� � l'�tranglement zf du faisceau ******/
             /****** Pour modifier la position du centre de masse 
                     du nuage, il suffit d'ajouter ou soustraire �
                     roo et zoo un certain rini et zini.       ****************/
             roo = rscale*1.0e-6*gasdev(&idum);
             zoo = zscale*1.0e-6*gasdev(&idum)+zf;
             
             /*********** Conditions initiales ********************************/
             phase = phaseo*Pi;
             ystart[1]=zoo;     // z_0
             ystart[2]=0.0;     // vz_0
             ystart[3]=roo;     // r_0
             ystart[4]=0.0;     // vr_0
             
             /************* Simulation ****************************************/
             odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
             printf("\n%s %3d\n","Point #",l);
             printf("\n%s %13s %3d\n","successful steps:"," ",nok);
             printf("%s %20s %3d\n","bad steps:"," ",nbad);
             printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
             printf("%s %3d\n\n","stored intermediate values:    ",kount);
            
             gam = 1/sqrt(1-yp[2][kount]*yp[2][kount]/(co*co)
                                      -yp[4][kount]*yp[4][kount]/(co*co));
             angle_degres = atan(yp[4][kount]/yp[2][kount])*360/2/Pi;
             
             /******** �criture du nuage initial ******************************/
             fprintf(fichier1, "%.12f \t %.12f \t %.12f \n",
                                       (zoo-zf)*norm,roo*norm,0);
                                       
             /******** �criture du nuage final ********************************/ 
             /******** Position z (relative � zfinal_onaxis), 
                       Position r,
                       Angle par rapport � l'axe z en degr�s
                       �nergie finale en MeV    *******************************/
             fprintf(fichier2, "%.12f \t %.12f \t %.12f \t %.12f\n",
                       (yp[1][kount]-zfinal_onaxis)*norm,
                        yp[3][kount]*norm,
                        angle_degres,
                        gam*m_mev);
             
             free_dmatrix(yp,1,N,1,kmax);
             free_dvector(xp,1,kmax);
             free_dvector(ystart,1,N);
        }
         
    printf("�criture des fichiers");
    fclose (fichier1);
    fclose (fichier2);
    notes();
    printf("\nTermin�");
    return (0);
}
