#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "Constants.hpp"
#include "LaserAccel.hpp"

/****************** Local functions *******************************************/
inline void ExecTime(clock_t dt)
{
    printf("Execution time : %d clock cycles (%f s).\n",
                              dt,(double)dt/CLOCKS_PER_SEC);
}

inline void Print_Optimal_Phases(double *OPT)
{
    printf("Optimal phase\n");
    printf("Phi_max = %f Pi rad, W_max = %f MeV.\n",OPT[0],OPT[1]);
    printf("Phi_min = %f Pi rad, W_min = %f MeV.\n",OPT[2],OPT[3]);
}

/******************************************************************************/
int main (void)
{    
    SolverParams sp;
    double *OPT = (double*)calloc(4,sizeof(double));
    
    // Parameters for the integration
    sp.t = 0.0, sp.tf = 30.0e-12; // Time box
    sp.h = 1.0e-15;               // Starting stepsize 
    sp.eps_abs = 2.0e-12;         // Absolute precision required

    // Number of points for RPLB_Phase_Scan()    
    int Npts = 30;
    
    // Laser Beam
    double P = 100.0e12;
    double wo = 3.0e-6;
    double T = 10.0e-15;
    double zf = 1.0;
    double dzo = -4.0*T*co;
    double lambda = 0.8e-6;
    double phio = 1.3*Pi;
    
    // Particle
    int N = 10; // Number of particle for RPLB_Axial_MultiParticle()
    double q = q_e;
    double m = me_kg;
    double r0 = 0.0; double z0 = zf+0.0;
    double vr0 = 0.0; double vz0 = 0.0;
    double cloud_width = 0.5*lambda; 
    
    //RPLB_Axial_Trajectory(P,wo,T,zf,dzo,lambda,phio,q,m,z0,vz0,sp);
    RPLB_Axial_MultiParticle(P,wo,T,zf,dzo,lambda,phio,q,m,N,z0,cloud_width,sp);
    //RPLB_3D_Trajectory(P,wo,T,zf,dzo,lambda,phio,q,m,r0,vr0,z0,vz0,sp);
    
    //OPT = RPLB_Phase_Scan(P,wo,T,zf,dzo,lambda,q,m,z0,vz0,sp,Npts);
    //Print_Optimal_Phases(OPT);
    
    //RPLB_Transverse_Distribution(200,3.5,"./data/RPLB_components.dat");
    free(OPT);
            
    /**************************************************************************/
    ExecTime(clock());
    return 0;
}

/****************** End of file ***********************************************/
