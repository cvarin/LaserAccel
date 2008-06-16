#include<cstdio>
#include<cmath>
#include<ctime>

#include"Constants.hpp"
#include"LaserAccel.hpp"

/****************** Local functions *******************************************/
inline void ExecTime(clock_t dt)
{
    printf("Execution time : %d clock cycles (%f s).\n",
                              dt,(double)dt/CLOCKS_PER_SEC);
}

/******************************************************************************/
int main (void)
{    
     SolverParams sp;
    
    // Parameters for the integration
     sp.t = 0.0, sp.tf = 40.0e-12; // Time box
     sp.h = 1.0e-15;               // Starting stepsize 
     sp.eps_abs = 1.0e-12;         // Absolute precision required
    
    // Particle
    double q = q_e;
    double m = me_kg;
    double Wo=me_MeV+0.0; //Initial Energy = masse energy + kinetic energy
    double z0=0.0;
    double v0=sqrt(1 - (me_MeV*me_MeV)/(Wo*Wo)); //Normalized initial velocity
    
    // Laser Beam
    double P = 100.0e12;
    double wo = 3.0e-6;
    double T = 10.0e-15;
    double dzo = -2.0*T*co;
    double lambda = 0.8e-6;
    double phio = 0.0*Pi;
    
    RPLB_Axial_Trajectory(P,wo,T,dzo,lambda,phio,q,m,z0,v0,sp);
    
    //RPLB_Transverse_Distribution(200,3.5,"./data/RPLB_components.dat");
        
    /**************************************************************************/
    ExecTime(clock());
    return 0;
}

/****************** End of file ***********************************************/
