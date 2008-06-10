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
    // Parameters for the integration
//    double t = 0.0, t1 = 100.0; // Time box
//    double h = 1e-6;            // Starting stepsize 
//    double eps_abs = 1e-6;      // Absolute precision required
    
    // Particle
//    double Wo=me_MeV+0.200; //Initial Energy = masse energy + kinetic energy
//    double z0=1.0;
//    double v0=sqrt(1 - (me_MeV*me_MeV)/(Wo*Wo)); //Normalized initial velocity
    
    RPLB_Transverse_Distribution(200,3.5,"./data/RPLB_components.dat");
    
    /**************************************************************************/
    ExecTime(clock());
    return 0;
}

/****************** End of file ***********************************************/
