#include<vector>
#include<cmath>

#include "Constants.hpp"
#include "SODE1D.hpp"
#include "GSLODE_IO.hpp"

/******************************************************************************/
int main (void)
{    
    motion Traj;
    GSL_ExampleParams Params;
    
    // Parameters for the integration
    int count=0;                // Counts the number of points calculated
    double t = 0.0, t1 = 100.0; // Time box
    double h = 1e-6;            // Starting stepsize 
    double eps_abs = 1e-6;      // Absolute precision required
    
    // Particle
    double Wo=me_MeV+0.200; //Initial Energy = masse energy + kinetic energy
    double z0=1.0;
    double v0=sqrt(1 - pow(me_MeV/Wo,2)); //Normalized initial velocity
    
    // Initializing the integration
    Read_Input_File(&Params);
    double y[2] = {z0,v0}; // Initial position and velocity
       
    // Integration of the ODEs
    SODE1D(t,t1,h,eps_abs,Params,y,Traj,count);
    
    // Output      
    WriteToFile(Traj.time,Traj.position,Traj.velocity,count);
    WriteToGnuplot(&Traj,count);
      
    return 0;
}

/****************** End of file ***********************************************/
