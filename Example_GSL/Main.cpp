#include<cstdio>
#include<vector>
#include<cmath>

#include"Constants.hpp"
#include"GSLODE_IO.hpp"
#include"LaserBeams.hpp"
#include"SODE1D.hpp"

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
    double v0=sqrt(1 - (me_MeV*me_MeV)/(Wo*Wo)); //Normalized initial velocity
    
    // Initializing the integration
    Read_Input_File(&Params,"ModelParams.ini");
    double y[2] = {1.0,1.0}; // Initial position and velocity
       
    // Integration of the ODEs
    SODE1D(t,t1,h,eps_abs,Params,y,Traj,count);
    
    // Output      
    WriteToFile(Traj.time,Traj.position,Traj.velocity,count);
    WriteToGnuplot(&Traj,count);
    
    //RPLB_Transverse_Distribution(200,3.5,"./data/RPLB_components.dat");
    
    return 0;
}

/****************** End of file ***********************************************/
