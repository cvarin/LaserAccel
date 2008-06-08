#include<cstdio>
#include<vector>
#include<cmath>

#include"Constants.hpp"
#include"GSLODE_IO.hpp"
#include"LaserBeams.hpp"
#include"SODE1D.hpp"

//    // Laser beam parameters
//    double P = 1e12;        // Laser power Watts
//    double lambda = 800e-9; // Wavelength in m
//    double wo = 3e-6;       // Beam waist in m
//    double T = 10e-15;      // Pulse duration in s
//    double phi_0 = 0.0;     // Field phase at beam waist
//    double zp0 = 0.0;       // Initial position of the center of the pulse
//    double z_Rayleigh = Pi*pow(wo,2)/2; // Rayleigh distance

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
    //WriteToFile(Traj.time,Traj.position,Traj.velocity,count);
    //WriteToGnuplot(&Traj,count);
    
    //Write_RPLB_Transverse_Distribution(300,5.0,"RPLB_components.dat");
    
    return 0;
}

/****************** End of file ***********************************************/
