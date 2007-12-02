#include<iostream>
#include<cmath>

#include "PhysConsts.hpp"
#include "SODE1D.hpp"
#include "WriteTo.hpp"

//// Laser
//double P = 1e12;       // Laser power [W]
//double lambda = 800e-9;// Wavelength [m]
//double wo = 3e-6;      // Beam waist [m]
//double T = 10e-15;     // Pulse duration [s]
//double phi_0 = 0.0;    // Field phase at beam waist [rad]
//double zp0 = 0.0;      // Initial position of the center of the pulse [m]
//double z_Rayleigh = Pi*pow(wo,2)/2; // Rayleigh distance [m]
//

int main (void)
{    
    Motion traj;   
    Parameters params;
    
    // Parameters for the integration
    int count=0;                // Counts the number of points calculated
    double t = 0.0, t1 = 100.0; // Time box
    double h = 1e-6;            // Starting stepsize 
    double eps_abs = 1e-6;      // Absolute precision required
    
    // Particle
    double Wo  = me_MeV + 0.200;     //Energy = mass + kinetic (all in [MeV])
    double ze0 = 1.0;             // Initial position of the electron [m]
    double v0  = sqrt( 1 - pow(me_MeV/Wo,2) ); //Normalized initial velocity
    
    // Setting the initial conditions
    double y[2] = {ze0,v0}; // Initial position and velocity
    params.mu = 10;
    
    // Filling the first item of the output vectors
    traj.time.push_back(t);
    traj.position.push_back(y[0]);
    traj.velocity.push_back(y[1]);
    
    // Integration of the ODEs
    SODE1D(t,t1,h,eps_abs,&params,y,&traj,count);
    
    // Graphical Output      
    WriteToGnuplot(traj,count);
    
    return 0;
}
