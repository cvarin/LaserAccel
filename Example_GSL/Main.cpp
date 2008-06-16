#include<cstdio>
#include<cmath>

#include"SODE1D.hpp"

/******************************************************************************/
int main (void)
{    
    // Parameters for the integration
    double mu = 10.0;
    double t = 0.0, t1 = 300.0; // Time box
    double h = 1e-6;            // Starting stepsize 
    double eps_abs = 1e-6;      // Absolute precision required
    
    // Initializing the integration
    double y[2] = {1.0,1.0}; // Initial position and velocity
    
    // Integration of the ODEs
    SODE1D(t,t1,h,eps_abs,mu,y);
    
    //RPLB_Transverse_Distribution(200,3.5,"./data/RPLB_components.dat");
    
    return 0;
}

/****************** End of file ***********************************************/
