/*******************************************************************************

    ODE integration routine. See the tutorial of the Gnu Scientific Library
    Weblink : http://www.gnu.org/software/gsl//manual/html_node/

*******************************************************************************/

#include <iostream>
#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "PhysConsts.hpp"
#include "SODE1D.hpp"

/******************************************************************************/
// The ordinary differential equations to solve
int func(double t, const double *y, double *dydt, void *params)
{
    Parameters *PtoParams = (Parameters *)params;
    double mu = PtoParams->mu;                      
        
    dydt[0] = y[1];
    dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}
    
/******************************************************************************/
// Ode Solver as from the Gnu Scientific Library
void SODE1D(double t, double t1, double h, double eps_abs,
             void *params, double *y, Motion *PtoTraj, int &count)
{ 
    Parameters *PtoParams = (Parameters *)params;
    
    std::ofstream OutFile;
    OutFile.open("Run.dat");
    
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * s            = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c         = gsl_odeiv_control_y_new (eps_abs, 0.0);
    gsl_odeiv_evolve * e          = gsl_odeiv_evolve_alloc (2);
    
    gsl_odeiv_system sys = {func,NULL,2,params};
    
    while (t < t1)
    {
        int status = gsl_odeiv_evolve_apply (e,c,s,&sys,&t,t1,&h,y);
        
        if (status != GSL_SUCCESS) break;
        
        // Filling the external motion vectors
        count++; // Counts the number of points    
        PtoTraj->time.push_back(t);
        PtoTraj->position.push_back(y[0]);
        PtoTraj->velocity.push_back(y[1]); 
        
        // Writing to a file
        OutFile << t    << "\t"
                << y[0] << "\t"
                << y[1] << std::endl;            
    }
    
    OutFile.close();
    
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);  
}

/****************** End of file ***********************************************/