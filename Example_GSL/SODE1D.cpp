/*******************************************************************************

       Integration follows the tutorial of the Gnu Scientific Library
       Weblink : http://www.gnu.org/software/gsl//manual/html_node/

*******************************************************************************/
#include <iostream>
#include <cstdio>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include"SODE1D.hpp"

/******************************************************************************/
// The ordinary differential equations to solve
int GSL_Example(double t, const double y[], double dydt[], void* params)
{
    GSL_ExampleParams *PtoParams = (GSL_ExampleParams*)params;
    double mu = PtoParams->mu;                      
        
    dydt[0] = y[1];
    dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

/******************************************************************************/
void SODE1D(double t, double t1, double h, double eps_abs,
             GSL_ExampleParams params, double y[], motion &traj, int &count)
{ 
    int status;
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * s            = gsl_odeiv_step_alloc(T,2);
    gsl_odeiv_control * c         = gsl_odeiv_control_y_new(eps_abs,0.0);
    gsl_odeiv_evolve * e          = gsl_odeiv_evolve_alloc(2);
    gsl_odeiv_system sys          = {GSL_Example,NULL,2,&params};
    
    // Filling the first item of the output vectors
    traj.time.push_back(t);
    traj.position.push_back(y[0]);
    traj.velocity.push_back(y[1]);
    // Time loop
    while(t < t1)
    {
       status = gsl_odeiv_evolve_apply (e,c,s,&sys,&t,t1,&h,y);   
       if(status != GSL_SUCCESS) break;
       traj.time.push_back(t);
       traj.position.push_back(y[0]);
       traj.velocity.push_back(y[1]);  
       count++;         
    }
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);  
}

/****************** End of file ***********************************************/
