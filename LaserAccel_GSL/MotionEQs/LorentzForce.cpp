/*******************************************************************************

         Various expressions of the Lorentz force equation

*******************************************************************************/
#include <cmath>
#include <gsl/gsl_errno.h>

#include "Constants.hpp"
#include "LorentzForce.hpp"

/****************** External functions implementation *************************/
int RPLB_3D_LorentzForce(double t, const double y[], double dydt[], void* p)
{
    /* 
         Note : r = y[0], z = y[1], vr = y[2], vz = y[3].
    */
    /************** Extract parameters from the external structure ************/ 
    RPLB_SimParams *pp = (RPLB_SimParams*)p;
    const double q_m = pp->q/pp->m;
    const double Er = pp->emf.Er;
    const double Ez = pp->emf.Ez;
    const double B_theta = pp->emf.B_theta;
    
    /************** Initialize some local variables ***************************/
    const double inv_gamma = sqrt(1 - y[2]*y[2]*inv_co_square
                                     - y[3]*y[3]*inv_co_square);
    const double vdotE = y[2]*Er + y[3]*Ez;
    
    /************** Equations of motion ***************************************/
    dydt[0] = y[2];
    dydt[1] = y[3];
    dydt[2] = q_m*inv_gamma*(Er - y[3]*B_theta - y[2]*vdotE*inv_co_square);
    dydt[3] = q_m*inv_gamma*(Ez + y[2]*B_theta - y[3]*vdotE*inv_co_square);
    
    return GSL_SUCCESS;
}

/******************************************************************************/
int RPLB_Axial_LorentzForce(double t, const double y[], double dydt[], void* p)
{
    /*
         Note : z = y[0], vz = y[1].
    */
    /************** Extract parameters from the external structure ************/ 
    RPLB_SimParams_Axial *pp = (RPLB_SimParams_Axial*)p;
    
    /************** Initialize some local variables ***************************/
    const double diff_v = 1 - y[1]*y[1]*inv_co_square;
    const double inv_gamma = sqrt(diff_v);
    
    /************** Equations of motion ***************************************/
    dydt[0] = y[1];
    dydt[1] = (pp->q/pp->m)*inv_gamma*diff_v*pp->Ez;
    
    return GSL_SUCCESS;
}

/****************** End of file ***********************************************/
