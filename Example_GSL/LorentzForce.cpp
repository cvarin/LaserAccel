/*******************************************************************************

         Various expressions of the Lorentz force equation

*******************************************************************************/
#include<cmath>
#include <gsl/gsl_errno.h>

#include"Constants.hpp"
#include"LorentzForce.hpp"

/****************** External functions implementation *************************/
int RPLB_3D_LorentzForce(const double t, double y[], 
                          double dydt[], void* p)
{
    /************** Extract parameters from the external structure ************/ 
    RPLB_SimParams *pp = (RPLB_SimParams*)p;
    const double q_m = pp->q/pp->m;
    const double Er = pp->emf.Er;
    const double Ez = pp->emf.Ez;
    const double B_theta = pp->emf.B_theta;
    
    /************** Initialize some local variables ***************************/
    const double r = y[0];
    const double z = y[1];
    const double vr = y[2];
    const double vz = y[3];
    const double inv_gamma = sqrt(1 - (vr*vr)*inv_co_square
                                     - (vz*vz)*inv_co_square);
    const double vdotE = vr*Er + vz*Ez;
    
    /************** Equations of motion ***************************************/
    dydt[0] = vr;
    dydt[1] = vz;
    dydt[2] = q_m*inv_gamma*(Er - vz*B_theta - vr*vdotE*inv_co_square);
    dydt[3] = q_m*inv_gamma*(Ez + vr*B_theta - vz*vdotE*inv_co_square);
    
    /************** Returns a positive value for r ****************************/
    if(r < 0.0) y[0] = -y[0];

    return GSL_SUCCESS;
}

/****************** End of file ***********************************************/
