/*******************************************************************************

         Various expressions of the Lorentz force equation

*******************************************************************************/
#include<cmath>
#include <gsl/gsl_errno.h>

#include"Constants.hpp"
#include"LorentzForce.hpp"

/****************** External functions implementation *************************/
int RPLB_3D_LorentzForce(const double t, const double y[], 
                          double dydt[], void* p)
{
    RPLB_SimParams *pp = (RPLB_SimParams*)p;
    const double q = pp->q;
    const double m = pp->m;
    
    const double r = y[0];
    const double z = y[1];
    const double vr = y[2];
    const double vz = y[3];
    const double inv_gamma = sqrt(1 - (vr*vr)*inv_co_square
                                     - (vz*vz)*inv_co_square);
    
    dydt[0] = vr;
    dydt[1] = vz;
    dydt[2] = 0.0;
    dydt[3] = 0.0;
    
    return GSL_SUCCESS;
}

/****************** End of file ***********************************************/
