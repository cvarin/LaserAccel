#ifndef _INCLUDE_LorentzForce_hpp
#define _INCLUDE_LorentzForce_hpp

#include "LaserBeams.hpp"

/****************** Structures ************************************************/
struct RPLB_SimParams
{
    double m;
    double q;   
    RPLB_EMfield emf;
};

/****************** Lorentz force equations ***********************************/
int RPLB_3D_LorentzForce(double t, const double y[], double dydt[], void* p);

#endif // End of file
