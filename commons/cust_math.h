#ifndef inc_cust_math_h
#define inc_cust_math_h

#include <complex.h>
#include <math.h>

/******************************************************************************/
inline complex double spherical_j0(complex double z)
{
     return csin(z)/z;    
}

/******************************************************************************/
inline complex double spherical_j1(complex double z)
{
     return csin(z)/(z*z) - ccos(z)/z;
}

/******************************************************************************/
inline complex double spherical_j2(complex double z)
{
     return (3.0/(z*z*z) - 1.0/z)*csin(z) - 3.0*ccos(z)/(z*z);
}

#endif /* inc_cust_math_h */