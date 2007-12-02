#ifndef INC_SODE1D_hpp
#define INC_SODE1D_hpp

#include <vector>

/****************** Structures ************************************************/
struct Motion
{
    std::vector<double> time;
    std::vector<double> position;
    std::vector<double> velocity; 
};
        
/******************************************************************************/
struct Parameters
{
    double mu;
};

/****************** Functions prototypes **************************************/
int func(double t, const double *y, double *dydt, void *params);

/******************************************************************************/
void SODE1D(double t, double t1, double h, double eps_abs,
             void *params, double *y, Motion *PtoTraj, int &count);

#endif // INC_SODE1D_hpp
