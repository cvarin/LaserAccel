using namespace std;

int func (double t, const double y[], double dydt[], void* params);

void SODE1D(double t, double t1, double h, double eps_abs,
     void* params, double y[], void* trajectory, int& count);
