/*******************************************************************************

       Integration follows the tutorial of the Gnu Scientific Library
       Weblink : http://www.gnu.org/software/gsl//manual/html_node/

*******************************************************************************/
#include <fstream>
#include <iostream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include"Gnuplot_i.hpp"
#include"SODE1D.hpp"

using namespace std;

/******************************************************************************/
/****************** Local functions prototypes ********************************/
/******************************************************************************/
void SODE1D_toFile(vector<double> t, vector<double> p, 
                    vector<double> v, int count);
void SODE1D_toGnuplot(const motion1D *m, const int count);

/******************************************************************************/
// The ordinary differential equations to solve
int GSL_Example(double t, const double y[], double dydt[], void *param)
{
    double *ptoparam = (double*)param;
    double mu = *ptoparam;                      
        
    dydt[0] = y[1];
    dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

/******************************************************************************/
void SODE1D(double t, double t1, double h, 
             double eps_abs, double mu, double y[])
{ 
    motion1D m;  
    int count = 0;       
    int status;
    int Neq = 2;
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s            = gsl_odeiv_step_alloc(T,Neq);
    gsl_odeiv_control *c         = gsl_odeiv_control_y_new(eps_abs,0.0);
    gsl_odeiv_evolve *e          = gsl_odeiv_evolve_alloc(Neq);
    gsl_odeiv_system sys         = {GSL_Example,NULL,Neq,&mu};
    
    // Filling the first item of the output vectors
    m.time.push_back(t);
    m.position.push_back(y[0]);
    m.velocity.push_back(y[1]);
    
    // Time loop
    while(t < t1)
    {
       status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,t1,&h,y);   
       if(status != GSL_SUCCESS) break;  
       m.time.push_back(t);
       m.position.push_back(y[0]);
       m.velocity.push_back(y[1]);
       count++;         
    }
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s); 
    
    // Output      
    SODE1D_toFile(m.time,m.position,m.velocity,count);
    SODE1D_toGnuplot(&m,count); 
}

/******************************************************************************/
/****************** Local functions implementations ***************************/
/******************************************************************************/
void SODE1D_toFile(vector<double> t, vector<double> p, 
                    vector<double> v, int count)
{      
     ofstream OutFile("./data/SODE1D.dat", ios::out);    
     for(int i=0 ; i<=count ; i++)
     {
        OutFile << t[i]  << "\t" 
                << p[i]  << "\t" 
                << v[i]  << "\n";
     }
     OutFile.close();
}

/******************************************************************************/
void SODE1D_toGnuplot(const motion1D *m, const int count)
{     
     Gnuplot Plot;
     
     Plot.set_style("lines");
     Plot.set_xlabel("Temps");
     Plot.set_ylabel("Amplitude");
     Plot.plot_xy(m->time,m->position,"Position");
     Plot.plot_xy(m->time,m->velocity,"Vitesse");
    
     puts("Press [ENTER] or [CTRL-C] to quit.");
     getchar();    
}

/****************** End of file ***********************************************/
