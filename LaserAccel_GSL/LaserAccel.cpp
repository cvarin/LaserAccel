/*******************************************************************************

       Routines for the equations of motion time integration 
       following the tutorial of the Gnu Scientific Library (GSL)
       Weblink : http://www.gnu.org/software/gsl//manual/html_node/

*******************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "Constants.hpp"
#include "LaserAccel.hpp"
#include "MotionEQs.hpp"

/******************************************************************************/
void RPLB_Axial_Trajectory(double P, double wo, double T, double dzo, 
                            double lambda, double phio, double q, double m,
                             double z0, double v0, SolverParams &sp)
{
    std::vector<double> time;
    std::vector<double> energy;                     
                             
    RPLB_Params bp;                         
    RPLB_SimParams eqp;
            
    int count = 0;
    int total = 0;
    int NEQ = 4; 
    int FileFreq = 800;
    int ScreenFreq = 500000;
    int status;
    
    /************** GSL Solver initialization *********************************/
    const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s               = gsl_odeiv_step_alloc(type,4);
    gsl_odeiv_control *c            = gsl_odeiv_control_y_new(sp.eps_abs,0.0);
    gsl_odeiv_evolve *e             = gsl_odeiv_evolve_alloc(4);
    gsl_odeiv_system sys            = {RPLB_3D_LorentzForce,NULL,NEQ,&eqp};
    
    /************** Simulation paramaters initialization **********************/
    double y[4] = {0.0,z0,0.0,v0};
    Set_RPLB_Params(P,wo,T,dzo,lambda,phio,&bp);
    eqp.q = q;
    eqp.m = m;
    eqp.emf = RPLB_field_components(y[0],y[1],sp.t,bp);

    /************** Time loop *************************************************/
    // Filling the first item of the output vectors
    time.push_back(sp.t);
    energy.push_back(me_MeV/sqrt(1-y[3]*y[3]*inv_co_square));
    while(sp.t < sp.tf)
    {
       /*********** Update laser beam and other parameters ********************/
       eqp.emf = RPLB_field_components(y[0],y[1],sp.t,bp);
       
       /*********** Accelerate and move particles *****************************/
       status = gsl_odeiv_evolve_apply(e,c,s,&sys,&sp.t,sp.tf,&sp.h,y);
       if(status != GSL_SUCCESS) break;
       
       /*********** Data are saved ********************************************/
       if(total%FileFreq == 0)
       {
          time.push_back(sp.t*1.0e15);
          energy.push_back(me_MeV/sqrt(1-y[3]*y[3]*inv_co_square));
          count++;
       }
       if(total%ScreenFreq == 0)
       {
          printf("Calculating t = %g fs (Max is %g fs)\n",sp.t*1e15,sp.tf*1e15);
       }
       total++;
    }
    //Puts the last item in the output vectors
    time.push_back(sp.t*1.0e15);
    energy.push_back(me_MeV/sqrt(1-y[3]*y[3]*inv_co_square));
    count++;
    
    /************** Data written to file **************************************/
    std::ofstream OutFile("./data/RPLBAxial.dat", std::ios::out);
    OutFile.precision(16);
    for(int i=0 ; i<=count ; i++)
    {
        OutFile << time[i]  << "\t" 
                << energy[i]  << "\n";
    }
    OutFile.close();
    
    /************** Memory is freed *******************************************/
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s); 
}

/****************** End of file ***********************************************/

