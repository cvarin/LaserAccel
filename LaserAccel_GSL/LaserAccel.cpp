/*******************************************************************************

       Routines for the equations of motion time integration 
       with functions from the Gnu Scientific Library (GSL)
       Weblink : http://www.gnu.org/software/gsl/

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
/****************** Local functions prototypes ********************************/
/******************************************************************************/
double *Get_Optimums(std::vector<double> time, 
                      std::vector<double> energy, int count);

/******************************************************************************/
/****************** Accessible functions implementation ***********************/
/******************************************************************************/
void RPLB_Axial_Trajectory(double P, double wo, double T, double zf, double dzo, 
                            double lambda, double phio, double q, double m,
                             double z0, double v0, SolverParams &sp)
{
    std::vector<double> time;
    std::vector<double> energy;                     
                             
    RPLB_Params bp;                         
    RPLB_SimParams_Axial eqp;
            
    int count = 0;
    int total = 0;
    int NEQ = 2; 
    int FileFreq = 800;
    int ScreenFreq = 500000;
    int status;
    double inv_fs = 1.0e15;
    
    /************** GSL Solver initialization *********************************/
    const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s               = gsl_odeiv_step_alloc(type,NEQ);
    gsl_odeiv_control *c            = gsl_odeiv_control_y_new(sp.eps_abs,0.0);
    gsl_odeiv_evolve *e             = gsl_odeiv_evolve_alloc(NEQ);
    gsl_odeiv_system sys            = {RPLB_Axial_LorentzForce,NULL,NEQ,&eqp};
    
    /************** Simulation paramaters initialization **********************/
    double y[2] = {z0,v0};
    Set_RPLB_Params(P,wo,T,zf,dzo,lambda,phio,&bp);
    eqp.q = q;
    eqp.m = m;
    eqp.Ez = RPLB_Axial_component(y[0],sp.t,bp);

    /************** Time loop *************************************************/
    // Filling the first item of the output vectors
    time.push_back(sp.t);
    energy.push_back(me_MeV/sqrt(1-y[1]*y[1]*inv_co_square));
    while(sp.t < sp.tf)
    {  
       /*********** Update laser beam and other parameters ********************/
       eqp.Ez = RPLB_Axial_component(y[0],sp.t,bp);
       
       /*********** Accelerate and move particles *****************************/
       status = gsl_odeiv_evolve_apply(e,c,s,&sys,&sp.t,sp.tf,&sp.h,y);
       if(status != GSL_SUCCESS) break;
       
       /*********** Data are saved ********************************************/
       if(total%FileFreq == 0)
       {
          time.push_back(sp.t*inv_fs);
          energy.push_back(me_MeV/sqrt(1-y[1]*y[1]*inv_co_square));
          count++;
       }
       if(total%ScreenFreq == 0)
       {
          printf("Calculating t = %g fs (Max is %g fs)\n",
          sp.t*inv_fs,sp.tf*inv_fs);
       }
       total++;
    }
    //Puts the last item in the output vectors
    time.push_back(sp.t*inv_fs);
    energy.push_back(me_MeV/sqrt(1-y[1]*y[1]*inv_co_square));
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

/******************************************************************************/
void RPLB_3D_Trajectory(double P, double wo, double T, double zf, double dzo, 
                         double lambda, double phio, double q, double m,
                          double r0, double vr0, double z0, double vz0, 
                           SolverParams &sp)
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
    double inv_fs = 1.0e15;
    
    /************** GSL Solver initialization *********************************/
    const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s               = gsl_odeiv_step_alloc(type,NEQ);
    gsl_odeiv_control *c            = gsl_odeiv_control_y_new(sp.eps_abs,0.0);
    gsl_odeiv_evolve *e             = gsl_odeiv_evolve_alloc(NEQ);
    gsl_odeiv_system sys            = {RPLB_3D_LorentzForce,NULL,NEQ,&eqp};
    
    /************** Simulation paramaters initialization **********************/
    double y[4] = {r0,z0,vr0,vz0};
    Set_RPLB_Params(P,wo,T,zf,dzo,lambda,phio,&bp);
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
       if(y[0] < 0.0) y[0] *= -1.0; // r is always positive in cyl. coord.
       eqp.emf = RPLB_field_components(y[0],y[1],sp.t,bp);
       
       /*********** Accelerate and move particles *****************************/
       status = gsl_odeiv_evolve_apply(e,c,s,&sys,&sp.t,sp.tf,&sp.h,y);
       if(status != GSL_SUCCESS) break;
       
       /*********** Data are saved ********************************************/
       if(total%FileFreq == 0)
       {
          time.push_back(sp.t*inv_fs);
          energy.push_back(me_MeV/sqrt(1-y[3]*y[3]*inv_co_square));
          count++;
       }
       if(total%ScreenFreq == 0)
       {
          printf("Calculating t = %g fs (Max is %g fs)\n",
          sp.t*inv_fs,sp.tf*inv_fs);
       }
       total++;
    }
    //Puts the last item in the output vectors
    time.push_back(sp.t*inv_fs);
    energy.push_back(me_MeV/sqrt(1-y[3]*y[3]*inv_co_square));
    count++;
    
    /************** Data written to file **************************************/
    std::ofstream OutFile("./data/RPLB3D.dat", std::ios::out);
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

/******************************************************************************/
double *RPLB_Phase_Scan(double P, double wo, double T, double zf, double dzo, 
                       double lambda, double q, double m, double z0, double v0)
{
    printf("Starting RPLB_Phase_Scan...\n");                   
    std::vector<double> phase;
    std::vector<double> energy;                     
                             
    RPLB_Params bp;                         
    RPLB_SimParams_Axial eqp;
    SolverParams sp;
    SolverParams sp_initial;
            
    int count = 0;
    int total = 0;
    int NEQ = 2; 
    int ScreenFreq = 500000;
    int status;
    double inv_fs = 1.0e15;
    
    /************** GSL Solver initialization *********************************/
    sp.t = 0.0, sp.tf = 20.0e-12; // Time box
    sp.h = 1.0e-15;               // Starting stepsize 
    sp.eps_abs = 2.0e-11;         // Absolute precision required
    sp_initial = sp;
        
    const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s               = gsl_odeiv_step_alloc(type,NEQ);
    gsl_odeiv_control *c            = gsl_odeiv_control_y_new(sp.eps_abs,0.0);
    gsl_odeiv_evolve *e             = gsl_odeiv_evolve_alloc(NEQ);
    gsl_odeiv_system sys            = {RPLB_Axial_LorentzForce,NULL,NEQ,&eqp};
    
    /************** Simulation paramaters initialization **********************/
    double y[2] = {z0,v0};
    double W;
    double phi_min = 0.0;
    double phi_max = 2.0;
    double phi = phi_min;
    int N = 100;
    eqp.q = q;
    eqp.m = m;

    /************** Phase scan ************************************************/
    std::ofstream OutFile("./data/RPLB_PhaseScan.dat", std::ios::out);
    OutFile.precision(16);
    while(phi < phi_max)
    {   
        printf("Simulating for phi = %.2f PI\n",phi);
        while(sp.t < sp.tf)
        {
            /*********** Update laser beam and other parameters ***************/
            Set_RPLB_Params(P,wo,T,zf,dzo,lambda,phi*Pi,&bp);
            eqp.Ez = RPLB_Axial_component(y[0],sp.t,bp);
            
            /*********** Accelerate and move particles ************************/
            status = gsl_odeiv_evolve_apply(e,c,s,&sys,&sp.t,sp.tf,&sp.h,y);
            if(status != GSL_SUCCESS) break;
            if(total%ScreenFreq == 0)
            {
               printf("\tCalculating t = %g fs (Max is %g fs)\n",
               sp.t*inv_fs,sp.tf*inv_fs);
            }
            total++;
        }
        
        /********** Prepare for the next step *********************************/
        W = me_MeV/sqrt(1-y[1]*y[1]*inv_co_square);
        phase.push_back(phi);
        energy.push_back(W);
        count++;
        OutFile << phi << "\t" << W << "\n";
        phi += (phi_max-phi_min)/N;
        sp = sp_initial;
        y[0] = z0; 
        y[1] = v0;
        total = 0;
    }
    OutFile.close();
    
    /************** Memory is freed *******************************************/
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    return Get_Optimums(phase,energy,count);
}

/******************************************************************************/
/****************** Local functions implementation ****************************/
/******************************************************************************/
double *Get_Optimums(std::vector<double> phase, 
                      std::vector<double> energy, int count)
{
     int i;
     double *opt = (double*)calloc(4,sizeof(double));
     opt[1] = energy[0];
     opt[3] = energy[0];
     for(i=1;i<count;i++)
     {
         if(energy[i] > opt[1]) opt[0] = phase[i],opt[1] = energy[i];
         if(energy[i] < opt[3]) opt[2] = phase[i],opt[3] = energy[i];
     }
     return opt;
     free(opt);
}

/****************** End of file ***********************************************/
