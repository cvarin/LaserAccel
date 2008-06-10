/*******************************************************************************

         Set of functions that gives the electromagnetic field 
         components of laser beams.

*******************************************************************************/
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "Constants.hpp"
#include "LaserBeams.hpp"

/******************************************************************************/
RPLB_EMfield RPLB_field_components(double r, double z, double t,
                                    const RPLB_Params bp)
{
    /************** Declarations **********************************************/ 
    RPLB_EMfield f;
    
    /************** Initializations *******************************************/
    const double TINY = 1.0e-30;
    const double exp_1_2 = exp(0.5);
    const double Eo = bp.Eo;
    const double ko = bp.k;
    const double omega = bp.omega;
    const double wo = bp.wo;
    const double T  = bp.T;
    const double zR = bp.z_Rayleigh;
    const double dzo = bp.dzo;
    const double phi_0 = bp.phi_0;
    const double w  = wo*sqrt(1 + (z*z)/(zR*zR));
    const double R  = z + (zR*zR)/z;
    const double Psi_G = atan(z/zR);
    const double Psi_C = (z <= TINY)? 0.0 : 0.5*ko*r*r/R;
    const double Psi = omega*t - ko*z + 2*Psi_G - Psi_C - phi_0;
    const double t_prime = t - (z - dzo)*inv_co;
    const double pulse_enveloppe = exp(-t_prime*t_prime/(T*T));
    const double spatial_enveloppe = (wo*wo)/(w*w)*exp(-(r*r)/(w*w));
    
    /************** Field components are calculated ***************************/
    f.Er = Eo*exp_1_2*(sqrt_2*r/wo)*spatial_enveloppe*pulse_enveloppe*cos(Psi);
    f.Ez = Eo*exp_1_2*(2*sqrt_2/(ko*wo))*spatial_enveloppe*pulse_enveloppe
              *((1 - r*r/(w*w))*sin(Psi) - Psi_C*cos(Psi));
    f.B_theta = f.Er*inv_co;
    
    return f;
}                                    

/******************************************************************************/
void Set_RPLB_Params(double P, double wo, double T, double dzo, double lambda, 
                      double phi_0, RPLB_Params *bp)
{
    /************** Provided values *******************************************/
    bp->P = P;            // Laser power [W]
    bp->wo = wo;          // Beam spot size at the waist [m]
    bp->T = T;            // Pulse duration [s]
    bp->dzo = dzo;        // Initial position of the center of the pulse [m]
    bp->lambda = lambda;  // Wavelength [m]
    bp->phi_0 = phi_0;    // Field phase at beam waist [rad]
    
    /************** Derived values ********************************************/
    bp->k = 2.0*Pi/lambda;            // Wavenumber [rad/m]
    bp->omega = bp->k*co;             // Angular frequency [rad/s]
    bp->z_Rayleigh = 0.5*bp->k*wo*wo; // Rayleigh distance [m]
    bp->I = 2*exp(-1.0)*P/(Pi*wo*wo); // Intensity [W/m^2]
    bp->W = 0.5*Pi*P*T;               // Pulse Energy [J]
    bp->Eo = sqrt(2*eta_0*bp->I); // Transverse electric component amplitude [V/m]
    bp->Ez_norm = 0.5*bp->k*wo*exp(-0.5)/sqrt_2;// Long. E-field normalization factor
    bp->B_theta_norm = co; // B-field normalization factor
    bp->ao = (q_e*bp->Eo)/(me_kg*co*bp->omega);// Normalized transv. E-field strength
    bp->az = bp->ao/bp->Ez_norm; // Normalized long. E-field strength
}

/******************************************************************************/
void RPLB_Transverse_Distribution(int N, double ro, const char *filename)
{
    RPLB_Params bp;
    RPLB_EMfield champ1;
    RPLB_EMfield champ2;
    double r = 0.0;
    const char *pathfile = "./Visualization/RPLB_components.path";
    
    /**************************************************************************/
    FILE *datapath = fopen(pathfile,"w");
    FILE *file = fopen(filename,"w");
    if(file==NULL || datapath==NULL)
    {
         if(file==NULL) printf("Can't open \'%s\'.\n",filename);
         if(datapath==NULL) printf("Can't open \'%s\'.\n",pathfile);
         getchar();
         exit(1);
    }
    
    /************** Export output file location *******************************/
    fprintf(datapath,"datafile = \'%s\'",filename);
    fclose(datapath);
    
    /************** Write data file *******************************************/
    Set_RPLB_Params(0.0,1.0,1.0,0.0,1.0,0.0,&bp);
    bp.Eo = 1.0;
    while(r < ro)
    {
         bp.phi_0 = 0.0;
         champ1 = RPLB_field_components(r,0.0,0.0,bp);
         bp.phi_0 = -0.5*Pi;
         champ2 = RPLB_field_components(r,0.0,0.0,bp);
         fprintf(file,"%e\t%e\t%e\t%e\n",
                  r, champ1.Er, champ2.Ez*bp.Ez_norm,
                    champ1.B_theta*bp.B_theta_norm);
         r+=2*ro/N;
    }
    fclose(file);
    system("gnuplot Visualization/RPLB_components.gp");
}

/****************** End of file ***********************************************/
