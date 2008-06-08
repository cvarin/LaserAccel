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
    const double ko = bp.ko;
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
void Set_RPLB_Params(RPLB_Params *beam_params)
{
     
}

/******************************************************************************/
void Write_RPLB_Transverse_Distribution(int N, double ro, const char *filename)
{
    RPLB_Params bp;
    RPLB_EMfield champ;
    double r = 0.0;
    
    bp.P = 0.0;          // Laser power [W]
    bp.I = 0.0;          // Intensity [W/cm^2]
    bp.W = 0.0;          // Pulse Energy [J]
    bp.Eo = 1.0;         // Transverse electric component amplitude
    bp.lambda = 1.0;     // Wavelength [m]
    bp.ko = 2*Pi;         // Wavenumber [rad/m]
    bp.omega = bp.ko*co;      // Angular frequency [rad/s]
    bp.wo = 1.0;         // Beam spot size at the waist [m]
    bp.T = 1.0;          // Pulse duration [s]
    bp.phi_0 = -0.0*Pi;      // Field phase at beam waist [rad]
    bp.dzo = 0.0;        // Initial position of the center of the pulse [m]
    bp.z_Rayleigh = 0.5*bp.ko*bp.wo*bp.wo; // Rayleigh distance [m]
    
    FILE *file = fopen(filename,"w");
    if(file==NULL)
    {
         printf("Can't open %s.\n",filename);
         getchar();
         exit(1);
    }
    while(r < ro)
    {
         champ = RPLB_field_components(r,0.0,0.0,bp);
         fprintf(file,"%e\t%e\t%e\t%e\n",r,champ.Er,champ.Ez,champ.B_theta);
         r+=2*ro/N;
    }
    fclose(file);
}

/****************** End of file ***********************************************/
