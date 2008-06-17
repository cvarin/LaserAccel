#ifndef _INCLUDE_RPLBeam_hpp
#define _INCLUDE_RPLBeam_hpp

/****************** Structures ************************************************/
struct RPLB_Params
{
    /************** Provided values *******************************************/
    double P;            // Laser power [W]
    double wo;           // Beam spot size at the waist [m]
    double T;            // Pulse duration [s]
    double dzo;          // Initial position of the center of the pulse [m]
    double lambda;       // Wavelength [m]
    double phi0;        // Field phase at beam waist [rad]
    
    /************** Derived values ********************************************/
    double k;            // Wavenumber [rad/m]
    double omega;        // Angular frequency [rad/s]
    double z_Rayleigh;   // Rayleigh distance [m]
    double I;            // Intensity [W/cm^2]
    double W;            // Pulse Energy [J]
    double Eo;           // Transverse electric component amplitude [V/m]
    double Ez_norm;      // Long. E-field normalization factor
    double B_theta_norm; // B-field normalization factor
    double ao;           // Normalized transv. E-field strength
    double az;           // Normalized long. E-field strength
};

/******************************************************************************/
struct Cartesian_Coordinates
{
    double x,y,z,t;
};

/******************************************************************************/
struct Cartesian_EMfield
{
    // Field components of a general laser beam in cartesian coordinates
    double Ex,Ey,Ez;
    double Bx,By,Bz;
};

/******************************************************************************/
struct Cylindric_Circular_Coordinates
{
    double r,theta,z,t;
};

/******************************************************************************/
struct RPLB_EMfield
{
    // Field components of a radially polarized laser beam
    double Er,Ez,B_theta;     
};

/****************** Accessible functions prototypes ***************************/
// Radially polarized laser beam (RPLB)
RPLB_EMfield RPLB_field_components(double r, double z, double t,
                                    const RPLB_Params bp);
double RPLB_Axial_component(double z, double t, const RPLB_Params bp);
void Set_RPLB_Params(double P, double wo, double T, double dzo, double lambda, 
                      double phi_0, RPLB_Params *bp);
void RPLB_Transverse_Distribution(int N, double ro, const char *filename);

#endif // End of file
