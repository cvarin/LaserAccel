#ifndef _INCLUDE_LaserBeams_hpp
#define _INCLUDE_LaserBeams_hpp

/****************** Structures ************************************************/
struct RPLB_Params
{
    double P;            // Laser power [W]
    double I;            // Intensity [W/cm^2]
    double W;            // Pulse Energy [J]
    double Eo;           // Transverse electric component amplitude
    double Er_norm;      // Normalization factor
    double Ez_norm;      // Normalization factor
    double B_theta_norm; // Normalization factor
    double lambda;       // Wavelength [m]
    double ko;           // Wavenumber [rad/m]
    double omega;        // Angular frequency [rad/s]
    double wo;           // Beam spot size at the waist [m]
    double T;            // Pulse duration [s]
    double phi_0;        // Field phase at beam waist [rad]
    double dzo;          // Initial position of the center of the pulse [m]
    double z_Rayleigh;   // Rayleigh distance [m]
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
RPLB_EMfield RPLB_field_components(double r, double z, double t,
                                    const RPLB_Params bp);
void Write_RPLB_Transverse_Distribution(int N, double ro, const char *filename);

#endif // End of file
