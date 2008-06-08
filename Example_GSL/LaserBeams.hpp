#ifndef _INCLUDE_LaserBeams_hpp
#define _INCLUDE_LaserBeams_hpp

/****************** Structures ************************************************/
struct LaserBeamParams
{
    double P;          // Laser power [W]
    double I;          // Intensity [W/cm^2]
    double W;          // Pulse Energy [J]
    double lambda;     // Wavelength [m]
    double omega;      // Angular frequency [rad/s]
    double wo;         // Beam spot size at the waist [m]
    double T;          // Pulse duration [s]
    double phi_0;      // Field phase at beam waist [rad]
    double dz0;        // Initial position of the center of the pulse [m]
    double z_Rayleigh; // Rayleigh distance [m]
};

/******************************************************************************/
struct Cartesian_EMfield
{
    // Field components of a general laser beam in cartesian coordinates
    double Ex,Ey,Ez;
    double Bx,By,Bz;
};

/******************************************************************************/
struct RPLB_EMfield
{
    // Field components of a radially polarized laser beam
    double Er,Ez,B_theta;     
};

/******************************************************************************/

#endif // End of file
