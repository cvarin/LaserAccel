#ifndef INC_Constants_hpp
#define INC_Constants_hpp

#define Pi      3.1415926535897932
#define sqrt_Pi 1.7724538509055160
#define sqrt_2  1.4142135623730950

/*******************************************************************************
   Physical constants from Codata, or derived from it.
   For more information, please see http://physics.nist.gov/cuu/Constants/ or 
   http://www.codata.org/ (Both accessed Sept. 13th 2007). See also
   P. J. Mohr and B. N. Taylor, "CODATA recommended values of the fundamental
   physical constants: 1998," Rev. Mod. Phys., Vol. 72, No. 2, 351-495 (2000).
*******************************************************************************/

#define co            2.99792458e8       // Speed of ligth in vacuum [m/s]
#define inv_co        3.335640952e-9     // 1/co [s/m]
#define inv_co_square 1.112650056e-17    // 1/co^2 [s^2/m^2]
#define mu_0          12.566370614e-7    // Magnetic constant in vacuum [N/A^2]
#define eps_0         8.854187817e-12    // Electric constant in vacuum [F/m]	
#define eta_0         376.730313461      // Impedance of vacuum [Ohm]
#define one_4Pieps0	  8.987551787e9	     // Coulomb's law constant [Nm^2/C^2]

#define mp_kg         1.67262158e-27     // Proton mass [kg]
#define mp_MeV        938.271998         // Proton mass energy (mc^2) [MeV]
#define q_p           1.602176462e-19    // Proton charge [C]
#define me_kg         9.10938188e-031    // Electron mass [kg]
#define me_MeV        0.510998902        // Electron mass energy (mc^2) [MeV]
#define q_e           -1.602176462e-19   // Electron charge [C]

#define uq            1.602176462e-19     // Unitary charge [C]

#endif // INC_Constants_hpp


