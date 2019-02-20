#ifndef EarthData_HH
#define EarthData_HH 
#include "globals.hh"

//mean orbit elements of Earth taken from Simon J.L. et al.,Numerical expressions for
// precession formulae and mean  elements for the Moon and the planets,
//Astronomy and Astrophysics, 282, 663-683(1994), see p675-676
//+trigonometric term see p.681
static const double a_earth[]={1.0000010178,0.};
static const double e_earth[]={1.67086342e-2, -4.203654e-4, -1.26734e-5};
static const double lm_earth[]={1.0046645683e2, 1.29597742283429e9,  -2.04411};
static const double mu_earth = 3.289005e5;
static const double Omega_earth[]={1.7487317577e2,  -8.67927034e3,1.534191e2};
static const double p_earth[]={1.0293734808e2, 1.161235290e4,   5.327577e2}; //omega+Omega
static const double i_earth[]={0.,  4.6997289e2,  -3.35053};
//trigonometric term
static const double pu_earth[]={16002L, 21863L, 32004L, 10931L, 14529L, 16368L,  15318L, 32794L,0L};
static const double Ca_earth[]={64L,   -152L,62L,-8L,32L,   -41L,19L,   -11L,0L};
static const double Sa_earth[]={ -150L,-46L,68L,54L,14L,24L,   -28L,22L,0L};
static const double qu_earth[]={ 10L,  16002L, 21863L, 10931L,  1473L, 32004L,   4387L,73L,  0L,   0L};
static const double Cl_earth[]={ -325L,  -322L,  -79L,  232L,  -52L,   97L,   55L,  -41L, 0L,0L};
static const double Sl_earth[]={ -105L,  -137L,  258L,   35L, -116L,  -88L, -112L,  -80L, 0L,0L};
//orientation and rotationof Earth in GEI2000
static const double alpha0_earth[]={0.0,-0.641};
static const double delta0_earth[]={90.,-0.557};
static const double W0_earth[]={190.16, 360.9856235}; //from Almanach 1992 and
							// Fränz and Harper paper 
//radius of Earth 
static const double radius_eq_earth=6378.137*km; 
static const double radius_pole_earth=6356.752*km;// Ellipsoid WGS84

#endif
