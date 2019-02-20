#ifndef MercuryNewData_HH
#define MercuryNewData_HH 
#include "globals.hh"

//mean orbit elements of MercuryNew taken from Simon J.L. et al.,Numerical expressions for
// precession formulae and mean  elements for the Moon and the planets,
//Astronomy and Astrophysics, 282, 663-683(1994), see p675-676
//+trigonometric term see p.681
static const double a_mercury[]={3.870983098e-1, 0.,  0.};
static const double e_mercury[]={2.056317526e-1,  2.040653e-4, -2.8349e-6};
static const double lm_mercury[]={2.5225090552e2, 5.38101628688982e9,  -1.92789e0};
static const double mu_mercury = 6.0236e6;
static const double Omega_mercury[]={4.833089304e1,  -4.51521727e3,  -3.179892e1};
static const double p_mercury[]={7.745611904e1,  5.71911590e3,   -4.83016e0}; //omega+Omega
static const double i_mercury[]={7.00498625e0, -2.1425629e2,   2.8977e-1};
//trigonometric term
static const double pu_mercury[]={69613, 75645L, 88306L, 59899L, 15746L, 71087L, 142173L,  3086L,0L};
static const double Ca_mercury[]={4L,-13L,11L,-9L,-9L,-3L,-1L, 4L,0L};
static const double Sa_mercury[]={-29L, -1L, 9L, 6L,-6L, 5L, 4L, 0L,0L};
static const double qu_mercury[]={3086L,  15746L, 69613L, 59899L, 75645L, 88306L,  12661L,  2658L,  0L,   0L};
static const double Cl_mercury[]={    21L,   -95L, -157L,   41L,   -5L,   42L,   23L,   30L, 0L,0L};
static const double Sl_mercury[]={ -342L,   136L,  -23L,   62L,   66L,  -52L,  -33L,   17L, 0L,0L};
//orientation and rotationof MercuryNew in GEI2000
static const double alpha0_mercury[]={281.01,-0.033};
static const double delta0_mercury[]={61.45,-0.005};
static const double W0_mercury[]={329.68, 6.1385025}; //from Almanach 1992 and
							// Fränz and Harper paper 
//flattening and radius of MercuryNew 
static const double radius_eq_mercury=2439.7*km; // in km
static const double flattening_mercury =0.;

#endif
