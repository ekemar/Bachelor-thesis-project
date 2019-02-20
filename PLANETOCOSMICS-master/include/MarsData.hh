#ifndef MarsData_HH
#define MarsData_HH 
#include "globals.hh"

//mean orbit elements of Mars taken from Simon J.L. et al.,Numerical expressions for
// precession formulae and mean  elements for the Moon and the planets,
//Astronomy and Astrophysics, 282, 663-683(1994), see p675-676
//+trigonometric term see p.681
static const double a_mars[]={1.5236793419e0,  3.e-10};
static const double e_mars[]={9.34006477e-2,  9.048438e-4, -8.0641e-6};
static const double lm_mars[]={3.5543299958e2,  6.8905077493988e8,   9.4264e-1};
static const double mu_mars = 3.09871e6;
static const double Omega_mars[]={4.955809321e1, -1.062090088e4, -2.3057416e2};
static const double p_mars[]={3.3606023395e2, 1.598045908e4,-6.2328e1}; //omega+Omega
static const double i_mars[]={1.84972648e0, -2.9331722e2,  -8.11830e0};
//trigonometric term
static const double pu_mars[]={6345L,   7818L, 15636L,  7077L,  8184L, 14163L,   1107L,  4872L,0L};
static const double Ca_mars[]={124L,621L,  -145L,   208L,54L,   -57L,30L,15L,0L};
static const double Sa_mars[]={-621L,532L,  -694L,   -20L,   192L,   -94L,71L,   -73L,0L};
static const double qu_mars[]={10L,   6345L,  7818L,  1107L, 15636L,  7077L,   8184L,   532L, 10L,   0L};
static const double Cl_mars[]={2268L,  -979L,  802L,  602L, -668L,  -33L,  345L,  201L,   -55L,0L};
static const double Sl_mars[]={854L,  -205L, -936L, -240L,  140L, -341L,  -97L, -232L,   536L,0L};
//orientation and rotationof Mars in GEI2000
/*Old values from /from Almanach 1992 and
static const double alpha0_mars[]={317.681,-0.108};
static const double delta0_mars[]={52.886,-0.061};
static const double W0_mars[]={176.901, 350.8919830}; // from pck00007.pc
//static const double W0_mars[]={176.868, 350.8919830}; //from Almanach 1992 and
							// Fränz and Harper paper 

*/
//New values from IAU2000
static const double alpha0_mars[]={317.68143,-0.1061};
static const double delta0_mars[]={52.88650,-0.0609};
static const double W0_mars[]={176.630, 350.89198226}; 
//flattening and radius of Mars 
static const double radius_eq_mars=3397.*km; // in km
static const double flattening_mars =6.4763e-3;

#endif
