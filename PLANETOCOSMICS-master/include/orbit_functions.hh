#ifndef MYFUNC_h
#define MYFUNC_h 1
#include <fstream>
#include <complex>
#include <vector>
#include "globals.hh"
#include "G4ThreeVector.hh"



namespace orbit_func{
void state_in_orbit_system(double mu, double a, double e, double landa, 
			    double omega_bar,double Omega,
			    double pos[3], double v[3]);
			    
void state_in_reference_system(double mu, double a, double e, double landa, 
			    double omega_bar, double i,double Omega,
			    double domega_bar_dt, double di_dt,double dOmega_dt,
			    double pos[3], double v[3]);
			    
G4ThreeVector pos_in_orbit_system(double mu, double a, double e, double landa, 
			 double omega_bar,double Omega);			    

//compute the Euler transformation of  a vector vec_in into a vector vec_out
//the euler rotation is defined by the Euler angles omega, theta and phi
//the euler rotation is considered here according to the definition in Appendix 1 of 
// M. Fraenz, and  D. Harper,Heliospheric coordinate systems, Corrected version
//of Plan. Space Science, 50, 217, 2002    
void euler_matrices(double omega, double theta,double  phi,
		    double domega_dt, double dtheta_dt,double  dphi_dt,	
		    double E[3][3], double dE_dt[3][3]);
void euler_rotation(double omega, double theta,double  phi,
		    const double vec_in[3], double vec_out[3]);
		    

//find the   solution of the kepler equation	
// M = E- e*sin(E)
//only 0<e<1 (ellipse) is considered
// in this case the function E-e*sin(E) is monotically increasing 
// indeded the derivate of this function is 1-e*cos(E) >0 if 0<e<1	    	
double solving_kepler_equation(double M, double e, double precision);		    		    			    




}
#endif
