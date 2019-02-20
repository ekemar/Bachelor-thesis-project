#include"orbit_functions.hh"
#include "globals.hh"
///////////////////////////////////////////////////////////////////////////////
//
void orbit_func::state_in_orbit_system(double mu, double a, double e, double landa, 
			               double omega_bar,double Omega,
			               double pos[3], double v[3])
{double M= landa-omega_bar-Omega; 
  
 double PI=2.*std::acos(0.); 
  //Mean anomaly we should probably take the modulo of it
 M=std::fmod(M,2.*PI);
 if (M<0) M=2.*PI+M;
 double E= solving_kepler_equation(M, e, 1.e-9); //Eccentric anomaly
 double cos_E=std::cos(E);
 double p=a*(1.-e*e);
 double cos_vu = (e-cos_E)/(e*cos_E-1);
 double sin_vu = std::sqrt(1. - cos_vu*cos_vu);
 if (M> PI) sin_vu= -sin_vu;
 double r = a*(1-e*cos_E);
 pos[0]=r*cos_vu;
 pos[1]=r*sin_vu;
 pos[2]=0.;
 
 double d = std::sqrt(mu/p);
 v[0]=-d*sin_vu;
 v[1]=d*(e+cos_vu);
 v[2]=0.;
  
 
}
////////////////////////////////////////////////////////////////////////////////
//
void orbit_func::state_in_reference_system(double mu, double a, double e, double landa, 
			                    double omega_bar, double i,double Omega,
			                    double domega_bar_dt, double di_dt,double dOmega_dt,
					    double pos[3], double v[3])
{double pos_orb_sys[3], v_orb_sys[3];
 state_in_orbit_system(mu, a, e, landa, omega_bar, Omega,
		       pos_orb_sys,  v_orb_sys);
 
 double E[3][3], dE_dt[3][3]; 
 euler_matrices(Omega, i, omega_bar, dOmega_dt, di_dt, domega_bar_dt,
		E, dE_dt);
 
 //The system S' is the orbital system
 //The system S is the reference system
 // We have V' =E.V and V=transpose(E).V'
 //We want V knowing V'm therefore we take transpose (E) 
 
 for (int i=0;i<3;i++) {
 	pos[i]=0;
	v[i]=0;
 	for (int j=0;j<3;j++){
		pos[i] += E[j][i]*pos_orb_sys[j];  //see remark above
		v[i] += dE_dt[j][i]*pos_orb_sys[j] + E[j][i] * v_orb_sys[j];
	}
 
 } 
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector orbit_func::pos_in_orbit_system(double , double a, double e, double landa, 
		    double omega_bar,double Omega)
{double M= landa-omega_bar-Omega; 
  
 //Mean anomaly we should probably take the modulo of it
 double PI=2.*std::acos(0.);
 M=std::fmod(M,2.*PI);
 if (M<0) M=2.*PI+M;
 double E= solving_kepler_equation(M, e, 1.e-9); //Eccentric anomaly
 double cos_E=std::cos(E);
 //double p=a*(1.-e*e);
 double cos_vu = (e-cos_E)/(e*cos_E-1);
 double sin_vu = std::sqrt(1. - cos_vu*cos_vu);
 if (M> PI) sin_vu= -sin_vu;
 double r = a*(1-e*cos_E);
 return r*G4ThreeVector(cos_vu,sin_vu,0.);
}		    
////////////////////////////////////////////////////////////////////////////////
//
void orbit_func::euler_matrices(double omega, double theta,double  phi,
		    		double domega_dt, double dtheta_dt,double  dphi_dt,	
		    		double E[3][3], double dE_dt[3][3])
{
 
 double cos_o = std::cos(omega);
 double sin_o = std::sin(omega);
 double cos_t = std::cos(theta);
 double sin_t = std::sin(theta);
 double cos_p = std::cos(phi);
 double sin_p = std::sin(phi);
 
 
 G4cout<<"omega :"<<omega<<std::endl;
 G4cout<<"theta :"<<theta<<std::endl;
 G4cout<<"phi :"<<phi<<std::endl;
 E[0][0]=cos_p*cos_o - sin_p * sin_o* cos_t;
 E[0][1]=cos_p*sin_o + sin_p * cos_o* cos_t;
 E[0][2]=sin_p * sin_t;
 
 E[1][0]=-sin_p*cos_o - cos_p * sin_o* cos_t;
 E[1][1]=-sin_p*sin_o + cos_p * cos_o* cos_t;
 E[1][2]=cos_p * sin_t;
 
 E[2][0]=sin_o*sin_t;
 E[2][1]=-cos_o*sin_t;
 E[2][2]=cos_t;
 
 double AE[3][3];
 
 AE[0][0] = E[1][0];
 AE[0][1] = E[1][1];
 AE[0][2] = E[1][2];
 
 AE[1][0] = -E[0][0];
 AE[1][1] = -E[0][1];
 AE[1][2] = -E[0][2]; 
 
 AE[2][0] = E[2][0];
 AE[2][1] = E[2][1];
 AE[2][2] = E[2][2];
 
 double EA[3][3];
 
 EA[0][0] = -E[0][1];
 EA[1][0] = -E[1][1];
 EA[2][0] = -E[2][1];
 
 EA[0][1] = E[0][0];
 EA[1][1] = E[1][0];
 EA[2][1] = E[2][0]; 
 
 EA[0][2] = E[0][2];
 EA[1][2] = E[1][2];
 EA[2][2] = E[2][2];
 
 double B[3][3];
 
 B[0][0]= cos_p*cos_p;
 B[0][1]= cos_p*sin_p;
 B[0][2]= sin_p;
 
 B[1][0]= -cos_p*sin_p;	
 B[1][1]= -sin_p*sin_p;
 B[1][2]= cos_p;
 
 B[2][0]= sin_p;
 B[2][1]= -cos_p;
 B[2][2]=0;
 
 double BE[3][3];
 
 for (int i=0;i<3;i++){
 	for (int j=0;j<3;j++){
 		BE[i][j]=0.;
		for (int k=0;k<3;k++){
			BE[i][j]+=B[i][k]*E[k][j];
		}
	}
 }
 
 for (int i=0;i<3;i++){
 	for (int j=0;j<3;j++){
		dE_dt[i][j] = AE[i][j] * domega_dt + BE[i][j] * dtheta_dt +
						   + EA[i][j] * dphi_dt;
	}
 }
 
 
 
 
 
 
 
}				
////////////////////////////////////////////////////////////////////////////////
//
void orbit_func::euler_rotation(double omega, double theta,double  phi,
		                const double vec_in[3], double vec_out[3])
{double cos_o = std::cos(omega);
 double sin_o = std::sin(omega);
 double cos_t = std::cos(theta);
 double sin_t = std::sin(theta);
 double cos_p = std::cos(phi);
 double sin_p = std::sin(phi);
 double E[3][3];
 E[0][0]=cos_p*cos_o - sin_p * sin_o* cos_t;
 E[0][1]=cos_p*sin_o + sin_p * cos_o* cos_t;
 E[0][2]=sin_p * sin_t;
 
 E[1][0]=-sin_p*cos_o - cos_p * sin_o* cos_t;
 E[1][1]=-sin_p*sin_o + cos_p * cos_o* cos_t;
 E[1][2]=cos_p * sin_t;
 
 E[2][0]=sin_o*sin_t;
 E[2][1]=-cos_o*sin_t;
 E[2][2]=cos_t;
 
 for (int i=0;i<3;i++) {
 	vec_out[i]=0;
 	for (int j=0;j<3;j++){
		vec_out[i]+= E[j][i]*vec_in[j];
	}
}
return; 
 
 
}
////////////////////////////////////////////////////////////////////////////////
//
double orbit_func::solving_kepler_equation(double M, double e,double precision)
{ double E1,E2,E;
  E1=M-e;
  E2=M+e;
  E= M;
  double M1=E-e*std::sin(E);
  while ( (1-E1/E) >precision && M1 !=M) {
	if (M1<M) E1= E;
	else E2 = E;
	E=(E1+E2)/2.;
	M1=E-e*std::sin(E);
	 
  }
  return E;
  
   
}
