//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineEquation.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
#include "BlineEquation.hh"


///////////////////////////////////////////////////////////////////////////
BlineEquation::BlineEquation( G4MagneticField* MagField )
                              : G4Mag_EqRhs( MagField ) 
{backward_direction=false;
 direction=1.;
}


/////////////////////////////////////////////////////////////////////////////
void BlineEquation::EvaluateRhsGivenB( const G4double y[],
			                        const G4double B[3],
				                      G4double dydx[] ) const
{ G4double Bmag=direction*std::sqrt(B[0] *B[0]+ B[1]*B[1]+ B[2]*B[2]);
   dydx[0] = B[0]/Bmag;       
   dydx[1] = B[1]/Bmag;       
   dydx[2] = B[2]/Bmag;
   
   dydx[3]=0. * y[0]; //y[0] is used to remove warning
   dydx[4]=0.;
   dydx[5]=0.; 
   return ;


}


//////////////////////////////////////////////////////////////////////
void BlineEquation::SetBackwardDirectionOfIntegration(G4bool abool)
{backward_direction=abool;
 direction=1.;
 if (backward_direction) direction= -1.;
}






