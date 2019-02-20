//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineEquation.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines the equation of motion needed to trace magnetic
// feild lines in the simulation. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
#ifndef BlineAndLorentzEQUATIONOFMOTION
#define BlineAndLorentzEQUATIONOFMOTION 
#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"
class BlineEquation : public G4Mag_EqRhs
{
   public:  // with description

     // Constructor and destructor.
     BlineEquation( G4MagneticField* MagField );
    ~BlineEquation() {;}
    
     // Given the value of the magnetic field B, this function 
     // calculates the value of the derivative dydx.

     void EvaluateRhsGivenB( const G4double y[],
			     const G4double B[3],
				   G4double dydx[] ) const;
     
     
				 			   			   		   
     
     //Set methods
     
     void SetBackwardDirectionOfIntegration(G4bool abool);  
     
    
     
    private:
     G4bool backward_direction;
     G4double direction;
     
     
};

#endif 
