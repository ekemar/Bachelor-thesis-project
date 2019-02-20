#ifndef MarsMAGNETICFIELD_HH
#define MarsMAGNETICFIELD_HH 




#include "globals.hh"
#include"G4ios.hh"

#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"
#include "PlanetMagneticField.hh"


class MarsMagneticFieldMessenger;
class G4ChordFinder;
class G4MagIntegratorStepper;

class MarsMagneticField : public PlanetMagneticField
{
public:
	 //constructor destructor	       
         MarsMagneticField() ;
	 
	 ~MarsMagneticField() ;
	 
	 
	 //Methods that should be provided for specific planet
	 void SetInternalField(G4String aString);
	 void SetExternalField(G4String aString);
	 void SetMagnetopauseModel(G4String aString);
	 void ComputeBfieldParametersAtReferenceTime(); 
	 bool OutsideMagnetosphere( G4ThreeVector pos) const; 
	     
         inline void SetPuruckerCutR(G4double r)
	 			      {Purucker_cut_r2=r*r;} 
	 			       
				       
	 								       				
	 
protected:

private:



        // messenger
        MarsMagneticFieldMessenger* theFieldMessenger;

    
	//Dipole momentum for radial equivalent dipole model of Purucker
	//
	
	std::vector<double > Purucker_dip_M;
	std::vector<G4ThreeVector > Purucker_dip_Pos;
	G4double Purucker_cut_r2;
	
	
	

			
//private methods       
private:     
        //pointers on the selected Internal and External field model
	G4ThreeVector GetInternalField (G4ThreeVector) const;
	G4ThreeVector GetExternalField (G4ThreeVector) const;
	
	G4ThreeVector (MarsMagneticField::* PInternalField)(G4ThreeVector)const;
	G4ThreeVector (MarsMagneticField::* PExternalField)(G4ThreeVector)const;
	
	//pointer on magnetopause model
	bool (MarsMagneticField::*POutsideMagnetosphere)
	                                             (G4ThreeVector pos) const;
        
	
	//Purucker model
	G4ThreeVector  GetPURUCKER(G4ThreeVector pos) const;
	G4ThreeVector  GetCAIN(G4ThreeVector pos) const;
	void ReadPuruckerData();  
	
	
	
		   
} ;

#endif
