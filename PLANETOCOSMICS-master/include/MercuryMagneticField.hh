#ifndef MercuryMAGNETICFIELD_HH
#define MercuryMAGNETICFIELD_HH 




#include "globals.hh"
#include"G4ios.hh"

#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"
#include "PlanetMagneticField.hh"


class MercuryMagneticFieldMessenger;
class G4ChordFinder;
class G4MagIntegratorStepper;

class MercuryMagneticField : public PlanetMagneticField
{
public:
	 //constructor destructor	       
         MercuryMagneticField() ;
	 
	 ~MercuryMagneticField() ;
	 
	 
	 //Methods that should be provided for specific planet
	 void SetInternalField(G4String aString);
	 void SetExternalField(G4String aString);
	 void SetMagnetopauseModel(G4String aString);
	 void ComputeBfieldParametersAtReferenceTime(); 
	 bool OutsideMagnetosphere( G4ThreeVector pos) const; 
	 bool TSY96MercuryOutsideMagnetosphere(G4ThreeVector pos) const;
       			       
	 								       				
	 
protected:

private:



        // messenger
       MercuryMagneticFieldMessenger* theFieldMessenger;
	
	

			
//private methods       
private:     
        //pointers on the selected Internal and External field model
	G4ThreeVector GetInternalField (G4ThreeVector) const;
	G4ThreeVector GetExternalField (G4ThreeVector) const;
	
	G4ThreeVector (MercuryMagneticField::* PInternalField)(G4ThreeVector)const;
	G4ThreeVector (MercuryMagneticField::* PExternalField)(G4ThreeVector)const;
	
	//pointer on magnetopause model
	bool (MercuryMagneticField::*POutsideMagnetosphere)
	                                             (G4ThreeVector pos) const;
        
	
	
	
		   
} ;

#endif
