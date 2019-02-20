#ifndef PLANETOCOSElectroMagneticPhysics_h
#define PLANETOCOSElectroMagneticPhysics_h 1
/////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VPhysicsConstructor.hh"

#include "G4MultipleScattering.hh"
class G4ElectroNuclearBuilder;
class G4MuNuclearInteraction;
class G4SynchrotronRadiation;

////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSElectroMagneticPhysics : public G4VPhysicsConstructor
{
  public: 
    PLANETOCOSElectroMagneticPhysics (const G4String& name , 
    				   const bool consider_emnuc,
				   const bool consider_munuc,
				   const bool consider_sync);
    virtual ~PLANETOCOSElectroMagneticPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess(); 
    
    inline void SetSteppingAlgorithmParameters(bool aBool, double fac){
    	SteppingAlgorithmMsc=aBool;
	facrange=fac;
    	
    }
    
    
    
    

  private:
  
   G4String mode;
   G4bool ConsiderEMNucPhysics;
   G4bool ConsiderMuonNucPhysics;
   G4bool ConsiderSyncPhysics;
   
   G4ElectroNuclearBuilder* theEMNucPhysics;
   G4MuNuclearInteraction* theMuMinusNucInt;
   G4MuNuclearInteraction* theMuPlusNucInt;
   G4SynchrotronRadiation* theElectronSyncRad;
   G4SynchrotronRadiation* thePositronSyncRad;
   
   
   // for test of stepping bounddary algorithm for multiple scattering
   G4bool SteppingAlgorithmMsc;
   G4double facrange;

  protected:
    
   

};
////////////////////////////////////////////////////////////////////////////////
#endif
