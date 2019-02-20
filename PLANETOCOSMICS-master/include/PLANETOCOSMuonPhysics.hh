#ifndef PLANETOCOSMuonPhysics_h
#define PLANETOCOSMuonPhysics_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"

#include "G4MuonMinusCaptureAtRest.hh"
////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSMuonPhysics : public G4VPhysicsConstructor
{
  public: 
    PLANETOCOSMuonPhysics(const G4String& name="muon");
    virtual ~PLANETOCOSMuonPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
  
 

};
////////////////////////////////////////////////////////////////////////////////
#endif
