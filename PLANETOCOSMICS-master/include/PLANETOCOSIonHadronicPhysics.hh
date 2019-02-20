#ifndef PLANETOCOSIonHadronicPhysics_h
#define PLANETOCOSIonHadronicPhysics_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4IonInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSIonHadronicPhysics : public G4VPhysicsConstructor
{
  public: 
    PLANETOCOSIonHadronicPhysics (const G4String& name="Ion");
    virtual ~PLANETOCOSIonHadronicPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle ();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess ();

  private:
   G4String mode;

  protected:
  
  //The Elasticprocess is now removed for alpha,tritium and deuterum 
  //For these particles  it is now defined with G4HadronElasticPhysics
  
  // Elastic Process
   G4LElastic*            theElasticModel;
   

   // Generic Ion physics
   G4HadronElasticProcess theIonElasticProcess;
   G4IonInelasticProcess  theIonInelasticProcess; 
  
   // Deuteron physics
  // G4HadronElasticProcess      theDElasticProcess;
   G4DeuteronInelasticProcess  theDeuteronInelasticProcess;
   G4LEDeuteronInelastic*      fDeuteronLEInelasticModel;

   // Triton physics
   //G4HadronElasticProcess      theTElasticProcess;
   G4TritonInelasticProcess    theTritonInelasticProcess;
   G4LETritonInelastic*        fTritonLEInelasticModel;
  
   // Alpha physics
   //G4HadronElasticProcess      theAElasticProcess;
   G4AlphaInelasticProcess     theAlphaInelasticProcess;
   G4LEAlphaInelastic*         fAlphaLEInelasticModel;

   // He3 physics
   G4HadronElasticProcess      theHe3ElasticProcess;
   G4IonInelasticProcess    theHe3InelasticProcess;
   
  
  

  

};
////////////////////////////////////////////////////////////////////////////////
#endif
