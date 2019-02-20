////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSGeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh" 
#include "G4StepLimiter.hh" 
// Particles
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh" 
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSGeneralPhysics::PLANETOCOSGeneralPhysics (const G4String& name)
  :G4VPhysicsConstructor(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSGeneralPhysics::~PLANETOCOSGeneralPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeneralPhysics::ConstructParticle ()
{
  
  
  //light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();   
  
 //boson
   G4BosonConstructor p1Constructor;
   p1Constructor.ConstructParticle();  
   
  //Lepton
   G4LeptonConstructor p2Constructor;
   p2Constructor.ConstructParticle(); 
   
  //Meson
   G4MesonConstructor p3Constructor;
   p3Constructor.ConstructParticle();   
   
   //Baryon
   G4BaryonConstructor p4Constructor;
   p4Constructor.ConstructParticle();  
   
   //ShortLived
   G4ShortLivedConstructor p5Constructor;
   p5Constructor.ConstructParticle(); 
   
   
   
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeneralPhysics::ConstructProcess()
{
  //StepLimiter process since version 7.0
  theParticleIterator->reset();
  G4StepLimiter* theStepLimiterProcess = new G4StepLimiter();
  while( (*theParticleIterator)() ){
    	G4ParticleDefinition* particle = theParticleIterator->value();
	G4ProcessManager* pmanager = particle->GetProcessManager();
    	G4String particleName = particle->GetParticleName();  		
	pmanager->AddDiscreteProcess(theStepLimiterProcess);	
  }	

  // Add Decay Process
  theParticleIterator->reset();
  while((*theParticleIterator)()){
    	G4ParticleDefinition* particle = theParticleIterator->value();
    	G4ProcessManager* pmanager = particle->GetProcessManager();
    	if (fDecayProcess.IsApplicable(*particle)) { 
      		pmanager ->AddProcess(&fDecayProcess);
      		// set ordering for PostStepDoIt and AtRestDoIt
      		pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      		pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
  	}
  }
}
////////////////////////////////////////////////////////////////////////////////

