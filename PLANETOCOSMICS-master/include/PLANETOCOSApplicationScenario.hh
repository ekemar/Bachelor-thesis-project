#ifndef PLANETOCOSApplicationScenario_h
#define PLANETOCOSApplicationScenario_h 1

#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4UserRunAction.hh"
#include <vector>
#include "G4Timer.hh"

namespace CLHEP {}
using namespace CLHEP;

class PLANETOCOSScenarioMessenger;
class PLANETOCOSPrimaryGeneratorAction;


class PLANETOCOSApplicationScenario : public G4UserRunAction 
{ public:
  
   PLANETOCOSApplicationScenario();
   virtual ~PLANETOCOSApplicationScenario();
  
  public:
   
   virtual void BeginOfRunAction(const G4Run* aRun);
   virtual void EndOfRunAction(const G4Run* aRun);
   void SetRandomSeed();
   inline void SetRandomSeedNeeded(G4bool aVal) {RandomSeedNeeded=aVal;}
   inline G4bool GetRandomSeedNeeded() {return RandomSeedNeeded;}
  
	       
  
 private:
   G4bool RandomSeedNeeded;
   PLANETOCOSScenarioMessenger* theMessenger;
   
  // run duration 
   bool TotalDurationReached; 
    
   
  
 
    
  
  
};












#endif
