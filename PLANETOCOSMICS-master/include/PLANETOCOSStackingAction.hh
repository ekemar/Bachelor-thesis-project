#ifndef PLANETOCOSStackingAction_h
#define PLANETOCOSStackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "PLANETOCOSSteppingAction.hh"

class G4Track;
class G4ParticleDefinition;
class PLANETOCOSStackingActionMessenger;


class PLANETOCOSStackingAction : public G4UserStackingAction
{
  public:
    PLANETOCOSStackingAction();
    virtual ~PLANETOCOSStackingAction();

  public:
    
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();
    inline void SetUntrackedParticleDic(PartDictionary* aDic)
                                                 {UntrackedParticleDic = aDic;}
    inline void SetStopAllParticles(G4bool aBool){stop_all_particles = aBool;} 					 
    inline void SetBugVerbose(G4int n){bug_verbose =n;}	
  private:
 
    PartDictionary* UntrackedParticleDic; 
    G4bool stop_all_particles;
    G4int bug_verbose;
       
};
  

#endif

