#ifndef MagneticShieldingToolSteppingAction_h
#define MagneticShieldingToolSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>


class MagneticShieldingToolEventAction;
class MagneticShieldingTool;
class PLANETOCOSSteppingAction;

class MagneticShieldingToolSteppingAction : public G4UserSteppingAction
{
  public:
    MagneticShieldingToolSteppingAction(MagneticShieldingTool* aMagneticShieldingTool);
    virtual ~MagneticShieldingToolSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);
    inline void SetUserSteppingAction(PLANETOCOSSteppingAction* anAction){ fUserSteppingAction = anAction;};
    inline void SetBlineMode(G4bool aBool){BlineMode =aBool;}
    
    private:
    MagneticShieldingTool* fMagneticShieldingTool;
    MagneticShieldingToolEventAction* eventAction;
    PLANETOCOSSteppingAction* fUserSteppingAction;
    G4bool BlineMode;
    
    
     
   
};

#endif
