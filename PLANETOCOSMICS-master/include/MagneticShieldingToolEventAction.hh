#ifndef MagneticShieldingToolEVENTACTION_HH
#define MagneticShieldingToolEVENTACTION_HH 
#include "G4UserEventAction.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

class G4Event;
class G4Polyline;
class G4Polymarker;
class MagneticShieldingTool;


class MagneticShieldingToolEventAction : public G4UserEventAction
{
  public:

   MagneticShieldingToolEventAction(MagneticShieldingTool* aMagneticShieldingTool);
   ~MagneticShieldingToolEventAction();

  public:
   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);
   inline void SetUserEventAction(G4UserEventAction* aUserEventAction)
   					{fUserEventAction = aUserEventAction;}
  private:
  
   MagneticShieldingTool* fMagneticShieldingTool;
   G4UserEventAction* fUserEventAction;
};

#endif

    
