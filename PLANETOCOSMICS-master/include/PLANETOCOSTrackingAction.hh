#ifndef PLANETOCOSTrackingAction_h
#define PLANETOCOSTrackingAction_h 1
#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "globals.hh"

////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSTrackingAction : public G4UserTrackingAction
{
  public:
    PLANETOCOSTrackingAction ();
    ~PLANETOCOSTrackingAction ();
   
    void PreUserTrackingAction (const G4Track* theTrack);
    void PostUserTrackingAction (const G4Track* theTrack);

  private:
    G4bool RegisterLastPoint;
    
  public:
    inline void SetRegisterLastPoint(G4bool abool)
                                      {RegisterLastPoint=abool;}  
};
////////////////////////////////////////////////////////////////////////////////
#endif
