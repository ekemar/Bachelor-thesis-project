//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolSteppingAction.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines the stteping actions used by the BlineTool
//    when tracing magnetic field lines 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifndef BlineToolSteppingAction_h
#define BlineToolSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>


class BlineToolEventAction;
class BlineTool;

class BlineToolSteppingAction : public G4UserSteppingAction
{
  public:
    BlineToolSteppingAction(BlineTool* aBlineTool);
    virtual ~BlineToolSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);
    void AddStopVolume(G4String aVolumeName);
    void RemoveStopVolume(G4String aVolumeName);
    void ClearStopVolumeVector();
    
  private:
    BlineTool* fBlineTool;
    std::vector<G4String> StopVolumeVector; 
   
};

#endif
