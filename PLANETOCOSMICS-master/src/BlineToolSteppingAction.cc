//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolSteppingAction.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
#include "BlineToolSteppingAction.hh"
#include "BlineToolEventAction.hh"
#include "G4SteppingManager.hh"
#include "BlineTool.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"

BlineToolSteppingAction::BlineToolSteppingAction(BlineTool* aBlineTool)
{ fBlineTool=aBlineTool;}

void BlineToolSteppingAction::UserSteppingAction(const G4Step* aStep )
{ 
  const G4VPhysicalVolume* currentVolume = aStep->GetPostStepPoint()
                ->GetPhysicalVolume();
  if (currentVolume){		
   	const G4String name = currentVolume->GetName();
   	for (unsigned int i=0;i<StopVolumeVector.size();i++){
  		if (name == StopVolumeVector[i]) {
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return;
		}	
	}
  }
  return;
    
   
}
////////////////////////////////////////////////////////////////////////////////
//
void BlineToolSteppingAction::AddStopVolume(G4String aVolumeName)
{ //check if it is not already in the list of volume where the field line tracing
  // should be stopped
  
  for (unsigned int i=0;i<StopVolumeVector.size();i++){
  	if (aVolumeName == StopVolumeVector[i]) {
		G4cout<<"This volume was already selected"<<std::endl;
		return;
	}
  }
  // if the volume exists its name is added into the StopVolumeVector object
  StopVolumeVector.push_back(aVolumeName);
  return;
 
 
 
}
////////////////////////////////////////////////////////////////////////////////
//
void BlineToolSteppingAction::RemoveStopVolume(G4String aVolumeName)
{ 
  
  for (unsigned int i=0;i<StopVolumeVector.size();i++){
  	if (aVolumeName == StopVolumeVector[i]) {
		StopVolumeVector.erase(StopVolumeVector.begin()+i);
		return;
	}
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void BlineToolSteppingAction::ClearStopVolumeVector()
{StopVolumeVector.clear();
}
