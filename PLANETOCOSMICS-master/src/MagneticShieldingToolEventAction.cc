#include "MagneticShieldingToolEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "MagneticShieldingTool.hh"
#include "DurationManager.hh"
#include "G4RunManager.hh"



MagneticShieldingToolEventAction::MagneticShieldingToolEventAction(MagneticShieldingTool* aMagneticShieldingTool)
{fMagneticShieldingTool=aMagneticShieldingTool;
}

MagneticShieldingToolEventAction::~MagneticShieldingToolEventAction()
{;
}

void MagneticShieldingToolEventAction::BeginOfEventAction(const G4Event* evt)
{ DurationManager* theDurationManager = DurationManager::GetInstance();
  if (!theDurationManager->CheckDurationAtBeginOfEvent()) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The run will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
	
  	
  }

 fMagneticShieldingTool->ResetMaxStepLength();


 G4ThreeVector pos = evt->GetPrimaryVertex()->GetPosition();

 fMagneticShieldingTool->CheckIfParticleShouldBeKilledAndSetMaxStepLength(pos);

 fUserEventAction->BeginOfEventAction(evt);

 
}

void MagneticShieldingToolEventAction::EndOfEventAction(const G4Event* evt)
{ fUserEventAction->EndOfEventAction(evt);
} 



