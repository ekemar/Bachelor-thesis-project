#include "MagneticShieldingToolSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "MagneticShieldingTool.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "G4RunManager.hh"
#include "DurationManager.hh"
#include "PLANETOCOSTrackInformation.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PlanetManager.hh"
MagneticShieldingToolSteppingAction::MagneticShieldingToolSteppingAction(MagneticShieldingTool* aMagneticShieldingTool)
{ fMagneticShieldingTool=aMagneticShieldingTool;
}

void MagneticShieldingToolSteppingAction::UserSteppingAction(const G4Step* aStep )
{ 

  DurationManager* theDurationManager = DurationManager::GetInstance();
  if (!theDurationManager->CheckDuration()) {
        PLANETOCOSAnalysisManager::GetInstance()->SwitchOffRegisterResults();
  	G4cout<< "The event will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
  }
 
  
  
  const G4VPhysicalVolume* currentVolume = aStep->GetPostStepPoint()
                					->GetPhysicalVolume();
 
  
  
  G4Track* aTrack=aStep->GetTrack();
  G4double StepLength = aTrack->GetStepLength(); 	
  G4int nStep =aTrack->GetCurrentStepNumber();	
  if (nStep > 10 && StepLength == 0.) {
  	aTrack->SetTrackStatus(fStopAndKill);
  	return;
  }
  
 
  if (BlineMode) {
  	
  	if (currentVolume) {
		G4String VName = currentVolume->GetName();
		G4String PlanetName = PlanetManager::GetInstance()->GetPlanetName();
		if ( VName == PlanetName || VName.contains("Soil"))
					aTrack->SetTrackStatus(fStopAndKill);
			
	
	}
	return;
  
  }
  
  	
  G4ThreeVector pre_position   = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
  if (fMagneticShieldingTool->CheckIfParticleShouldBeKilledAndSetMaxStepLength(position)){
  	aTrack->SetTrackStatus(fStopAndKill);
	return;
  }
  
  PLANETOCOSTrackInformation* aTrackInfo 
  			= dynamic_cast<PLANETOCOSTrackInformation*>
				            (aTrack->GetUserInformation());
  PLANETOCOSGeometryConstruction* theGeometry =
   			dynamic_cast< PLANETOCOSGeometryConstruction* >
   			  	(const_cast< G4VUserDetectorConstruction* >(
				    G4RunManager::GetRunManager()
						->GetUserDetectorConstruction()));
  
   G4String geometry_type =  theGeometry->GetGeometryType();
//  G4cout<<"test2"<<std::endl;
  if (geometry_type=="SPHERICAL" && currentVolume) {
  	
	if (pre_position.z()*position.z()<0) aTrackInfo->AddOneCrossingOfEquator(); 
   	G4ThreeVector vec1,vec2;
   	vec1=G4ThreeVector(pre_position.x(),pre_position.y(),0 );
   	vec2=G4ThreeVector(position.x(),position.y(),0 );
   	G4double dTurn= std::abs(vec1.angle(vec2))/2./3.1415926;
   	G4ThreeVector vec3 = vec1.cross(vec2);
   	if (vec3.z() <0) dTurn=-dTurn;
   	aTrackInfo->AddNbOfPlanetTurn(dTurn);
	
	G4String VName = currentVolume->GetName();
	if ( VName == "MagnetospherePV" || VName.contains("DetectionLayer"))
		     aTrackInfo->SetWasAlreadyInMagnetosphere(true);
					
	if (std::abs(aTrackInfo->GetNbOfPlanetTurn())>= 
	            fUserSteppingAction->GetMaxNumberOfTurnAroundThePlanet()){
        	aTrack->SetTrackStatus(fStopAndKill);
 	}
  }
  	  	
  G4UserLimits* theUserLimits =0;
  if (currentVolume){
  	theUserLimits = currentVolume->GetLogicalVolume()->GetUserLimits();
  }					 
      
 // G4cout<<"Stepping action 4"<<std::endl; 
  if (theUserLimits) {  
  	G4double max_track_length = theUserLimits->GetUserMaxTrackLength(*aTrack);
  
 	 if (aTrack->GetTrackLength() >= max_track_length){
           		aTrack->SetTrackStatus(fStopAndKill);
 			
	  		return;
  	} 
			 
			 
  	//stop the particle if the proper time of the particle is greater than the maximum allowed
  	// proper time of a trajectory
      
  	G4double max_local_time =
                     		theUserLimits->GetUserMaxTime(*aTrack);
 	 if (aTrack->GetLocalTime() >= max_local_time){
           	aTrack->SetTrackStatus(fStopAndKill);
	   	G4cout<<"Stop  at local time ="<<aTrack->GetLocalTime()/s<< "second"<<std::endl;
	   	return;
  
  	}
  }
  	 			 					    			 
}
