#include "PLANETOCOSTrackingAction.hh"
#include "G4ThreeVector.hh"
#include "PLANETOCOSApplicationScenario.hh"
#include "PlanetMagneticField.hh"
#include"G4TransportationManager.hh"
#include"G4FieldManager.hh"
#include"G4RunManager.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSTrackInformation.hh"
#include "PLANETOCOSSD.hh"
#include "G4SDManager.hh"
#include "MagneticShieldingTool.hh"
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSTrackingAction::PLANETOCOSTrackingAction ()
{RegisterLastPoint=false;}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSTrackingAction::~PLANETOCOSTrackingAction ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSTrackingAction::PreUserTrackingAction (const G4Track* aTrack)
{ /*std::cout<<aTrack->GetVertexPosition()<<std::endl;
  std::cout<<aTrack->GetPosition()<<std::endl;
  std::cout<<aTrack->GetVolume()<<std::endl;
  std::cout<<aTrack->GetNextVolume()<<std::endl;*/
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSTrackingAction::PostUserTrackingAction (const G4Track* aTrack)
{ 
  //case for magnetic shieldimng analysis mode  
  
  if (RegisterLastPoint) {
	G4ThreeVector LastPosition = aTrack->GetPosition();
   	G4double aCharge =std::abs(aTrack->GetDynamicParticle()
                                          ->GetDefinition()->GetPDGCharge());
   	G4ThreeVector LastMomentumOnCharge = aTrack->GetMomentum() / aCharge;
   
   	//check if particle come from outside the magnetosphere
   
   	PlanetMagneticField* theField = (PlanetMagneticField*)
     		G4TransportationManager::GetTransportationManager()->
                              			GetFieldManager()->GetDetectorField();
   	G4double FilterValue =0;
   	if (theField->OutsideMagnetosphere(LastPosition)) FilterValue=1;
   
   
  	// send information to PLANETOCOSApplicationScenario
   
  	MagneticShieldingTool::GetInstance()->RegisterTrackLastPoint(LastPosition,LastMomentumOnCharge,FilterValue);
  }
  
  
  // Pseudo Trapping Analysis if needed
  //////////////////////// 
  if (PLANETOCOSAnalysisManager::GetInstance()
  			->GetPseudoTrappingAnalyser()
					->GetDetectPseudoTrapping()){
  	//G4cout<<"Test11"<<std::endl;
  	G4double weight  = aTrack->GetWeight();
  	G4int PDGCode = aTrack->GetDefinition()->GetPDGEncoding();
  	G4double lifeTime = aTrack->GetLocalTime();
  	G4double stop_energy = aTrack->GetKineticEnergy();
  	G4double start_energy = aTrack->GetVertexKineticEnergy();
  	PLANETOCOSTrackInformation* theTrackInfo = dynamic_cast<PLANETOCOSTrackInformation*>
		            				(aTrack->GetUserInformation()); 
  
  	G4int NbOfEquatorCrossing = theTrackInfo->GetNbOfEquatorCrossing();
 	if ( theTrackInfo->GetWasAlreadyInMagnetosphere() &&
	      ( !theTrackInfo->GetWasProducedInTheCusp() || NbOfEquatorCrossing >0 )){
  		G4int NbOfOutwardCrossing = theTrackInfo->GetNbOfOutwardCrossingOfLastDetector(); 
  		G4int NbOfInwardCrossing = theTrackInfo->GetNbOfInwardCrossingOfLastDetector();
  
  		G4double NbOfPlanetTurn = theTrackInfo->GetNbOfPlanetTurn();
       		if (NbOfPlanetTurn <0) NbOfPlanetTurn = -NbOfPlanetTurn;
  		//G4cout<<NbOfPlanetTurn<<std::endl;
		PLANETOCOSPostTrackHit* aHit =new  PLANETOCOSPostTrackHit(weight, PDGCode,
  		                        			  start_energy,
		                        			  stop_energy,
		                        			  lifeTime,
                                        			  NbOfPlanetTurn,
								  NbOfEquatorCrossing,
		                        			  NbOfOutwardCrossing, 
		                        			  NbOfInwardCrossing);
  
  		PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                       			(G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
  	        sd->GetPostTrackHitCollection()->insert(aHit);
		//G4cout<<"test hit"<<std::endl;
   	}
  }
  
 /* PLANETOCOSFluxDetectionAnalyser* theFluxAnalyser =PLANETOCOSAnalysisManager::GetInstance()->GetFluxDetectionAnalyser();
 if (theFluxAnalyser->GetDetectLowestAltitudeNeeded()) {
  	PLANETOCOSTrackInformation* theTrackInfo = dynamic_cast<PLANETOCOSTrackInformation*>
		            				(aTrack->GetUserInformation()); 
  	//G4cout<<"Test1"<<std::endl;
	if (theTrackInfo) {
		if (theTrackInfo->GetHasBeenAlreadyUpwardDetected()){
			//G4cout<<"Test2"<<std::endl;
			G4int PDGcode = aTrack->GetDefinition()->GetPDGEncoding();
			G4double energy =theTrackInfo->GetDetectionEnergy();
			G4double lowest_altitude_needed = 
					theTrackInfo->GetLowestAltNeededForThisTrack();
*/		
			/*G4cout<<"energy "<<energy/MeV<<std::endl;
			G4cout<<"code "<<PDGcode<<std::endl;
			G4cout<<"lowest_altitude_needed "<<lowest_altitude_needed/km<<std::endl;*/
		
/*			theFluxAnalyser->RegisterLowestAltitudeNeeded(PDGcode, 
							      energy,
							      lowest_altitude_needed);	
	
	
		}
	 	
        //G4cout<<"Test2"<<std::endl;
  	}
 }*/							
  
 
}
