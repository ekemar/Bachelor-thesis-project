#include "PLANETOCOSStackingAction.hh"

#include "PLANETOCOSSteppingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "PLANETOCOSTrackInformation.hh"
#include "PlanetManager.hh"
#include "PLANETOCOSSD.hh"


PLANETOCOSStackingAction::PLANETOCOSStackingAction()
{ bug_verbose =1;}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSStackingAction::~PLANETOCOSStackingAction()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4ClassificationOfNewTrack 
PLANETOCOSStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{  

   
   G4ClassificationOfNewTrack classification = fUrgent;
   
   if (stop_all_particles) return fKill;
  
  //Stop particle for which E^2 != p^2+m^2
  //bug occuring with BinaryCascade for protons around 8.8 GeV in atmosphere 
  //--------------------------------------------
 const G4DynamicParticle*  thePart =aTrack->GetDynamicParticle();
 double E =thePart->GetTotalEnergy();
 double p =thePart->GetMomentum().mag();
 double Ekin = thePart->GetKineticEnergy();
 double m0= thePart->GetDefinition()->GetPDGMass();
 if (std::abs(E-m0-Ekin) >100.*eV ) {
 		G4bool proton_or_neutron =
			(thePart->GetDefinition()->GetParticleName() == "neutron" || thePart->GetDefinition()->GetParticleName() == "proton");
	  	if (proton_or_neutron) {
			PLANETOCOSAnalysisManager::GetInstance()->SwitchOffRegisterResults();    
 			if (bug_verbose >0){
				G4cout<<"The particle will be stopped because E!=Ekin+m0";
				G4cout<<"particle: "<<
		   			thePart->GetDefinition()->GetParticleName()<<std::endl;
					G4cout<<"Creator process:"
		      				<<aTrack->GetCreatorProcess()
		        			->GetProcessName()<<std::endl;
			   
					G4cout<<"E[MeV] Ekin P[MeV/c] m0 E-Ekin-E0: "<<E/MeV<<'\t'
		                                    <<Ekin/MeV<<'\t'
					            <<p/MeV<<'\t'
						    <<m0/MeV<<'\t'
						   <<(E-Ekin-m0)/MeV<<std::endl;
						   
			}
			classification=fKill;
			stop_all_particles = true;
		}
		else {
			if (bug_verbose >0){
				G4cout<<"Warning a particle has been produced with E!=Ekin+m0";
				G4cout<<"particle: "<<
		   			thePart->GetDefinition()->GetParticleName()<<std::endl;
					G4cout<<"Creator process:"
		      				<<aTrack->GetCreatorProcess()
		        			->GetProcessName()<<std::endl;
			   
					G4cout<<"E[MeV] Ekin P[MeV/c] m0 E-Ekin-E0: "<<E/MeV<<'\t'
		                                    <<Ekin/MeV<<'\t'
					            <<p/MeV<<'\t'
						    <<m0/MeV<<'\t'
						   <<(E-Ekin-m0)/MeV<<std::endl;
						   
			}
		}					        
 }
 
  G4String aParticleName=aTrack->GetDefinition()->GetParticleName();

   
  if (aTrack->GetParentID() != 0){
	G4String procName =aTrack->GetCreatorProcess()->GetProcessName();
	G4String volumeName ="";
  	if (aTrack->GetVolume())
  		volumeName = aTrack->GetVolume()->GetName();

	//Register ion produced by Inelastic reaction in the atmosphere 
 	//-------------------------------------------------------------
	if (aParticleName== "alpha" || 
	    aParticleName== "He3" ||
	    aParticleName== "triton" ||
	    aParticleName== "deuteron" || aParticleName.contains("[")) {
	
		if (volumeName.contains("Atmos") &&  !procName.contains("LElast")){
			G4int Z= int(aTrack->GetDefinition()->GetPDGCharge()/eplus); 
        		G4int A=  aTrack->GetDefinition()->GetBaryonNumber(); 
        		G4int N=A-Z;
			G4double Weight = aTrack->GetWeight();
			PLANETOCOSAnalysisManager::GetInstance()
						->DetectNuclideProducedInAtmosphere(N,Z, Weight);
		}
		
	}
	
	
      //kill particle with energy greater than stop energy
      //---------------------------------------------------
	
	if ( UntrackedParticleDic->find(aParticleName) 
                               != UntrackedParticleDic->end() ||
	    (aParticleName.contains("[") && UntrackedParticleDic->find("GenericIon") 
                                  != UntrackedParticleDic->end())){
     		
		G4double emin=(*UntrackedParticleDic)[aParticleName];
		
		
		
      		if (Ekin <= emin) {
			
     			// pass the track to the sensitive detector
       			G4SDManager* SDman   = G4SDManager::GetSDMpointer();
       			PLANETOCOSSD* sd = dynamic_cast<PLANETOCOSSD*>
                               		(SDman->FindSensitiveDetector("atmoSD") );
                	//G4cout<<"Test stack1"<<std::endl;
			//G4cout<<sd<<std::endl;
       			if ( !aParticleName.contains("nu_")){
	   			sd->RegisterEkinOfKilledParticle( aTrack);
	   		}
			//G4cout<<"Test stack2"<<std::endl; 
       		classification=fKill;
		}
		
      	}
	
	
  }
  
  if (classification !=fKill){
  	
//	G4cout<<"Stacking 1"<<std::endl;
	
  	const PLANETOCOSGeometryConstruction* theGeometry = 
     		dynamic_cast< const PLANETOCOSGeometryConstruction* > 
                	       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  	G4Track* modTrack =  const_cast<G4Track*> (aTrack);
	PLANETOCOSTrackInformation* aTrackInformation = 
   				new PLANETOCOSTrackInformation(false);
	modTrack->SetUserInformation(aTrackInformation);			
	
  	if (theGeometry->GetGeometryType() =="SPHERICAL" ){
  	
		G4ThreeVector FirstPosition = aTrack->GetPosition();
   		G4double Theta = FirstPosition.getTheta()/degree; 
   		G4bool WasProducedInCusp = (Theta<20. || Theta >160. ); 
		aTrackInformation->SetWasProducedInTheCusp(WasProducedInCusp);
		
   		
		G4double LowestAltitude=	
			FirstPosition.mag()-PlanetManager::GetInstance()->GetRplanet();
  		
		//set lowest altitude
		if (modTrack->GetTrackID() >1){
		
			const PLANETOCOSSteppingAction* 
			    theSteppingAction = 
		               dynamic_cast<const PLANETOCOSSteppingAction*>
				  (G4RunManager::GetRunManager()->GetUserSteppingAction());
				       
			G4double LastLowestAltitude =theSteppingAction->GetLastLowestAltitude();
			if (LastLowestAltitude < LowestAltitude){
				LowestAltitude = LastLowestAltitude;

			}
		}
		aTrackInformation->SetLowestAltNeededForThisTrack(LowestAltitude);
  	
	}
	
	
//	G4cout<<"Stacking 2"<<std::endl;
  }	
  
  // check if ionisation e- and is produce on boundary
  
  if (aTrack->GetParentID() != 0 &&  aParticleName == "e-"){
  	//G4String aName = aTrack->GetCreatorProcess()->GetProcessName();
  	
	//std::cout<<aTrack->GetVolume()<<std::endl;
	//std::cout<<aTrack->GetNextVolume()<<std::endl;
	
  
  
  }
  
   
  
  
 		     
  return classification;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSStackingAction::NewStage()
{ stackManager->ReClassify();
}
////////////////////////////////////////////////////////////////////////////////
//    
void PLANETOCOSStackingAction::PrepareNewEvent()
{;   
}



