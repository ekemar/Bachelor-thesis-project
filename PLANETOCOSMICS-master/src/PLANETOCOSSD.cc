#include "PLANETOCOSSD.hh"
//#include "PLANETOCOSAnalysisManager.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "myfunctions.hh"
#include "G4SDManager.hh"
#include "PLANETOCOSSoilSD.hh"

//#include "G4ProcessTable.hh"
class G4Step;

PLANETOCOSSD::PLANETOCOSSD(G4String name, PLANETOCOSGeometryConstruction* det)
:G4VSensitiveDetector(name), fDetector(det)
{ collectionName.insert("primaryCol"); 
  collectionName.insert("detCol"); //equivalent to push_back
  collectionName.insert("edepCol"); //equivalent to push_back
  collectionName.insert("post_trackCol");

 
  DetectEdepvsAltitude = false;
  DetectEdepvsDepth = false;
  primaryHit=0;
  pAltitudes=0;

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSD::~PLANETOCOSSD()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::Initialize(G4HCofThisEvent*HCE)
{ static int HCID = -1;
  PrimaryFluxHitCollection = new PLANETOCOSPrimaryHitsCollection
                                  (SensitiveDetectorName,collectionName[0]);
  DetectorFluxHitCollection = new PLANETOCOSFluxHitsCollection
                                  (SensitiveDetectorName,collectionName[1]);
  
  EdepHitCollection = new PLANETOCOSEdepHitsCollection
                                  (SensitiveDetectorName,collectionName[2]);
  
  PostTrackHitCollection = new PLANETOCOSPostTrackHitsCollection
                                  (SensitiveDetectorName,collectionName[3]);				  
  HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,PrimaryFluxHitCollection);
  
  
  HCID = GetCollectionID(1);
  HCE->AddHitsCollection(HCID,DetectorFluxHitCollection);
 
  
  HCID = GetCollectionID(2);
  HCE->AddHitsCollection(HCID,EdepHitCollection);
  
 
  HCID = GetCollectionID(3);
  HCE->AddHitsCollection(HCID,PostTrackHitCollection); 
  

}
////////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSSD::ProcessHits(G4Step*aStep,G4TouchableHistory* )
{G4int nlayer;
 nlayer= -1;
 G4double edep = aStep->GetTotalEnergyDeposit();
 
 if (edep > 0. && (DetectEdepvsAltitude ||  DetectEdepvsDepth) ){
 	nlayer = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo();
    	nlayer += -nb_magneto_layers;
	PLANETOCOSEdepHit* aHit = new PLANETOCOSEdepHit();
    	aHit->SetEdep(edep);
    	aHit->SetWeight(aStep->GetTrack()->GetWeight());
    	G4ThreeVector pos1 = aStep->GetPreStepPoint()->GetPosition();
    	G4ThreeVector pos2 = aStep->GetPostStepPoint()->GetPosition();
    	G4double alt1,alt2;
    	/*G4int PDGCode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    	aHit->SetPDGCode(PDGCode);*/
    	if ( GeometryType == "SPHERICAL"){
      		alt1 = pos1.mag();
       		alt2 = pos2.mag(); 
      	}
    	else { 
      		alt1 = pos1.z();
       		alt2 = pos2.z(); 
      	}	
      
    	alt1 -= ZSeaLevel;
    	alt2 -= ZSeaLevel;  
    	aHit->SetAltitude1(alt1);
    	aHit->SetAltitude2(alt2);
  
    	if (DetectEdepvsDepth) {
      		G4double alt_up,alt_low,depth_up,depth_low;
		/*G4cout<<pAltitudes<<std::endl;
		G4cout<<pDepths<<std::endl;*/
       		/*G4cout<<"test2 test2 test2 test2"<<std::endl;
      
       	        G4cout<<"alt1 "<<alt1/km<<std::endl;
		G4cout<<"alt2 "<<alt1/km<<std::endl;
       		G4cout<<(*pAltitudes)[nlayer]/km<<std::endl;
       		G4cout<<(*pAltitudes)[nlayer+1]/km<<std::endl;*/
		alt_up = (*pAltitudes)[nlayer];
       		alt_low = (*pAltitudes)[nlayer+1];
       		depth_up = (*pDepths)[nlayer];
       		depth_low = (*pDepths)[nlayer+1];
       		
		
       		G4double depth1,depth2;
       		depth1 = depth_up + (alt1-alt_up)*(depth_low -depth_up) / (alt_low -alt_up);
       		depth2 = depth_up + (alt2-alt_up)*(depth_low -depth_up) / (alt_low -alt_up);
       		aHit->SetDepth1(depth1);
       		aHit->SetDepth2(depth2);
#ifdef  TEST_EDEP  
		G4ParticleDefinition* aDef = aStep->GetTrack()->GetDefinition();
		G4String aParticleName=aDef->GetParticleName();
		G4String aParticleType =aDef->GetParticleType();
		G4int ParticleTypeCode=-1;
		if (aParticleName == "gamma"){
			ParticleTypeCode=0;
		}
		else if (aParticleType == "lepton"){
			if (aParticleName == "e-"){
				ParticleTypeCode=2;
			}
			else if (aParticleName == "e+"){
				ParticleTypeCode=4;
			}
			else if (aParticleName ==("mu-") || aParticleName ==("mu+")){
				ParticleTypeCode=6;
			
			}
			else if (aParticleName ==("tau-") || aParticleName ==("tau+")){
				ParticleTypeCode=8;
			
			}
		}
		else if (aParticleType == "baryon"){
			if (aParticleName == "proton"){
				ParticleTypeCode=10;
			}
			else if (aParticleName == "neutron"){
				ParticleTypeCode=12;
			}
			else ParticleTypeCode=14;
			
		}
		else if (aParticleType == "nucleus"){
			ParticleTypeCode=16;
		}
		else if (aParticleType == "meson"){
			ParticleTypeCode=18;
		}
		aHit->SetParticleTypeCode(ParticleTypeCode);
			
#endif		
      	}
	
      	EdepHitCollection->insert(aHit);
	
   
    
  }
 return true;   
}  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::EndOfEvent(G4HCofThisEvent*)
{if (primaryHit)  PrimaryFluxHitCollection->insert(primaryHit);
 primaryHit = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::clear()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::DrawAll()
{
} 
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::PrintAll()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSD::RegisterEkinOfKilledParticle(const G4Track* aTrack)
{ 


  
  //Check if the particle in an Atmospheric layer
  if (aTrack->GetStep()){
  	const G4VPhysicalVolume* currentVolume = aTrack->GetStep()->GetPostStepPoint()
                                                       ->GetPhysicalVolume();
  	//G4cout<<"Test SD1"<<std::endl;
	if (currentVolume){
		if (!currentVolume->GetName().contains("Atmo")){
			G4SDManager* SDman   = G4SDManager::GetSDMpointer();
       			PLANETOCOSSoilSD* sd = dynamic_cast<PLANETOCOSSoilSD*>
                               			(SDman->FindSensitiveDetector("soilSD") );
			//G4cout<<"Test SD2"<<std::endl;
			return sd->RegisterEkinOfKilledParticle(aTrack);
			//G4cout<<"Test SD3"<<std::endl;				
		
		} 
		
	}
	//G4cout<<"Test SD4"<<std::endl;		 
  }
  else {
  	const G4LogicalVolume* theLogicalVolume = aTrack->GetLogicalVolumeAtVertex();
	G4VPhysicalVolume* thePhysicalVolume =aTrack->GetVolume(); 
	G4String volName;
	
	if (theLogicalVolume) volName=theLogicalVolume->GetName();
	else if (thePhysicalVolume) volName=thePhysicalVolume->GetName();
	//G4cout<<"Test SD5"<<std::endl;
	
  	if (theLogicalVolume || thePhysicalVolume){
		if (!volName.contains("Atmo")) {
			G4SDManager* SDman   = G4SDManager::GetSDMpointer();
       			PLANETOCOSSoilSD* sd = dynamic_cast<PLANETOCOSSoilSD*>
                               			(SDman->FindSensitiveDetector("soilSD") );
			//G4cout<<"Test SD6"<<std::endl;
			return sd->RegisterEkinOfKilledParticle(aTrack);
			//G4cout<<"Test SD7"<<std::endl;
		}
  	}
	//G4cout<<"Test SD8"<<std::endl;
	
  		
  }
  
  if (!pAltitudes) return;
  
  G4double ekin = aTrack->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetPosition();
  G4double alt;
  if ( GeometryType == "SPHERICAL") alt = pos.mag();
  else alt = pos.z();
  alt -= ZSeaLevel;
 
  G4bool  out_of_range;
  G4int nlayer  = myfunc::locate(*pAltitudes, alt,out_of_range);
  
  
  
  if (ekin > 0. && (DetectEdepvsAltitude ||  DetectEdepvsDepth) && !out_of_range ){
  
	PLANETOCOSEdepHit* aHit = new PLANETOCOSEdepHit();
     	aHit->SetEdep(ekin);
    	aHit->SetWeight(aTrack->GetWeight());
     	G4ThreeVector pos = aTrack->GetPosition();
     	aHit->SetAltitude1(alt);
     	aHit->SetAltitude2(alt);
     
     	if (DetectEdepvsDepth){
      		G4double alt_up,alt_low,depth_up,depth_low;
       		/*G4cout<<"test2 test2 test2 test2"<<std::endl;
      
       	        G4cout<<"alt "<<alt/km<<std::endl;
       		G4cout<<(*pAltitudes)[nlayer]/km<<std::endl;
       		G4cout<<(*pAltitudes)[nlayer+1]/km<<std::endl;*/
       
       		alt_up = (*pAltitudes)[nlayer];
       		alt_low = (*pAltitudes)[nlayer+1];
      		depth_up = (*pDepths)[nlayer];
       		depth_low = (*pDepths)[nlayer+1];
		
       
       		G4double depth;
       		depth = depth_up +(alt-alt_up)*(depth_low -depth_up) / (alt_low -alt_up);
       		aHit->SetDepth1(depth);
       		aHit->SetDepth2(depth);

#ifdef  TEST_EDEP  
		G4ParticleDefinition* aDef = aTrack->GetDefinition();
		G4String aParticleName = aDef->GetParticleName();
		G4String aParticleType = aDef->GetParticleType();
		G4int ParticleTypeCode=-1;
		if (aParticleName == "gamma"){
			ParticleTypeCode=1;
		}
		else if (aParticleType == "lepton"){
			if (aParticleName == "e-"){
				ParticleTypeCode=3;
			}
			else if (aParticleName == "e+"){
				ParticleTypeCode=5;
			}
			else if (aParticleName ==("mu-") || aParticleName ==("mu+")){
				ParticleTypeCode=7;
			
			}
			else if (aParticleName ==("tau-") || aParticleName ==("tau+")){
				ParticleTypeCode=9;
			
			}
		}
		else if (aParticleType == "baryon"){
			if (aParticleName == "proton"){
				ParticleTypeCode=11;
			}
			else if (aParticleName == "neutron"){
				ParticleTypeCode=13;
			}
			else ParticleTypeCode=15;
			
		}
		else if (aParticleType == "nucleus"){
			ParticleTypeCode=17;
		}
		else if (aParticleType == "meson"){
			ParticleTypeCode=19;
		}
		aHit->SetParticleTypeCode(ParticleTypeCode);
			
#endif
      	}
     
     	EdepHitCollection->insert(aHit);
  }
 
 
  
}					
