#include "PLANETOCOSSoilSD.hh"
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

//#include "G4ProcessTable.hh"
class G4Step;

PLANETOCOSSoilSD::PLANETOCOSSoilSD(G4String name, PLANETOCOSGeometryConstruction* det)
:G4VSensitiveDetector(name), fDetector(det)
{ 
  collectionName.insert("edepSoilCol"); //equivalent to push_back
  DetectEdepVsThickness = false;
  DetectEdepVsDepth = false;
  pDepthsInLengthUnit =0;


}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSoilSD::~PLANETOCOSSoilSD()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::Initialize(G4HCofThisEvent*HCE)
{ static int HCID = -1;
  
  SoilEdepHitCollection = new PLANETOCOSEdepHitsCollection
                                  (SensitiveDetectorName,collectionName[0]);
  
  HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,SoilEdepHitCollection);
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSSoilSD::ProcessHits(G4Step*aStep,G4TouchableHistory* )
{G4int nlayer;
 nlayer= -1;
 G4double edep = aStep->GetTotalEnergyDeposit();
 
 if (edep > 0. && (DetectEdepVsThickness ||  DetectEdepVsDepth) ){
 	nlayer = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo();
	//G4cout<<"Nlayer "<<nlayer<<std::endl;
    	nlayer += -nb_layers_above_soil;
	PLANETOCOSEdepHit* aHit = new PLANETOCOSEdepHit();
    	aHit->SetEdep(edep);
    	aHit->SetWeight(aStep->GetTrack()->GetWeight());
    	G4ThreeVector pos1 = aStep->GetPreStepPoint()->GetPosition();
    	G4ThreeVector pos2 = aStep->GetPostStepPoint()->GetPosition();
    	G4double h1,h2;
    	
    	if ( GeometryType == "SPHERICAL"){
      		h1 = pos1.mag();
       		h2 = pos2.mag(); 
      	}
    	else { 
      		h1 = pos1.z();
       		h2 = pos2.z(); 
      	}	
      
    	h1 -= ZGroundLevel;
    	h2 -= ZGroundLevel;

	//here h represents the thickness of soil above the given position
	h1 =-h1;
	h2=-h2;  
	//G4cout<<h1/m<<'\t'<<h2/m<<std::endl;
    	aHit->SetAltitude1(h1);
    	aHit->SetAltitude2(h2);
  
    	if (DetectEdepVsDepth) {
      		G4double h_up,h_low,depth_up,depth_low;
		
		h_up = (*pDepthsInLengthUnit)[nlayer];
       		h_low = (*pDepthsInLengthUnit)[nlayer+1];
       		depth_up = (*pDepthsInDepthUnit)[nlayer];
       		depth_low = (*pDepthsInDepthUnit)[nlayer+1];
		//G4cout<<h_up/m<<'\t'<<h_low/m<<std::endl;
       		//G4cout<<depth_up*cm2/g<<'\t'<<depth_low*cm2/g<<std::endl;
		
       		G4double depth1,depth2;
       		depth1 = depth_up + (h1-h_up)*(depth_low -depth_up) / (h_low -h_up);
       		depth2 = depth_up + (h2-h_up)*(depth_low -depth_up) / (h_low -h_up);
       		aHit->SetDepth1(depth1);
       		aHit->SetDepth2(depth2);
		//G4cout<<depth1*cm2/g<<'\t'<<depth2*cm2/g<<std::endl;
		
	
      	}
	
      	SoilEdepHitCollection->insert(aHit);
	
   
    
  }
 return true;   
}  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::EndOfEvent(G4HCofThisEvent*)
{;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::clear()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::DrawAll()
{
} 
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::PrintAll()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSoilSD::RegisterEkinOfKilledParticle(const G4Track* aTrack)
{ 
  //Check if the particle in an Atmospheric layer
  if (aTrack->GetStep()){
  	const G4VPhysicalVolume* currentVolume = aTrack->GetStep()->GetPostStepPoint()
                                                       ->GetPhysicalVolume();
  	if (currentVolume){
		if (!currentVolume->GetName().contains("Soil")) return;
	}		 
  }
  else {
  	const G4LogicalVolume* theLogicalVolume = aTrack->GetLogicalVolumeAtVertex();
	G4VPhysicalVolume* thePhysicalVolume =aTrack->GetVolume(); 
	G4String volName;
	
	if (theLogicalVolume) volName=theLogicalVolume->GetName();
	else if (thePhysicalVolume) volName=thePhysicalVolume->GetName();
	
  	if (theLogicalVolume || thePhysicalVolume){
  	
		if (!volName.contains("Soil")) return;
  	}
	
  		
  }
  
  if (!pDepthsInLengthUnit) return;

  G4double ekin = aTrack->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetPosition();
  G4double h;
  if ( GeometryType == "SPHERICAL") h = pos.mag();
  else h = pos.z();
  h -= ZGroundLevel;
  h = -h;
 
  G4bool  out_of_range;
  G4int nlayer  = myfunc::locate(*pDepthsInLengthUnit, h,out_of_range);
  

  if (ekin > 0. && (DetectEdepVsThickness ||  DetectEdepVsDepth) && !out_of_range ){
   	PLANETOCOSEdepHit* aHit = new PLANETOCOSEdepHit();
     	aHit->SetEdep(ekin);
    	aHit->SetWeight(aTrack->GetWeight());
     	G4ThreeVector pos = aTrack->GetPosition();
     	aHit->SetAltitude1(h);
     	aHit->SetAltitude2(h);
     
     	if (DetectEdepVsDepth){
      		G4double h_up,h_low,depth_up,depth_low;
       		h_up = (*pDepthsInLengthUnit)[nlayer];
       		h_low = (*pDepthsInLengthUnit)[nlayer+1];
      		depth_up = (*pDepthsInDepthUnit)[nlayer];
       		depth_low = (*pDepthsInDepthUnit)[nlayer+1];
		
       
       		G4double depth;
       		depth = depth_up +(h-h_up)*(depth_low -depth_up) / (h_low -h_up);
       		aHit->SetDepth1(depth);
       		aHit->SetDepth2(depth);


      	}
     
     	SoilEdepHitCollection->insert(aHit);
  }
 
 
  
}					
