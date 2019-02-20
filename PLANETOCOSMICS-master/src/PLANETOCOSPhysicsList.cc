////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSPhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "PLANETOCOSGeometryConstruction.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ios.hh"
#include <iomanip>

#include "PLANETOCOSPhysicsListMessenger.hh"


//#include "GeneralPhysics.hh"
//#include "EMPhysics.hh"

#include "PLANETOCOSGeneralPhysics.hh"
#include "PLANETOCOSElectroMagneticPhysics.hh"
#include "PLANETOCOSMuonPhysics.hh"
#include "PLANETOCOSHadronicPhysics.hh"
#include "PLANETOCOSIonHadronicPhysics.hh"
#include "G4RunManager.hh"
#include "G4MiscLHEPBuilder.hh"
#include "G4StoppingHadronBuilder.hh" 
//#include "G4EmStandardPhysics.hh"
//#include "G4HadronQEDBuilder.hh"



////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPhysicsList::PLANETOCOSPhysicsList(): PLANETOCOSVModularPhysicsList()
{ 
  physListMessenger = new PLANETOCOSPhysicsListMessenger(this);
//
//
// The default cut value is set to 1.0mm.
//
  defaultCutValue  = 10. *m;
  cutForGamma      = defaultCutValue;
  cutForElectron   = defaultCutValue;
  cutForPositron   = defaultCutValue;


  
 
 
  //Default physics list
  HadronicPhysicsType="QGSP_BIC_HP";
  theHadronicPhysics = new PLANETOCOSHadronicPhysics(HadronicPhysicsType);
  theGeneralPhysics = new PLANETOCOSGeneralPhysics("general");
  RegisterPhysics(theGeneralPhysics);
  LightIonHadronicPhysicsType="BIC"; 
  EMPhysicsType="STANDARD";
  ConsiderElectroNucPhysics =false;
  ConsiderMuonNucPhysics =false;
  ConsiderSyncPhysics= false;
  
  //Msc stepping algorithm
  //----------------------
  
  SteppingAlgorithmMsc =true;
  facrange =0.02;
  
  SetVerboseLevel(0);
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPhysicsList::~PLANETOCOSPhysicsList()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPhysicsList::BuildList()
{
  
  

  
  if (EMPhysicsType != "NONE"){
  	
	
	//Muon physics
  	//RegisterPhysics(new PLANETOCOSMuonPhysics("Muon"));  
	
	//EM physics
	PLANETOCOSElectroMagneticPhysics* theEMPhysics;
  	if (EMPhysicsType == "LOWENERGY") {
		theEMPhysics = new PLANETOCOSElectroMagneticPhysics("LowEnergyEM",ConsiderElectroNucPhysics,
										  ConsiderMuonNucPhysics,
										  ConsiderSyncPhysics);
  	}
	else {
  		theEMPhysics = new PLANETOCOSElectroMagneticPhysics("StandardEM",ConsiderElectroNucPhysics,
										  ConsiderMuonNucPhysics,
										  ConsiderSyncPhysics);
  	}
	/*G4cout<<"Stepping algorithm "<<SteppingAlgorithmMsc<<'\t'<<facrange<<std::endl;
	theEMPhysics->SetSteppingAlgorithmParameters(SteppingAlgorithmMsc,facrange);*/
	RegisterPhysics(theEMPhysics);
  
  	
 	//RegisterPhysics(new PLANETOCOSIonHadronicPhysics(LightIonHadronicPhysicsType));
  
  	//hadronic physics
  	if (HadronicPhysicsType != "NOHADRONIC"){
	        theHadronicPhysics->SetHadronicPhysicsType(HadronicPhysicsType);
  		RegisterPhysics(theHadronicPhysics);
	
	//Light ion hadronic physics
	
		RegisterPhysics(new PLANETOCOSIonHadronicPhysics(LightIonHadronicPhysicsType));
  	}
 }	
  
  
  
//General physics already regsitered in constructor for the definition of particles (Chamnges needed for Geant4.8.0)
  

  if (verboseLevel > 1) {
   	G4cout <<"The registered physics modules are: " <<G4endl;
    	G4PhysConstVector::iterator itr;
    	for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
      	G4cout <<"       " <<(*itr)->GetPhysicsName() <<G4endl;
    	}
    	theParticleIterator->reset();
    	while ((*theParticleIterator)()) {
      		G4ParticleDefinition* particle = theParticleIterator->value();
      		G4cout <<particle->GetParticleName() <<G4endl;
      		G4ProcessManager* pmanager = particle->GetProcessManager();
      		G4ProcessVector* alist     = pmanager->GetProcessList();
      		G4cout <<"Number of process = " <<alist->size() <<G4endl;
//      for (G4int i = 0; i < alist->size(); i++) {
//        G4cout << "  " << (alist[i])->GetProcessName() <<G4endl;
//      }
    } 
  }  
}
////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPhysicsList::SetCutInDepthForAllLayers(G4double aDepth)
{
 G4RunManager *runManager = G4RunManager::GetRunManager();
 geometry = (PLANETOCOSGeometryConstruction*)(runManager->GetUserDetectorConstruction());
 std::vector<G4LogicalVolume*> theAtmosphereLogicalVolumes =
                           geometry->GetAtmosphereLogicalVolumeVector();
 unsigned int nsize = production_cuts.size(); 		   
 
 for (unsigned int i =0; i<theAtmosphereLogicalVolumes.size(); i++){
	
	G4double density = theAtmosphereLogicalVolumes[i]->GetMaterial()->GetDensity();
	G4String name =theAtmosphereLogicalVolumes[i]->GetName();
	G4double cut = aDepth/density;
	
	G4Region* aRegion=0;
	
   	if (nsize  < 1){
	       
   		aRegion = new G4Region(name);
		theAtmosphereLogicalVolumes[i]->SetRegion(aRegion);
        	aRegion->AddRootLogicalVolume(theAtmosphereLogicalVolumes[i]);
		
		production_cuts.push_back(new G4ProductionCuts());
		
		production_cuts[i]->SetProductionCut(cut);
		
		aRegion->SetProductionCuts(production_cuts[i]);
		
   	}
   	else  {    
	    	production_cuts[i]->SetProductionCut(cut);
		if (aRegion->GetNumberOfRootVolumes() <1) {
			aRegion->AddRootLogicalVolume(theAtmosphereLogicalVolumes[i]);
		}
	}    
        	
 }
    
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPhysicsList::SetCutInDepthForASpecificLayer(G4double aDepth , G4String volumeName)
{  G4LogicalVolumeStore* theVolumeStore = G4LogicalVolumeStore::GetInstance();
   G4int index=-1;
   for (unsigned int i=0; i<theVolumeStore->size();i++){
   	G4cout<<(*theVolumeStore)[i]->GetName()<<std::endl;
	if (volumeName == (*theVolumeStore)[i]->GetName()){
		index= int(i);
		i=theVolumeStore->size();
		
	}
   
   }
   G4cout<<index<<std::endl;
   if (index != -1) {
   	G4double density = (*theVolumeStore)[index]->GetMaterial()->GetDensity();
	G4double cut =aDepth/density;
	G4cout<<cut/cm<<std::endl;
	G4Region* aRegion = (*theVolumeStore)[index]->GetRegion();
	if (aRegion->GetName() != volumeName) aRegion = new G4Region(volumeName);
	(*theVolumeStore)[index]->SetRegion(aRegion);
	aRegion->SetProductionCuts(new G4ProductionCuts());
	aRegion->AddRootLogicalVolume((*theVolumeStore)[index]);		
	aRegion->GetProductionCuts()->SetProductionCut(cut);
   
   
   }
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPhysicsList::SetCutInLengthForASpecificLayer(G4double aCut , G4String volumeName)
{  G4LogicalVolumeStore* theVolumeStore = G4LogicalVolumeStore::GetInstance();
   G4int index=-1;
   for (unsigned int i=0; i<theVolumeStore->size();i++){
   	G4cout<<(*theVolumeStore)[i]->GetName()<<std::endl;
   	if (volumeName == (*theVolumeStore)[i]->GetName()){
		index= int(i);
		i=theVolumeStore->size();
		
	}
   
   }
   G4cout<<index<<std::endl;
   if (index != -1) {
	G4Region* aRegion = (*theVolumeStore)[index]->GetRegion();
	if (aRegion->GetName() != volumeName) aRegion = new G4Region(volumeName);
	(*theVolumeStore)[index]->SetRegion(aRegion);
	aRegion->SetProductionCuts(new G4ProductionCuts());
	aRegion->AddRootLogicalVolume((*theVolumeStore)[index]);		
	aRegion->GetProductionCuts()->SetProductionCut(aCut);
	
   
   
   }
}
////////////////////////////////////////////////////////////////////////////////
void PLANETOCOSPhysicsList::SetCutInDepthForAllLayersAndForParticle(G4double , G4String )
{;
}
////////////////////////////////////////////////////////////////////////////////
void PLANETOCOSPhysicsList::DeleteAllRegions()
{
 
 G4RunManager *runManager = G4RunManager::GetRunManager();
 geometry = (PLANETOCOSGeometryConstruction*)(runManager->GetUserDetectorConstruction());
 std::vector<G4LogicalVolume*> theAtmosphereLogicalVolumes =
                           geometry->GetAtmosphereLogicalVolumeVector();
 
 
 for (unsigned int i =0; i<theAtmosphereLogicalVolumes.size(); i++){
	G4Region* aRegion;
	aRegion = G4RegionStore::GetInstance()
	            ->GetRegion(theAtmosphereLogicalVolumes[i]->GetName());
	if (aRegion) 
	       aRegion->RemoveRootLogicalVolume(theAtmosphereLogicalVolumes[i]); 	    
 }
 
 //G4RegionStore::GetInstance()->Clean();
 //production_cuts.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPhysicsList::SetGlobalCuts()
{
  if (verboseLevel>0) {
    G4cout <<"PLANETOCOSPhysicsList::SetGlobalCuts:" <<G4endl;
    G4cout <<"  GlobalCutLength :   "
	   <<G4BestUnit(defaultCutValue,"Length") <<G4endl;
    G4cout <<"  GammaCutLength :    "
	   <<G4BestUnit(cutForGamma,"Length") <<G4endl;
    G4cout <<"  ElectronCutLength : "
	   <<G4BestUnit(cutForElectron,"Length") <<G4endl;
    G4cout <<"  PositronCutLength :   "
	   <<G4BestUnit(cutForPositron,"Length") <<G4endl;
  } 
  //
  //
  // Set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma.
  //
  //if (cutForPositron == cutForElectron == cutForGamma == defaultCutValue) {
  //  SetCutsWithDefault();
  //} else {
    SetCutValue(cutForGamma,"gamma");
    SetCutValue(cutForElectron,"e-");
    SetCutValue(cutForPositron,"e+");
  //}
  //  if (verboseLevel>0) DumpCutValuesTable();
}
////////////////////////////////////////////////////////////////////////////////
//


void PLANETOCOSPhysicsList::SetCuts()
{

  // set global values first
  SetGlobalCuts() ;
  // deal with regions which have specified cuts
 // SetCutsByRegion();
  //
  if (verboseLevel>0)   DumpCutValuesTable();
}

