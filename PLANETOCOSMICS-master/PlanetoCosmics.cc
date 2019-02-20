#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UIGAG.hh"
#include "G4UItcsh.hh"
#include "G4UIXm.hh"
#include "G4UIXaw.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"


#include "PlanetManager.hh"
#include "PlanetUnits.hh"


#include "PLANETOCOSGeometryConstruction.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSStackingAction.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PLANETOCOSApplicationScenario.hh"
#include "PLANETOCOSEventAction.hh"
#include "PLANETOCOSTrackingAction.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSVisManager.hh"
#include "PLANETOCOSPhysicsList.hh"
#include "BlineTool.hh"
#include "MagneticShieldingTool.hh"
#include "DurationManager.hh"




int
main(int argc,char** argv) {




  //Planet Manager
  PlanetManager * thePlanetManager = PlanetManager::GetInstance();
#ifndef DEBUG_MODE
  if(argc==1) {
  	G4cout<<"You should give the name of the planet that you want to ";
	G4cout<<" consider in the simulation as first argument of the ";
	G4cout<<"PLANETOCOSMICS code"<<std::endl;  
	
  }
  else if (G4String(argv[1]) == "Earth"|| 
           G4String(argv[1]) == "Mars" || 
	   G4String(argv[1]) == "Mercury" 
#ifdef USE_JUPITER
	  || G4String(argv[1]) == "Jupiter" 	
#endif	   
	){
 
  	G4String planet_name = G4String(argv[1]);
#else
	
	G4String planet_name = G4String("Earth");
	argc=2;
#endif	
//Selection of the planet
//--------------------------	
	//the unit table is also define when selecting the planet    
	
	thePlanetManager->SelectPlanet(planet_name);
   
 //Duration manager
  DurationManager* theDurationManager;
  theDurationManager = DurationManager::GetInstance();
 



// Run manager
  	G4RunManager * runManager = new G4RunManager;
  

// UserInitialization classes
//--------------------------

  	PLANETOCOSGeometryConstruction* theGeometry = new PLANETOCOSGeometryConstruction;  
  	runManager->SetUserInitialization(theGeometry);  
	runManager->SetUserInitialization(new PLANETOCOSPhysicsList);
 
#ifdef G4VIS_USE

// Visualization, if you choose to have it!

  	G4VisManager* visManager = new PLANETOCOSVisManager;
 	visManager->Initialize();
	
#endif
  
  // UserAction classes
  //-----------------------
       
  	PLANETOCOSPrimaryGeneratorAction* thePrimaryAction  = new PLANETOCOSPrimaryGeneratorAction();
  	
	runManager->SetUserAction(thePrimaryAction);
  	PLANETOCOSAnalysisManager::GetInstance()->SetPrimaryAction(thePrimaryAction);
         
  	PLANETOCOSSteppingAction* theSteppingAction = new PLANETOCOSSteppingAction();
  	PLANETOCOSStackingAction* theStackingAction = new PLANETOCOSStackingAction();
  
  	theStackingAction->SetUntrackedParticleDic(theSteppingAction->GetUntrackedParticleDic());
  
  	PLANETOCOSEventAction* theEventAction = new PLANETOCOSEventAction();
	theEventAction->SetSteppingAction(theSteppingAction);
  	PLANETOCOSApplicationScenario* theRunAction = new PLANETOCOSApplicationScenario();
  
  	runManager->SetUserAction(theSteppingAction);
  	runManager->SetUserAction(theStackingAction);
  	runManager->SetUserAction(theEventAction);
  	runManager->SetUserAction(theRunAction);  
	runManager->SetUserAction(new PLANETOCOSTrackingAction());   
   			      
//Initialize G4 kernel
//-----------------------
//runManager->Initialize(); 
  
  
//Analysis manager
//---------------------
  
  	PLANETOCOSAnalysisManager* theAnalysisManager = 
        				PLANETOCOSAnalysisManager::GetInstance();
	theAnalysisManager->SetEventAction(theEventAction);
	theAnalysisManager->SetSteppingAction(theSteppingAction);	  
// User interface
//--------------------
  	G4UImanager * UI = G4UImanager::GetUIpointer();
	
//Magnetic shielding tool
//-----------------------
	MagneticShieldingTool* theMagneticShieldingTool;
	theMagneticShieldingTool = MagneticShieldingTool::GetInstance();	
  
  
  
  
//Bline tool
//------------------------ 
  	/*BlineTool* theBlineTool;
  	theBlineTool = new BlineTool(thePrimaryAction,
  			       theStackingAction,
			       theSteppingAction,
			       0,
			       theEventAction,
			       theRunAction);
  	theBlineTool->AddStopVolume(planet_name);*/			       
  
  	if(argc==2){ //Interactive mode	
		G4UIsession * session = new G4UIterminal(new G4UItcsh);
    		session->SessionStart(); 
    		delete session;
  	}
  	else { // Batch mode
    		G4String command = "/control/execute ";
    		G4String fileName = argv[2];
    		UI->ApplyCommand(command+fileName);
  	}

#ifdef G4VIS_USE
  	delete visManager;
#endif
 	delete runManager; 
  	delete theAnalysisManager;
  	delete thePlanetManager;
#ifndef DEBUG_MODE 
 }
  else {
  	G4cout<<argv[1]<<" is not in the list of planet available"<<std::endl;
  
  
  }
#endif   


 

  return 0;
}


