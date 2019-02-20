#include"PLANETOCOSScenarioMessenger.hh"
#include"PLANETOCOSApplicationScenario.hh"
#include"PLANETOCOSPrimaryGeneratorAction.hh"
#include"G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//Planet
//------

#include "PlanetMagneticField.hh"
#include "PlanetManager.hh"

//units
#include "G4UnitsTable.hh"


////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSScenarioMessenger::PLANETOCOSScenarioMessenger (PLANETOCOSApplicationScenario* ApplicationScenario )
{ 
  theApplicationScenario = ApplicationScenario;
  
  
  G4String candidates, cmd_name, guidance;
  
  
  
  //General run commands
  
  cmd_name = "/PLANETOCOS/RANDOM/FixTheSeedRandomly";
  SetRandomSeedCmd= new G4UIcmdWithoutParameter(cmd_name,this);
  SetRandomSeedCmd->SetGuidance("Compute a new seed randomly depending on time");
  SetRandomSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/RANDOM/SetRandomSeedAtRunStart";
  SetRandomSeedAtRunStartCmd= new G4UIcmdWithABool(cmd_name,this);
  SetRandomSeedAtRunStartCmd->
      SetGuidance("If true the seed is recomputed randomly at beginning of a run ");
  SetRandomSeedAtRunStartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}  
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSScenarioMessenger::~PLANETOCOSScenarioMessenger()
{
}		  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSScenarioMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if (command == SetRandomSeedCmd) theApplicationScenario->SetRandomSeed();	  
   
  else if (command == SetRandomSeedAtRunStartCmd) 
    	theApplicationScenario->SetRandomSeedNeeded(SetRandomSeedAtRunStartCmd->GetNewBoolValue(newValues));  
  
      
}	  

	    				      












