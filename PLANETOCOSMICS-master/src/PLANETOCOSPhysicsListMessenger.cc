////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSPhysicsListMessenger.hh"

#include "G4ios.hh"
#include "G4Tokenizer.hh"
#include <iomanip>               
#include <sstream>
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "PLANETOCOSAnalysisManager.hh"
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPhysicsListMessenger::PLANETOCOSPhysicsListMessenger(PLANETOCOSPhysicsList* pPhys)
  :pPhysicsList(pPhys)
{

  PhysDir = new G4UIdirectory("/PLANETOCOS/PHYSICS/");
  PhysDir->SetGuidance(" Control of the physics processes");
  
  CutInRangeDir = new G4UIdirectory("/PLANETOCOS/CUT/");
  CutInRangeDir->SetGuidance("Define cut in range by regions");
 
  //parameters
  
  G4UIparameter* depth_param = new G4UIparameter("Depth",'d',false);
  depth_param->SetParameterRange("Depth >0");
  
  G4UIparameter* depth_unit_param = new G4UIparameter("DepthUnit",'s',false);
  depth_unit_param->SetParameterCandidates("g/cm2 kg/cm2 g/m2 kg/m2");
  
  G4UIparameter* length_param = new G4UIparameter("Length",'d',false);
  length_param->SetParameterRange("Length >0");
  
  G4UIparameter* length_unit_param = new G4UIparameter("LengthUnit",'s',false);
  length_unit_param->SetParameterCandidates("cm m km mm mum");
  
  G4UIparameter* name_param = new G4UIparameter("Name",'s',false);
  
  G4UIparameter* MaxA = new G4UIparameter("MaxA",'i',false);
  MaxA->SetParameterRange("MaxA >0");
  
  G4UIparameter* MaxZ = new G4UIparameter("MaxZ",'i',false);
  MaxZ->SetParameterRange("MaxZ >0");
 
  G4UIparameter* facrange_param = new G4UIparameter("Facrange",'d',false);
  facrange_param->SetParameterRange("Facrange > 0");
  
  G4UIparameter* bool_param = new G4UIparameter("aBool",'b',false);
  
  
 
  G4String cmd_title = "/PLANETOCOS/PHYSICS/SelectTypeOfEMPhysics";
  SelectTypeOfEMPhysicsCmd = new G4UIcmdWithAString(cmd_title,this);
  SelectTypeOfEMPhysicsCmd->SetGuidance(" select the type of electromagnetic physics model");
  SelectTypeOfEMPhysicsCmd->SetGuidance("            standard: EM interactions using standard G4EM ");
  SelectTypeOfEMPhysicsCmd->SetGuidance("            low: EM interactions using low energy G4EM ");
  SelectTypeOfEMPhysicsCmd->SetGuidance("            none: no EM  and no hadronic interactions ");
  
  G4String candidates =" standard STANDARD lowenergy LOWENERGY none NONE ";
  SelectTypeOfEMPhysicsCmd->SetCandidates(candidates);
  SelectTypeOfEMPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title = "/PLANETOCOS/PHYSICS/SelectTypeOfHadronicPhysics";
  SelectTypeOfHadronicPhysicsCmd = new G4UIcmdWithAString(cmd_title,this);
  SelectTypeOfHadronicPhysicsCmd->SetGuidance(" select the type of hadronic physics model");
  
  candidates =" NOHADRONIC LHEP  LHEP_HP LHEP_BIC LHEP_BIC_HP LHEP_BERT LHEP_BERT_HP ";
  candidates +="QGSP  QGSP_HP QGSP_BIC QGSP_BIC_HP QGSP_BERT QGSP_BERT_HP ";
  candidates +="QGSC  QGSC_HP QGSC_BIC QGSC_BIC_HP QGSC_BERT QGSC_BERT_HP ";
  candidates +="FTFP  FTFP_HP FTFP_BIC FTFP_BIC_HP FTFP_BERT FTFP_BERT_HP ";
  candidates +="FTFC  FTFC_HP FTFC_BIC FTFC_BIC_HP FTFC_BERT FTFC_BERT_HP ";
  
   
  SelectTypeOfHadronicPhysicsCmd->SetCandidates(candidates);
  SelectTypeOfHadronicPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title = "/PLANETOCOS/PHYSICS/SelectTypeOfIonHadronicPhysics";
  SelectTypeOfIonHadronicPhysicsCmd = new G4UIcmdWithAString(cmd_title,this);
  SelectTypeOfIonHadronicPhysicsCmd->SetGuidance(" select the type of hadronic physics model for ions");
  candidates =" NONE LE BIC ABRASION"; 
  SelectTypeOfIonHadronicPhysicsCmd->SetCandidates(candidates);
  SelectTypeOfIonHadronicPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title = "/PLANETOCOS/PHYSICS/ConsiderElectromagneticNuclearPhysics";
  SetConsiderElectroNuclearPhysicsCmd = new G4UIcmdWithABool(cmd_title,this);
  SetConsiderElectroNuclearPhysicsCmd->SetGuidance("If true (false) the electromagnetic nuclear processes are (are not) taken into account.");
  SetConsiderElectroNuclearPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title = "/PLANETOCOS/PHYSICS/ConsiderMuonNuclearPhysics";
  SetConsiderMuonNuclearPhysicsCmd = new G4UIcmdWithABool(cmd_title,this);
  SetConsiderMuonNuclearPhysicsCmd->SetGuidance("If true (false) the Muon nuclear interactions are (are not) taken into account.");
  SetConsiderMuonNuclearPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title = "/PLANETOCOS/PHYSICS/ConsiderSynchrotronRadiation";
  SetConsiderSyncPhysicsCmd = new G4UIcmdWithABool(cmd_title,this);
  SetConsiderSyncPhysicsCmd->SetGuidance("If true (false) the  synchrotron radaition is (are not) taken into account.");
  SetConsiderSyncPhysicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  cmd_title = "/PLANETOCOS/Initialise";
  BuildListCmd =  new  G4UIcmdWithoutParameter(cmd_title, this);
  BuildListCmd->SetGuidance("Construct the geometry and build and construct the physics list and ini");
  BuildListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 

  //cut in range commands
  //-------------------
  cmd_title = "/PLANETOCOS/CUT/SetCutInDepthForAllAtmosphericLayers";
  SetCutInDepthForAllLayersCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
  SetCutInDepthForAllLayersCmd->SetGuidance("Define the cut in depth for all layers");
  SetCutInDepthForAllLayersCmd->SetUnitCategory("Depth");
  SetCutInDepthForAllLayersCmd->SetParameterName("depth", false);
  SetCutInDepthForAllLayersCmd->SetRange(" depth > 0. ");
  SetCutInDepthForAllLayersCmd->AvailableForStates(G4State_Idle);

  cmd_title = "/PLANETOCOS/CUT/RemoveRegions";
  DeleteAllRegionsCmd =  new  G4UIcmdWithoutParameter(cmd_title, this);
  DeleteAllRegionsCmd->SetGuidance("Deselect the cut per region mode");
  DeleteAllRegionsCmd->AvailableForStates(G4State_Idle);
  
  cmd_title= "/PLANETOCOS/CUT/SetCutInDepthForAGivenVolume";
  SetCutInDepthForASpecificLayerCmd = new G4UIcommand(cmd_title,this);
  SetCutInDepthForASpecificLayerCmd->SetGuidance("Define the cut in range of a given volume by a depth");
  SetCutInDepthForASpecificLayerCmd->SetParameter(name_param);
  SetCutInDepthForASpecificLayerCmd->SetParameter(depth_param);
  SetCutInDepthForASpecificLayerCmd->SetParameter(depth_unit_param);
  SetCutInDepthForASpecificLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title= "/PLANETOCOS/CUT/SetCutInLengthForAGivenVolume";
  SetCutInLengthForASpecificLayerCmd = new G4UIcommand(cmd_title,this);
  SetCutInLengthForASpecificLayerCmd->SetGuidance("Define the cut in range of a given volume by a length");
  SetCutInLengthForASpecificLayerCmd->SetParameter(name_param);
  SetCutInLengthForASpecificLayerCmd->SetParameter(length_param);
  SetCutInLengthForASpecificLayerCmd->SetParameter(length_unit_param);
  SetCutInLengthForASpecificLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title= "/PLANETOCOS/PHYSICS/SetMinEForMultriFrag";
  SetMinEForMultiFragCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
  SetMinEForMultiFragCmd
    ->SetGuidance("Define the minimum energy per nucleon for the multi fragmentation model");
  SetMinEForMultiFragCmd->SetUnitCategory("Energy");
  SetMinEForMultiFragCmd->SetParameterName("Emin", false);
  SetMinEForMultiFragCmd->SetRange(" Emin > 0. ");
  SetMinEForMultiFragCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_title= "/PLANETOCOS/PHYSICS/SetMaxAandZForFermiBreakUp";
  SetMaxAandZForFermiBreakUpCmd = new G4UIcommand(cmd_title,this);
  SetMaxAandZForFermiBreakUpCmd->SetParameter(MaxA);
  SetMaxAandZForFermiBreakUpCmd->SetParameter(MaxZ);
  SetMaxAandZForFermiBreakUpCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  cmd_title= "/PLANETOCOS/PHYSICS/SetSteppingMscAlgorithmParameters";
  SetSteppingAlgorithmParametersForMscCmd = new G4UIcommand(cmd_title,this);
  SetSteppingAlgorithmParametersForMscCmd->SetParameter(bool_param);
  SetSteppingAlgorithmParametersForMscCmd->SetParameter(facrange_param);
  SetSteppingAlgorithmParametersForMscCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPhysicsListMessenger::~PLANETOCOSPhysicsListMessenger()
{
  delete PhysDir;
  delete SelectTypeOfEMPhysicsCmd;
  delete SelectTypeOfHadronicPhysicsCmd;
  delete ShowTypeOfPhysicsCmd;
  delete SetCutInDepthForAllLayersCmd;
  delete DeleteAllRegionsCmd;
  delete SetMinEForMultiFragCmd;
  delete SetMaxAandZForFermiBreakUpCmd;
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == SelectTypeOfEMPhysicsCmd) {
  	newValue.toUpper();
  	pPhysicsList->SelectEMPhysicsType(newValue) ; 
  }
  else if (command == SelectTypeOfHadronicPhysicsCmd){
  	newValue.toUpper();
    	pPhysicsList->SelectHadronicPhysicsType(newValue); 
  }
  else if (command == SelectTypeOfIonHadronicPhysicsCmd){
  	newValue.toUpper();
    	pPhysicsList->SelectLightIonHadronicPhysicsType(newValue); 
  }
  else if (command == SetConsiderElectroNuclearPhysicsCmd){
  	G4bool aBool =SetConsiderElectroNuclearPhysicsCmd->GetNewBoolValue(newValue);
	pPhysicsList->SetConsiderElectromagneticNuclearPhysics(aBool); 
  }
  else if (command == SetConsiderMuonNuclearPhysicsCmd){
  	G4bool aBool =SetConsiderMuonNuclearPhysicsCmd->GetNewBoolValue(newValue);
	pPhysicsList->SetConsiderMuonNuclearPhysics(aBool); 
  }
  else if (command == SetConsiderSyncPhysicsCmd){
  	G4bool aBool =SetConsiderSyncPhysicsCmd->GetNewBoolValue(newValue);
	pPhysicsList->SetConsiderSynchrotronRadiation(aBool); 
  }
  
  
  
  else if (command == BuildListCmd){
  	pPhysicsList->BuildList();
	G4RunManager::GetRunManager()->Initialize();
	
	
  } 
  else if (command == SetCutInDepthForAllLayersCmd){
        G4double val = SetCutInDepthForAllLayersCmd->GetNewDoubleValue(newValue);
   	pPhysicsList->SetCutInDepthForAllLayers(val);
  }
  else if (command == DeleteAllRegionsCmd){
  	pPhysicsList->DeleteAllRegions();
  }
  else if (command == SetCutInDepthForASpecificLayerCmd){
  	const char* paramString=newValue;
     	G4double depth;
	G4String depth_unit,volume_name;
     	std::istringstream is((char*)paramString);
     	is>>volume_name>>depth>>depth_unit; 
	G4cout<<depth_unit<<std::endl;
	depth*=G4UnitDefinition::GetValueOf(depth_unit);
	G4cout<<depth<<std::endl; 
	pPhysicsList->SetCutInDepthForASpecificLayer(depth,volume_name);
  }
  else if (command == SetCutInLengthForASpecificLayerCmd){
  	const char* paramString=newValue;
     	G4double length;
	G4String length_unit,volume_name;
     	std::istringstream is((char*)paramString);
     	is>>volume_name>>length>>length_unit;
	G4cout<<length_unit<<std::endl; 
	length*=G4UnitDefinition::GetValueOf(length_unit);
	G4cout<<length<<std::endl; 
	pPhysicsList->SetCutInLengthForASpecificLayer(length,volume_name);
  }
  
  else if (command == SetMinEForMultiFragCmd){
  	G4double val = SetMinEForMultiFragCmd->GetNewDoubleValue(newValue);
  	pPhysicsList->SetMinEForMultiFrag(val);
  }
  else if (command == SetMaxAandZForFermiBreakUpCmd){
  	const char* paramString=newValue;
     	G4int MaxA, MaxZ;
     	std::istringstream is((char*)paramString);
     	is>>MaxA>>MaxZ;
	pPhysicsList->SetMaxAandZForFermiBreakUp(MaxA,MaxZ);
  }
  else if (command == SetSteppingAlgorithmParametersForMscCmd){
  	const char* paramString=newValue;
     	G4bool aBool;
	G4double facrange;
     	std::istringstream is((char*)paramString);
     	is>>aBool>>facrange;
	G4cout<<"Stepping algorithm "<<aBool<<'\t'<<facrange<<std::endl;
	
	G4cout<<"test "<<std::endl; 
	pPhysicsList->SetSteppingAlgorithmParameters(aBool,facrange);
  }
  
 	 
   
}
////////////////////////////////////////////////////////////////////////////////
