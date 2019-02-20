#include "PlanetMessenger.hh"
#include "PlanetManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserStackingAction.hh"
#include "PLANETOCOSStackingAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinatePlanet.hh"


PlanetMessenger::PlanetMessenger(PlanetManager* theManager)
                                :pPlanetManager(theManager)
{
  
  
  //command directory 
  //-------------------
  
  planetDir = new G4UIdirectory("/PLANETOCOS/");
  planetDir->SetGuidance("Planetocosmics code control");
  
  testDir = new G4UIdirectory("/PLANETOCOS/SPACECOORDINATE/");
  testDir->SetGuidance("Control the space coordinate convertor");
  
  //parameters
  //---------
   G4UIparameter* input_file_param = new G4UIparameter("InputFile name",'s',false);
 
   G4UIparameter* output_file_param = new G4UIparameter("OutputFile name",'s',false);
 
  
  //commands
  //--------
 /* G4String title;
  SelectPlanetCmd =  new G4UIcmdWithAString("/PLANETOCOS/SelectPlanet",this);
  SelectPlanetCmd->SetGuidance("Select a planet for the simulation");
  SelectPlanetCmd->SetGuidance("Earth Mars Mercury"); 
  SelectPlanetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
  
  //Space ccordinate test
  //--------------------

  G4String cmd_name;

//#ifdef USE_SPICE  
  cmd_name = "/PLANETOCOS/SPACECOORDINATE/UseSpice";
  SetUseSpiceCmd =  new G4UIcmdWithABool(cmd_name,this);
  SetUseSpiceCmd->SetGuidance("If true (false) the Spice library will (not) be used for computing coordinate transformation");  
  SetUseSpiceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
//#endif  
   
  cmd_name = "/PLANETOCOS/SPACECOORDINATE/MarinerOrbitTest";
  TestMarinerOrbitCmd = new G4UIcommand(cmd_name,this);
  TestMarinerOrbitCmd->SetGuidance("Mariner Orbit Test"); 
  TestMarinerOrbitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  TestMarinerOrbitCmd->SetParameter(input_file_param);
  TestMarinerOrbitCmd->SetParameter(output_file_param);
  

  cmd_name = "/PLANETOCOS/SPACECOORDINATE/MarsOdysseyOrbitTest";
  TestMarsOdysseyCmd = new G4UIcommand(cmd_name,this);
  TestMarsOdysseyCmd->SetGuidance("Mariner Orbit Test"); 
  TestMarsOdysseyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  TestMarsOdysseyCmd->SetParameter(input_file_param);
  TestMarsOdysseyCmd->SetParameter(output_file_param);
  
  
  cmd_name = "/PLANETOCOS/SPACECOORDINATE/ImageOrbitTest";
  TestImageOrbitCmd = new G4UIcommand(cmd_name,this);
  TestImageOrbitCmd->SetGuidance("Image Orbit Test"); 
  TestImageOrbitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  TestImageOrbitCmd->SetParameter(input_file_param);
  TestImageOrbitCmd->SetParameter(output_file_param);
  
  
  cmd_name = "/PLANETOCOS/SPACECOORDINATE/WindOrbitTest";
  TestWindOrbitCmd = new G4UIcommand(cmd_name,this);
  TestWindOrbitCmd->SetGuidance("Wind Orbit Test"); 
  TestWindOrbitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  TestWindOrbitCmd->SetParameter(input_file_param);
  TestWindOrbitCmd->SetParameter(output_file_param);
  

}
////////////////////////////////////////////////////////////////////////////////
//
PlanetMessenger::~PlanetMessenger()
{delete planetDir;
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
 /*if (command == SelectPlanetCmd)  
     			pPlanetManager->SelectPlanet(newValues);*/
  
 if ( command == TestMarinerOrbitCmd ){
  	const char* paramString=newValues;
        G4String  InputFile, OutputFile;
        std::istringstream is((char*)paramString);
        is >> InputFile >> OutputFile;
	SpaceCoordinatePlanet::GetInstance()->TestMarinerOrbit(InputFile,OutputFile );      
 }
 else if ( command == TestMarsOdysseyCmd ){
  	const char* paramString=newValues;
        G4String  InputFile, OutputFile;
        std::istringstream is((char*)paramString);
        is >> InputFile >> OutputFile;
	SpaceCoordinatePlanet::GetInstance()->TestMarsOdysseyOrbit(InputFile,OutputFile );      
 }
 else if ( command == TestImageOrbitCmd ){
  	const char* paramString=newValues;
        G4String  InputFile, OutputFile;
        std::istringstream is((char*)paramString);
        is >> InputFile >> OutputFile;
	SpaceCoordinatePlanet::GetInstance()->TestImageOrbit(InputFile,OutputFile );      
 }
 else if ( command == TestWindOrbitCmd ){
  	const char* paramString=newValues;
        G4String  InputFile, OutputFile;
        std::istringstream is((char*)paramString);
        is >> InputFile >> OutputFile;
	SpaceCoordinatePlanet::GetInstance()->TestWindOrbit(InputFile,OutputFile );      
 }

//#ifdef USE_SPICE 
 else if ( command == SetUseSpiceCmd ){
  	G4bool aBool = SetUseSpiceCmd->GetNewBoolValue(newValues);
	SpaceCoordinatePlanet::GetInstance()->SetUseSpice(aBool);      
 }
//#endif 			      			      

}

