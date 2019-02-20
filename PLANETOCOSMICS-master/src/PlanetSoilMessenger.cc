#include"PlanetSoil.hh"
#include"PlanetSoilMessenger.hh"
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

//units
#include "G4UnitsTable.hh"

PlanetSoilMessenger::PlanetSoilMessenger (PlanetSoil* aSoil )
{ 
  theSoil = aSoil;

  G4String candidates;
  G4String guidance;
  G4String cmd_name;
 
//directories
///////////////
 
  SoilDir= new G4UIdirectory("/PLANETOCOS/SOIL/");
  SoilDir->SetGuidance("Definition of the soil");
  
 
//parameters
/////////////
  
  
  G4UIparameter* el_name_or_symbol_param = new G4UIparameter("Name",'s',false);
 
  
  
  G4UIparameter* density_unit_param = new G4UIparameter("Density Unit",'s',false);
  density_unit_param->SetParameterCandidates("g/cm3 mg/cm3 kg/m3");

  G4UIparameter* depth_unit_param = new G4UIparameter("DepthOrLength Unit",'s',false);
  depth_unit_param->SetParameterCandidates("g/cm2 g/m2 kg/m2 kg/cm2 cm dm km m");

  G4UIparameter* density_param = new G4UIparameter("Density",'d',false);
  density_param->SetParameterRange("Density > 0");
  
  G4UIparameter* depth_param = new G4UIparameter("DepthOrThickness",'d',false);
  depth_param->SetParameterRange("DepthOrThickness > 0");
  
  G4UIparameter* concentration_param = new G4UIparameter("Concentration",'d',false);
  concentration_param->SetParameterRange("Concentration > 0 && Concentration<=1.");
  
  G4UIparameter* nel_param = new G4UIparameter("nb_elements",'i',false);
  nel_param->SetParameterRange("nb_elements >0");
  
  G4UIparameter* nb_sub_layers_param = new G4UIparameter("nb_sub_layers",'i',true);
  nb_sub_layers_param->SetParameterRange("nb_sub_layers >0");
  nb_sub_layers_param->SetDefaultValue(1);
  
    
// commands
/////////////////////////////////

 
  cmd_name = "/PLANETOCOS/SOIL/AddMonoElementLayer";
  AddMonoElementLayerCmd = new G4UIcommand(cmd_name,this);
  AddMonoElementLayerCmd->SetGuidance("Add a soil layer made of one element"); 
  AddMonoElementLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  AddMonoElementLayerCmd->SetParameter(el_name_or_symbol_param);
  AddMonoElementLayerCmd->SetParameter(density_param);
  AddMonoElementLayerCmd->SetParameter(density_unit_param); 
  AddMonoElementLayerCmd->SetParameter(depth_param);
  AddMonoElementLayerCmd->SetParameter(depth_unit_param);
  AddMonoElementLayerCmd->SetParameter(nb_sub_layers_param);
  
  cmd_name = "/PLANETOCOS/SOIL/AddLayer";
  AddLayerCmd = new G4UIcommand(cmd_name,this);
  AddLayerCmd->SetGuidance("Add a soil layer made of several elements"); 
  AddLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  AddLayerCmd->SetParameter(nel_param);
  AddLayerCmd->SetParameter(density_param);
  AddLayerCmd->SetParameter(density_unit_param); 
  AddLayerCmd->SetParameter(depth_param);
  AddLayerCmd->SetParameter(depth_unit_param);
  AddLayerCmd->SetParameter(nb_sub_layers_param);
  
  cmd_name = "/PLANETOCOS/SOIL/AddElementToLayer";
  AddElementToLayerCmd = new G4UIcommand(cmd_name,this);
  AddElementToLayerCmd->SetGuidance("Add an element to the soil layer being defined");
  AddElementToLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  AddElementToLayerCmd->SetParameter(el_name_or_symbol_param);
  AddElementToLayerCmd->SetParameter(concentration_param);
  
  cmd_name = "/PLANETOCOS/SOIL/ResetLayers";
  ResetLayersCmd  = new G4UIcmdWithoutParameter(cmd_name,this);
  ResetLayersCmd->SetGuidance("Remove all soil layers");
  ResetLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
}
////////////////////////////////////////////////////////////////////////////////
//
PlanetSoilMessenger::~PlanetSoilMessenger()
{ delete AddMonoElementLayerCmd;
  
}		  
////////////////////////////////////////////////////////////////////////////////
//
void PlanetSoilMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
   
  if (command == AddMonoElementLayerCmd){
  	const char* paramString=newValues;
	G4String el_name_or_symbol, density_unit, depth_or_thickness_unit;
	G4double density, depth_or_thickness;
	G4int nb_sub_layers;
	std::istringstream is((char*)paramString);
        is >> el_name_or_symbol >> density >> density_unit 
	   >> depth_or_thickness >> depth_or_thickness_unit>>nb_sub_layers;
	density *= G4UnitDefinition::GetValueOf(density_unit);
	depth_or_thickness *= G4UnitDefinition::GetValueOf(depth_or_thickness_unit);   
  	if (G4UnitDefinition::GetCategory(depth_or_thickness_unit) == "Depth") {
		theSoil->AddMonoElementLayerAndSetDepth(el_name_or_symbol, density,depth_or_thickness,nb_sub_layers);
	}
	else if (G4UnitDefinition::GetCategory(depth_or_thickness_unit) == "Length"){
		theSoil->AddMonoElementLayerAndSetThickness(el_name_or_symbol, density,depth_or_thickness,nb_sub_layers);
	} 
	else G4cout<< depth_or_thickness_unit+" is not a Length or Depth unit."<<std::endl; 
  } 
  else  if (command == AddLayerCmd){
  	const char* paramString=newValues;
	G4int nb_elements, nb_sub_layers;
	G4String  density_unit, depth_or_thickness_unit;
	G4double density, depth_or_thickness;
	std::istringstream is((char*)paramString);
	is >> nb_elements>> density >> density_unit 
	   >> depth_or_thickness >> depth_or_thickness_unit>>nb_sub_layers;
	density *= G4UnitDefinition::GetValueOf(density_unit);
	depth_or_thickness *= G4UnitDefinition::GetValueOf(depth_or_thickness_unit);   
  	
	if (G4UnitDefinition::GetCategory(depth_or_thickness_unit) == "Depth") {
		theSoil->AddLayerAndSetDepth(nb_elements, density,depth_or_thickness,nb_sub_layers);
	}
	else if (G4UnitDefinition::GetCategory(depth_or_thickness_unit) == "Length"){
		theSoil->AddLayerAndSetThickness(nb_elements, density,depth_or_thickness,nb_sub_layers);
	} 
	else G4cout<< depth_or_thickness_unit+" is not a Length or Depth unit."<<std::endl; 
  } 
  else if (command == AddElementToLayerCmd){
  	const char* paramString=newValues;
	G4String el_name_or_symbol;
	G4double concentration;
	std::istringstream is((char*)paramString);
        is >> el_name_or_symbol >> concentration;
	theSoil->AddElementToLayer(el_name_or_symbol, concentration);
  	
  }	 
  else if (command == ResetLayersCmd){
  	theSoil->ResetLayers();
  }
 
  
  //bool result=true; 
   
 /* if (command == SetEpsilonCmd) 
            theMotherField->SetEpsilon(SetEpsilonCmd->GetNewDoubleValue(newValues)); 

  else if (command == SetG4DeltaChordCmd) 
            theMotherField->SetDeltaChord(SetG4DeltaChordCmd->GetNewDoubleValue(newValues));   
   
 
  else if (command == SetDeltaIntersectionCmd) 
            theMotherField->SetDeltaIntersection(SetDeltaIntersectionCmd->GetNewDoubleValue(newValues));   

  else if (command == ResetIntegrationParametersCmd) 
                            theMotherField->ResetIntegrationParameters();  
 
 
  else if ( command == SetStepperCmd )
            theMotherField->SetStepper(newValues);
	    			     
  // magnetic field model parameters 
  
  else if( command == SetTimeOfBCmd)
            theMotherField->SetTimeOfB(SetTimeOfBCmd->GetNewDoubleValue(newValues));
   
  else if ( command == SetStartDateCmd ){
  	const char* paramString=newValues;
        G4int  Year,Month,Day,Hour,Minute,Second;
        std::istringstream is((char*)paramString);
        is >> Year >> Month >> Day >> Hour >> Minute >> Second;
        theMotherField->SetStartDate(Year,Month,Day,Hour,Minute,Second);
  }
  
  
  else if( command == Setnmax_GaussCmd ) 
          theMotherField->Setnm_gauss(Setnmax_GaussCmd->GetNewIntValue(newValues));
	  
  else if ( command == SetInternalFieldCmd )
            theMotherField->SetInternalField(newValues);
 
  else if ( command == SetExternalFieldCmd )
            theMotherField->SetExternalField(newValues);	    
  
  else if( command == SetDipoleB0Cmd ) 
            theMotherField->SetDipoleB0(SetDipoleB0Cmd->GetNewDoubleValue(newValues));
  
  else if( command == SetDipoleShiftCmd){
  	G4ThreeVector vec=SetDipoleShiftCmd->GetNew3VectorValue(newValues);
  	theMotherField->SetDipoleShift(vec);
  }
  
  else if ( command == SetDipoleAxisCmd ){
   	const char* paramString=newValues;
        G4double  Theta,Phi;
	G4String AngleUnit;
        std::istringstream is((char*)paramString);
        is >> Theta >> Phi;
	is >>AngleUnit;
	G4double angle_unit = G4UnitDefinition::GetValueOf(AngleUnit);
        Theta *=angle_unit;
	Phi *=angle_unit;
        theMotherField->SetDipoleAxis(Theta,Phi);
  } 
       
  
  
  else if (command ==ComputeBfieldAtDifferentPositions){
  	const char* paramString=newValues;
       	G4double  alt0,dalt,lat0,dlat,long0,dlong;
	G4int   nlat,nlong, nalt ;
	G4String OutputFile,CoordSys,AngleUnit,LengthUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>>alt0>>dalt>>nalt>>LengthUnit>>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	AngleUnit >> OutputFile;
	alt0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dalt*=G4UnitDefinition::GetValueOf(LengthUnit);
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMotherField->ComputeBfieldAtDifferentPosition
	                 (CoordSys,alt0,dalt,nalt,
			  lat0,dlat,nlat,
			  long0,dlong,nlong,
		          OutputFile);
  
  
  
  }
  else if (command ==ComputeBfieldAtDifferentPositions1){
  	const char* paramString=newValues;
       	G4double  alt0,dalt,x0,dx,y0,dy;
	G4int   nalt,nx, ny ;
	G4String OutputFile,LengthUnit;
        std::istringstream is((char*)paramString);
        is >>alt0>>dalt>>nalt>>
	x0>>dx>>nx>>y0>>dy>>ny>>LengthUnit>>OutputFile;
	alt0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dalt*=G4UnitDefinition::GetValueOf(LengthUnit);
	x0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dx*=G4UnitDefinition::GetValueOf(LengthUnit);
	y0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dy*=G4UnitDefinition::GetValueOf(LengthUnit);
        theMotherField->ComputeBfieldAtDifferentPosition
	                 (alt0,dalt,nalt,
			  x0,dx,nx,
			  y0,dy,ny,
		          OutputFile);
  
  
  
  }
  else if (command ==TraceBlineAtDifferentPositions){
  	const char* paramString=newValues;
       	G4double  alt,lat0,dlat,long0,dlong;
	G4int   nlat,nlong;
	G4String CoordSys,AngleUnit,LengthUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>>alt>>LengthUnit>>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	AngleUnit;
	alt*=G4UnitDefinition::GetValueOf(LengthUnit);
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMotherField->TraceBlineFromDifferentPositions
	                 (CoordSys,alt,
			  lat0,dlat,nlat,
			  long0,dlong,nlong);
  
  
  
  }
  else if (command ==TraceBlineAtDifferentPositions1){
  	const char* paramString=newValues;
       	G4double  alt,x0,dx,y0,dy;
	G4int   nx,ny;
	G4String LengthUnit;
        std::istringstream is((char*)paramString);
        is >>alt>>
	x0>>dx>>nx>>y0>>dy>>ny>>LengthUnit;
	alt*=G4UnitDefinition::GetValueOf(LengthUnit);
	x0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dx*=G4UnitDefinition::GetValueOf(LengthUnit);
	y0*=G4UnitDefinition::GetValueOf(LengthUnit);
	dy*=G4UnitDefinition::GetValueOf(LengthUnit);
        theMotherField->TraceBlineFromDifferentPositions
	                 (alt,
			  x0,dx,nx,
			  y0,dy,ny);
  }
   else if (command == SetHomogeneousFieldCmd){
  	theMotherField->SetBfieldConst(SetHomogeneousFieldCmd->
					GetNew3VectorValue(newValues));
  }
  
  else if (command == SetHomogeneousFieldFromSelectedModelCmd){
  	const char* paramString=newValues;
	G4double  altitude,latitude,longitude;
	G4String coor_sys, length_unit, model;
	std::istringstream is((char*)paramString);
	is >> altitude >> length_unit >> latitude >> longitude >>coor_sys
		>>model;
	altitude *= G4UnitDefinition::GetValueOf(length_unit);
	latitude *=degree;
	longitude *=degree;
	theMotherField->ComputeOrientationAndBfieldConst(altitude,
				  latitude, longitude,
	                          coor_sys,
			           model);
	theMotherField->SetInternalField("CONST");
	
  }
  
  else if (command == SetWorldCenterPositionCmd){
  	const char* paramString=newValues;
	G4double latitude,longitude;
	G4String coor_sys;
	std::istringstream is((char*)paramString);
	is >> latitude >> longitude >>coor_sys;
	latitude *=degree;
	longitude *=degree;
	theMotherField->SetWorldCenterPosition(latitude,longitude,coor_sys);	
  }
  
  
  
  else if( command ==  SwitchOnFieldCmd) theMotherField->SwitchOn(); 
   
  else if( command ==  SwitchOffFieldCmd) theMotherField->SwitchOff();  
  
  else result = false;
  */
  return;

}












