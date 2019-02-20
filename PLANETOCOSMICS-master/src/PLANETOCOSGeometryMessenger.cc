#include "PLANETOCOSGeometryMessenger.hh"
#include "PLANETOCOSGeometryConstruction.hh"
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
#include "EarthAtmosphereModel.hh"
#include "DateAndTime.hh"


PLANETOCOSGeometryMessenger::PLANETOCOSGeometryMessenger(PLANETOCOSGeometryConstruction * myDet)
:myGeometry(myDet)
{  
   //parameters
   G4UIparameter* month = new G4UIparameter("month",'i',false);
   month->SetParameterRange(" month > 0 && month < 13");
   
   //G4UIparameter* year = new G4UIparameter("year",'i',false);
   
   G4UIparameter* day = new G4UIparameter("day",'i',false);
   day->SetParameterRange(" day > 0 && day < 32");
   
   G4UIparameter* hour = new G4UIparameter("hour",'i',true);
   hour->SetParameterRange(" hour >= 0 && hour < 24");
   hour->SetDefaultValue(0);
   
   G4UIparameter* minute = new G4UIparameter("minute",'i',true);
   minute->SetParameterRange(" minute >= 0 && minute < 60");
   minute->SetDefaultValue(0);
   
   
   G4UIparameter* second = new G4UIparameter("second",'i',true);
   second->SetParameterRange(" second >= 0 && second < 60");
   second->SetDefaultValue(0);
   
   G4UIparameter* longitude = new G4UIparameter("longitude",'d',false);
   longitude->SetParameterRange(" longitude >= -360 && longitude <= 360");
   
   G4UIparameter* latitude = new G4UIparameter("latitude",'d',false);
   latitude->SetParameterRange(" latitude >= -90 && latitude <= 90");
   
   G4UIparameter* coord_sys = new G4UIparameter("coord_sys",'s',false);
   coord_sys->SetParameterCandidates("GEODETIC GEO GEOID");
   
   G4UIparameter* length_unit = new G4UIparameter("length_unit",'s',false);
   length_unit->SetParameterCandidates("m km");
   
   G4UIparameter* depth_unit = new G4UIparameter("depth_unit",'s',false);
   depth_unit->SetParameterCandidates("g/cm2 g/m2 kg/m2 kg/cm2"); 
   
   
   
   
   G4String guidance;
   G4String cmd_title;
  
   //command directories 
   atmocosmicsDir = new G4UIdirectory("/PLANETOCOS/");
   atmocosmicsDir->SetGuidance("Planetocosmics control"); 
  
   myGeometryDir = new G4UIdirectory("/PLANETOCOS/GEOMETRY/");
   myGeometryDir->SetGuidance("Geometry  control");
  
  
   myUserlimitDir = new G4UIdirectory("/PLANETOCOS/USERLIMIT/");
   myUserlimitDir->SetGuidance("User limit  control");
   
   //Geometry UI commands
   
   cmd_title="/PLANETOCOS/USERLIMIT/SetAtmoMaxStepLength";
   AtmosphereMaxStepLengthCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set the Atmosphere MaxStepLength";
   AtmosphereMaxStepLengthCmd->SetGuidance(guidance);
   AtmosphereMaxStepLengthCmd->SetParameterName("MaxStepLength",false);
   AtmosphereMaxStepLengthCmd->AvailableForStates(G4State_Idle);
   AtmosphereMaxStepLengthCmd->SetUnitCategory("Length");
   
   cmd_title="/PLANETOCOS/USERLIMIT/SetMagnetoMaxStepLength";
   MagnetosphereMaxStepLengthCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set the Magnetosphere MaxStepLength ";
   MagnetosphereMaxStepLengthCmd->SetGuidance(guidance);
   MagnetosphereMaxStepLengthCmd->SetParameterName("MaxStepLength",false);
   MagnetosphereMaxStepLengthCmd->AvailableForStates(G4State_Idle);
   MagnetosphereMaxStepLengthCmd->SetUnitCategory("Length");
   
   cmd_title="/PLANETOCOS/STOPCONDITION/SetMaxTrackLength";
   MagnetosphereMaxTrackLengthCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set the maximum track length.";
   MagnetosphereMaxTrackLengthCmd->SetGuidance(guidance);
   MagnetosphereMaxTrackLengthCmd->SetParameterName("MaxTrackLength",false);
   MagnetosphereMaxTrackLengthCmd->AvailableForStates(G4State_Idle);
   MagnetosphereMaxTrackLengthCmd->SetUnitCategory("Length");
   
   cmd_title="/PLANETOCOS/STOPCONDITION/SetMaxTrackDuration";
   MagnetosphereMaxTrackDurationCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set the maximum track duration.";
   MagnetosphereMaxTrackDurationCmd->SetGuidance(guidance);
   MagnetosphereMaxTrackDurationCmd->SetParameterName("MaxTrackDuration",false);
   MagnetosphereMaxTrackDurationCmd->AvailableForStates(G4State_Idle);
   MagnetosphereMaxTrackDurationCmd->SetUnitCategory("Time");
   
   
   
   cmd_title="/PLANETOCOS/GEOMETRY/SetAtmosphereTop";
   SetAtmosphereHmaxCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Defines the altitude of the top of the planet atmosphere";                         
   SetAtmosphereHmaxCmd->SetGuidance(guidance);
   SetAtmosphereHmaxCmd->SetUnitCandidates("km m");
   SetAtmosphereHmaxCmd->SetParameterName("MaxAltitude",false);
   SetAtmosphereHmaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   cmd_title = "/PLANETOCOS/GEOMETRY/SetGroundAltitude";
   SetAtmosphereHminCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance ="Defines the altitude of the bottom of the planet atmosphere";                          
   SetAtmosphereHminCmd->SetGuidance(guidance);
   SetAtmosphereHminCmd->SetUnitCandidates("km m");
   SetAtmosphereHminCmd->SetParameterName("MinAltitude",false);
   SetAtmosphereHminCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   cmd_title = "/PLANETOCOS/GEOMETRY/SetHeigthOfWorldAboveAtmosphere";
   SetMagnetosphereHCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance ="Defines the heigth of the worl volume above the atmosphere";                          
   SetMagnetosphereHCmd->SetGuidance(guidance);
   SetMagnetosphereHCmd->SetUnitCategory("Length");
   SetMagnetosphereHCmd->SetParameterName("HWorld",false);
   SetMagnetosphereHCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   
   cmd_title = "/PLANETOCOS/GEOMETRY/SetPlanetCoreThickness";
   SetPlanetHCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance ="Defines the surface thickness";                          
   SetPlanetHCmd->SetGuidance(guidance);
   SetPlanetHCmd->SetUnitCandidates("km m");
   SetPlanetHCmd->SetParameterName("CoreThickess",false);
   SetPlanetHCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   cmd_title = "/PLANETOCOS/GEOMETRY/verbose";
   SetVerbosityCmd = new  G4UIcmdWithAnInteger(cmd_title,this);
   guidance ="If the geometry verbosity > 0 the structure ";
   guidance +=" of the atmosphere is printed when the geometry is built."; 
   SetVerbosityCmd->SetGuidance(guidance);
   SetVerbosityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   cmd_title = "/PLANETOCOS/GEOMETRY/SetMaxLayerThickness";
   SetMaxThicknessCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Defines the maximum thickness of an atmospheric layer";                         
   SetMaxThicknessCmd->SetGuidance(guidance);
   SetMaxThicknessCmd->SetUnitCandidates("km m");
   SetMaxThicknessCmd->SetParameterName("MaxThickness",false);
   SetMaxThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
   cmd_title = "/PLANETOCOS/GEOMETRY/SetMinLayerThickness";
   SetMinThicknessCmd  =  new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Defines the minimum thickness of an atmospheric layer";                         
   SetMinThicknessCmd->SetGuidance(guidance);
   SetMinThicknessCmd->SetUnitCandidates("m km");
   SetMinThicknessCmd->SetParameterName("MinThickness",false);
   SetMinThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   cmd_title = "/PLANETOCOS/GEOMETRY/SetLayerLength";
   SetHalfLengthCmd  = new  G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance="Defines the length of the atmospheric layers";                         
   SetHalfLengthCmd->SetGuidance(guidance);
   SetHalfLengthCmd->SetUnitCategory("Length");
   SetHalfLengthCmd->SetParameterName("HalfLength",false);
   SetHalfLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   cmd_title ="/PLANETOCOS/GEOMETRY/SetPercentOfDepth";
   SetDepthPercentCmd  =  new  G4UIcmdWithADouble(cmd_title,this);
   guidance="Defines the percent of total depth contained in one layer ";
   SetDepthPercentCmd->SetGuidance(guidance);
   SetDepthPercentCmd->SetParameterName("DepthPercent",false);
   SetDepthPercentCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   cmd_title = "/PLANETOCOS/GEOMETRY/SetType";
   SetGeometryTypeCmd = new G4UIcmdWithAString(cmd_title,this);
   guidance = "Define the type of geometry: spherical or flat";
   SetGeometryTypeCmd->SetGuidance(guidance);
   SetGeometryTypeCmd->SetParameterName("GeometryType",false);
   SetGeometryTypeCmd->SetCandidates("Flat flat FLAT spherical Spherical SPHERICAL");
   SetGeometryTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  /* cmd_title = "/PLANETOCOS/GEOMETRY/Construct";
   UpdateGeometryCmd = new  G4UIcmdWithoutParameter(cmd_title,this);
   UpdateGeometryCmd->SetGuidance("Construct the geometry. ");
   guidance = "Should be called before the definition of the pysics list ";
   guidance += " and the /run initialization cmd";
   UpdateGeometryCmd->SetGuidance(guidance); 
   UpdateGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
   cmd_title = "/PLANETOCOS/GEOMETRY/RemoveAllDetectors";
   RemoveAllDetectorsCmd =  new  G4UIcmdWithoutParameter(cmd_title,this);
   guidance = "Clear the vector of altitudes and depths at which shower ";
   guidance += " particle flux can be detected during the simulation."; 
   RemoveAllDetectorsCmd->SetGuidance(guidance);
   RemoveAllDetectorsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   cmd_title = "/PLANETOCOS/GEOMETRY/DetectorAtAltitude";
   AddDetectorAtAltitudeCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set a new detector at a given altitude ";                         
   AddDetectorAtAltitudeCmd->SetGuidance(guidance);
   AddDetectorAtAltitudeCmd->SetUnitCandidates("km m");
   AddDetectorAtAltitudeCmd->SetParameterName("Altitude",false);
   AddDetectorAtAltitudeCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
        
   cmd_title = "/PLANETOCOS/GEOMETRY/DetectorAtDepth";
   AddDetectorAtDepthCmd = new G4UIcmdWithADoubleAndUnit(cmd_title,this);
   guidance = "Set a new detector at a given depth ";                         
   AddDetectorAtDepthCmd->SetGuidance(guidance);
   AddDetectorAtDepthCmd->SetUnitCategory("Depth");
   AddDetectorAtDepthCmd->SetParameterName("Depth",false);
   AddDetectorAtDepthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  
  
  //Atmosphere model cmd
  if (myGeometry->GetAtmosphericModel()){
   	SetAtmosphereModelCmd = new
         	G4UIcmdWithAString("/PLANETOCOS/GEOMETRY/SetAtmosphereModel",this);
   	guidance ="Define the atmospheric composition model"; 
   	SetAtmosphereModelCmd->SetGuidance(guidance);
   	SetAtmosphereModelCmd->SetParameterName("AtmosphereModel",false);
   	std::vector<G4String > ListOfModels = myGeometry->GetAtmosphericModel()
						->GetListOfAtmosphericModels();
   	G4String candidates= ListOfModels[0];
	for (unsigned int i=1; i<ListOfModels.size();i++){
		candidates+=" "+ListOfModels[i];
	
	}

   	SetAtmosphereModelCmd->SetCandidates(candidates);
   	SetAtmosphereModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   	
	
	cmd_title = "/PLANETOCOS/GEOMETRY/ReadAtmosphereCompositionTable";
   	ReadAtmosphereModelCmd = new G4UIcmdWithAString(cmd_title,this);
   	guidance = "read a table where the atmospheric composition is defined";
   	ReadAtmosphereModelCmd->SetGuidance(guidance);
   	ReadAtmosphereModelCmd->SetParameterName("file_name",false);
   	ReadAtmosphereModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
        cmd_title = "/PLANETOCOS/GEOMETRY/SetConsiderAtmosphere";
        SetConsiderAtmosphereCmd = new  G4UIcmdWithABool(cmd_title,this);
        guidance ="If true the atmosphere (if any) will be considered. \n";
        guidance +="If false the atmosphere is not considered. "; 
        SetConsiderAtmosphereCmd->SetGuidance(guidance);
        SetConsiderAtmosphereCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
	
	
	cmd_title = "/PLANETOCOS/GEOMETRY/SetDetectInMiddleOfAtmosphereMode";
        SetDetectInMiddleOfAtmosphereModeCmd = new  G4UIcmdWithABool(cmd_title,this);
        guidance ="If true the atmosphere (if any) will be considered. \n";
        guidance +="If false the atmosphere is not considered. "; 
        SetDetectInMiddleOfAtmosphereModeCmd->SetGuidance(guidance);
        SetDetectInMiddleOfAtmosphereModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
   
   
   }
    
    
    cmd_title = "/PLANETOCOS/GEOMETRY/SetDetectionBelowSoilLayers";
    SetDetectionBelowSoilLayersCmd = new  G4UIcmdWithABool(cmd_title,this);
    guidance ="If true the flux will be detected below detection layer \n"; 
    SetDetectionBelowSoilLayersCmd->SetGuidance(guidance);
    SetDetectionBelowSoilLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  
   
  
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSGeometryMessenger::~PLANETOCOSGeometryMessenger()
{ delete myGeometryDir;
  delete AtmosphereMaxStepLengthCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ if( command == AtmosphereMaxStepLengthCmd ) 
  	myGeometry->SetAtmosphereMaxStepLength(AtmosphereMaxStepLengthCmd
                              ->GetNewDoubleValue(newValues));
  else if (command == MagnetosphereMaxStepLengthCmd ) 
  	myGeometry->SetMagnetosphereMaxStepLength(MagnetosphereMaxStepLengthCmd
                              ->GetNewDoubleValue(newValues));
  else if (command == MagnetosphereMaxTrackLengthCmd ) 
  	myGeometry->SetMagnetosphereMaxTrackLength(MagnetosphereMaxTrackLengthCmd
                              ->GetNewDoubleValue(newValues));
  else if (command == MagnetosphereMaxTrackDurationCmd ) 
  	myGeometry->SetMagnetosphereMaxTrackDuration(MagnetosphereMaxTrackDurationCmd
                              ->GetNewDoubleValue(newValues));			      			      
  else if (command == SetAtmosphereHmaxCmd)
     	myGeometry->SetAtmosphereHmax(SetAtmosphereHmaxCmd->GetNewDoubleValue(newValues)); 
  
  else if (command == SetAtmosphereHminCmd)
     	myGeometry->SetAtmosphereHmin(SetAtmosphereHminCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetMagnetosphereHCmd)
     	myGeometry->SetMagnetosphereH(SetMagnetosphereHCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetPlanetHCmd)
     	myGeometry->SetPlanetH(SetPlanetHCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetVerbosityCmd)
  	myGeometry->SetGeometryVerbosity(
	                        SetVerbosityCmd->GetNewIntValue(newValues));
  
  else if (command == SetMaxThicknessCmd)
     	myGeometry->SetMaxThickness(SetMaxThicknessCmd->GetNewDoubleValue(newValues)); 
  
  else if (command == SetMinThicknessCmd)
     	myGeometry->SetMinThickness(SetMinThicknessCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetHalfLengthCmd)
     	myGeometry->SetHalfLength(SetHalfLengthCmd->GetNewDoubleValue(newValues)/2.);
  
  else if (command == SetDepthPercentCmd)
     	myGeometry->SetDepthPercent(SetDepthPercentCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetGeometryTypeCmd){
  	G4String type=newValues; 
      	type.toUpper();
      	myGeometry->SetGeometryType(type);
  }
  
  //else if (command == UpdateGeometryCmd) myGeometry->UpdateGeometry();
  
  else if (command == RemoveAllDetectorsCmd) myGeometry->RemoveAllDetectors();
  
  else if (command == AddDetectorAtAltitudeCmd)
     	myGeometry->AddDetectorAtAltitude
                    (AddDetectorAtAltitudeCmd->GetNewDoubleValue(newValues));
  
  else if (command == AddDetectorAtDepthCmd)
     	myGeometry->AddDetectorAtDepth
                    (AddDetectorAtDepthCmd->GetNewDoubleValue(newValues));
  else if (command ==  SetDetectionBelowSoilLayersCmd){
		myGeometry->SetDetectionBelowSoilLayers(
				SetDetectionBelowSoilLayersCmd->GetNewBoolValue(newValues));
		
  }		    
  else if (myGeometry->GetAtmosphericModel()){		    
  	if (command == SetAtmosphereModelCmd)
     		myGeometry->GetAtmosphericModel()->SetAtmosphericModel(newValues);
  
  	else if (command == ReadAtmosphereModelCmd)
     		myGeometry->GetAtmosphericModel()->ReadAtmosphereComposition(newValues);
  	else if (command ==  SetConsiderAtmosphereCmd){
		myGeometry->SetConsiderAtmosphere(
				SetConsiderAtmosphereCmd->GetNewBoolValue(newValues));
		
	}
	else if (command ==  SetDetectInMiddleOfAtmosphereModeCmd){
		myGeometry->SetDetectionInAtmosphereInMiddleOfALayer(
				SetDetectInMiddleOfAtmosphereModeCmd->GetNewBoolValue(newValues));
		
	}	
		
  }
  
  
  
  
}
