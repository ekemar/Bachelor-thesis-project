#include "EarthAtmosphereMessenger.hh"
#include "EarthAtmosphereModel.hh"
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
#include <sstream>


EarthAtmosphereMessenger::EarthAtmosphereMessenger(EarthAtmosphereModel * myMod)
:myModel(myMod)
{  
   //parameters
   G4UIparameter* month = new G4UIparameter("month",'i',false);
   month->SetParameterRange(" month > 0 && month < 13");
   
   G4UIparameter* year = new G4UIparameter("year",'i',false);
   
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
   coord_sys->SetParameterCandidates("GEODETIC GEO ");
   
   G4UIparameter* length_unit = new G4UIparameter("length_unit",'s',false);
   length_unit->SetParameterCandidates("m km");
   
   G4UIparameter* depth_unit = new G4UIparameter("depth_unit",'s',false);
   depth_unit->SetParameterCandidates("g/cm2 g/m2 kg/m2 kg/cm2"); 
   
   
   
   
   G4String guidance;
   G4String cmd_title;
  
  
  
  
   //Earth atmosphere model command
  
   SetAtmosphereReferenceDateCmd = new 
           G4UIcommand("/PLANETOCOS/GEOMETRY/SetReferenceDate",this);
   guidance ="Define the reference date for the atmospheric model";
   SetAtmosphereReferenceDateCmd->SetGuidance(guidance);	   
   SetAtmosphereReferenceDateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   SetAtmosphereReferenceDateCmd->SetParameter(year);
   SetAtmosphereReferenceDateCmd->SetParameter(month);
   SetAtmosphereReferenceDateCmd->SetParameter(day);
   SetAtmosphereReferenceDateCmd->SetParameter(hour);
   SetAtmosphereReferenceDateCmd->SetParameter(minute);
   SetAtmosphereReferenceDateCmd->SetParameter(second);
  
  
   SetApCmd= 
       new  G4UIcmdWithADouble("/PLANETOCOS/GEOMETRY/SetAp",this);
   guidance="Defines the dayly Ap value used in MSIS atmospheric models ";
   SetApCmd->SetGuidance(guidance);
   SetApCmd->SetParameterName("Ap",false);
   SetApCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   SetF107Cmd= 
       new  G4UIcmdWithADouble("/PLANETOCOS/GEOMETRY/SetF107",this);
   guidance="Defines the dayly F107 value used in MSIS atmospheric models ";
   SetF107Cmd->SetGuidance(guidance);
   SetF107Cmd->SetParameterName("F107",false);
   SetF107Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   SetF107ACmd= 
       new  G4UIcmdWithADouble("/PLANETOCOS/GEOMETRY/SetF107A",this);
   guidance="Defines the dayly F107A value used in MSIS atmospheric models ";
   SetF107ACmd->SetGuidance(guidance);
   SetF107ACmd->SetParameterName("F107A",false);
   SetF107ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   SetPositionCmd = new 
           G4UIcommand("/PLANETOCOS/GEOMETRY/SetReferencePosition",this);
   guidance ="Define the position for the atmosphere model";
   SetPositionCmd->SetGuidance(guidance);	   
   SetPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   SetPositionCmd->SetParameter(latitude);
   SetPositionCmd->SetParameter(longitude);
 
   
}
////////////////////////////////////////////////////////////////////////////////
//
EarthAtmosphereMessenger::~EarthAtmosphereMessenger()
{ ;
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthAtmosphereMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ if (command == SetAtmosphereReferenceDateCmd){ 
       	const char* paramString=newValues;
        G4int year, month, day, hour, minute, second;
	std::istringstream is((char*)paramString);
        is >> year >> month >> day >> hour >> minute >> second;
	DateAndTime date = DateAndTime(year,month,day,hour,minute,second);
	myModel->SetReferenceDate(date);
  }
  else if (command == SetApCmd)
     	myModel->SetAp(SetApCmd->GetNewDoubleValue(newValues));    
  else if (command == SetF107Cmd)
      	myModel->SetF107(SetF107Cmd->GetNewDoubleValue(newValues));   
  else if (command == SetF107ACmd)
     	myModel->SetF107A(SetF107ACmd->GetNewDoubleValue(newValues));
  else if (command == SetPositionCmd){ 
    	const char* paramString=newValues;
        G4double longitude, latitude;
	std::istringstream is((char*)paramString);
        is >> latitude>>longitude;
	myModel->SetGEODETICLongitude(longitude);
	myModel->SetGEODETICLatitude(latitude);
  }
	
  
}
