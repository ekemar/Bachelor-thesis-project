#include"EarthMagneticFieldMessenger.hh"
#include"EarthMagneticField.hh"
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

EarthMagneticFieldMessenger::EarthMagneticFieldMessenger (EarthMagneticField* aField )
:PlanetMagneticFieldMessenger(aField) 
{ 
 theField = aField;
 G4String cmd_name;
 G4String guidance;
 
 
 
 /*SetTiltedDipoleParameterFromIGRFCmd= 
     new G4UIcmdWithoutParameter("/PLANETOCOS/BFIELD/SetNonshiftedGeodipoleFromIGRF",this);
 SetTiltedDipoleParameterFromIGRFCmd
            ->SetGuidance("Compute the parameters of the tilted geomagnetic dipole from the IGRF coefficients");
 SetTiltedDipoleParameterFromIGRFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
 SetEccentricDipoleParameterFromIGRFCmd= 
    new G4UIcmdWithoutParameter("/PLANETOCOS/BFIELD/SetShiftedGeodipoleFromIGRF",this);
 SetEccentricDipoleParameterFromIGRFCmd
     ->SetGuidance("Compute the parameters of the eccentric geomagnetic dipole from the IGRF coefficients");
 SetEccentricDipoleParameterFromIGRFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
 
 //Magnetic activity command
 ///////////////////////////
 
  cmd_name ="/PLANETOCOS/BFIELD/SetIopt";
  SetIoptCmd = new G4UIcmdWithAnInteger(cmd_name,this);
  SetIoptCmd->SetGuidance("Set the Iopt (kp level) parameter for the TSY89 model");
  SetIoptCmd->SetParameterName("Iopt",false);
  SetIoptCmd->SetRange("Iopt >0 Iopt<8");
  SetIoptCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name ="/PLANETOCOS/BFIELD/SetPdyn";
  SetPdynCmd = new G4UIcmdWithADouble(cmd_name,this);
  SetPdynCmd->SetGuidance("Set the Pdyn (sw dynamic pressure in nPa) parameter for TSY96, TSY2001, TSY2004");
  SetPdynCmd->SetParameterName("Pdyn",false);
  SetPdynCmd->SetRange("Pdyn > 0");
  SetPdynCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name ="/PLANETOCOS/BFIELD/SetTiltAngle";
  SetTiltAngleCmd = new G4UIcmdWithADouble(cmd_name,this);
  SetTiltAngleCmd->SetGuidance("Set the TiltAngle (rad) parameter for TSY96, TSY2001, TSY2004");
  SetTiltAngleCmd->SetParameterName("TiltAngle",false);
  SetTiltAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);   
  
  cmd_name ="/PLANETOCOS/BFIELD/SetDst";
  SetDstCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetDstCmd->SetGuidance("Set the Dst parameter for the Tsy96, TSY2001, TSY2004 models");
  SetDstCmd->SetParameterName("Dst",false);
  SetDstCmd->SetUnitCategory("Magnetic flux density");
  SetDstCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  cmd_name ="/PLANETOCOS/BFIELD/SetImfBy";
  SetImfyCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetImfyCmd->SetGuidance("Set the ImfY parameter for the Tsy96, TSY2001, TSY2004 models");
  SetImfyCmd->SetParameterName("Imfy",false);
  SetImfyCmd->SetUnitCategory("Magnetic flux density");
  SetImfyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name ="/PLANETOCOS/BFIELD/SetImfBz";
  SetImfzCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetImfzCmd->SetGuidance("Set the Imfz parameter for the Tsy96, TSY2001, TSY2004 models");
  SetImfzCmd->SetParameterName("Imfz",false);
  SetImfzCmd->SetUnitCategory("Magnetic flux density");
  SetImfzCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetG1";
  SetG1Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetG1Cmd->SetGuidance("Set the the G1  parameter for the TSY2001 model");
  SetG1Cmd->SetParameterName("G1",false);
  SetG1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name = "/PLANETOCOS/BFIELD/SetG2";
  SetG2Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetG2Cmd->SetGuidance("Set the G2 parameter for the TSY2001 model");
  SetG2Cmd->SetParameterName("G2",false);
  SetG2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name = "/PLANETOCOS/BFIELD/ReadTSY2001Parameters"; 
  ReadTSY2001ParameterCmd = new G4UIcmdWithAString(cmd_name,this);
  guidance = "Read in a file the values of the magnetic activity and solar ";
  guidance +="wind parameters in function of time (TSY2001)";
  ReadTSY2001ParameterCmd->SetGuidance(guidance);
  ReadTSY2001ParameterCmd->SetParameterName("filename",true);
  ReadTSY2001ParameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  cmd_name = "/PLANETOCOS/BFIELD/PrintTSY2001Parameters";
  PrintTSY2001ParameterCmd=  new G4UIcmdWithoutParameter(cmd_name,this);
  guidance = "Print the actual magnetic activity and solar wind parameters ";
  guidance +="used in  the Tsyganenko 2001 model";
  PrintTSY2001ParameterCmd->SetGuidance(guidance);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetW1";
  SetW1Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW1Cmd->SetGuidance("Set the the W1  parameter for the TSY2004 model");
  SetW1Cmd->SetParameterName("W1",false);
  SetW1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);   

  cmd_name = "/PLANETOCOS/BFIELD/SetW2";
  SetW2Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW2Cmd->SetGuidance("Set the the W2  parameter for the TSY2004 model");
  SetW2Cmd->SetParameterName("W2",false);
  SetW2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  cmd_name = "/PLANETOCOS/BFIELD/SetW3";
  SetW3Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW3Cmd->SetGuidance("Set the the W3  parameter for the TSY2004 model");
  SetW3Cmd->SetParameterName("W3",false);
  SetW3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  cmd_name = "/PLANETOCOS/BFIELD/SetW4";
  SetW4Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW4Cmd->SetGuidance("Set the the W4  parameter for the TSY2004 model");
  SetW4Cmd->SetParameterName("W4",false);
  SetW4Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  cmd_name = "/PLANETOCOS/BFIELD/SetW5";
  SetW5Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW5Cmd->SetGuidance("Set the the W5  parameter for the TSY2004 model");
  SetW5Cmd->SetParameterName("W5",false);
  SetW5Cmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  cmd_name = "/PLANETOCOS/BFIELD/SetW6";
  SetW6Cmd = new G4UIcmdWithADouble(cmd_name,this);
  SetW6Cmd->SetGuidance("Set the the W6  parameter for the TSY2004 model");
  SetW6Cmd->SetParameterName("W6",false);
  SetW6Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);   
  
  cmd_name = "/PLANETOCOS/BFIELD/ReadTSY2004Parameters"; 
  ReadTSY2004ParameterCmd = new G4UIcmdWithAString(cmd_name,this);
  guidance = "Read in a file the values of the magnetic activity and solar ";
  guidance +="wind parameters in function of time (TSY2004)";
  ReadTSY2004ParameterCmd->SetGuidance(guidance);
  ReadTSY2004ParameterCmd->SetParameterName("filename",true);
  ReadTSY2004ParameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  cmd_name = "/PLANETOCOS/BFIELD/PrintTSY2004Parameters";
  PrintTSY2004ParameterCmd=  new G4UIcmdWithoutParameter(cmd_name,this);
  guidance = "Print the actual magnetic activity and solar wind parameters ";
  guidance +="used in  the Tsyganenko 2004 model";
  PrintTSY2004ParameterCmd->SetGuidance(guidance);  
 
}
////////////////////////////////////////////////////////////////////////////////
//
EarthMagneticFieldMessenger::~EarthMagneticFieldMessenger()
{  delete SetIoptCmd;
   delete SetPdynCmd;
   delete SetTiltAngleCmd;
   delete SetDstCmd;
   delete SetImfyCmd;
   delete SetImfzCmd;
   delete SetG1Cmd;
   delete SetG2Cmd;
   delete ReadTSY2001ParameterCmd;
   delete SetW1Cmd;   
   delete SetW2Cmd;   
   delete SetW3Cmd;   
   delete SetW4Cmd;   
   delete SetW5Cmd;   
   delete SetW6Cmd;   
   delete ReadTSY2004ParameterCmd;
} 
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticFieldMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if (SetMotherNewValue(command, newValues)) return;
 
 // magnetic activity command
  
  
  /*else if( command == SetTiltedDipoleParameterFromIGRFCmd)
                   theField-> SetTiltedDipoleParameterFromIGRF(); 
  
  else if( command == SetEccentricDipoleParameterFromIGRFCmd)
               theField-> SetEccentricDipoleParameterFromIGRF(); */
  
  else if( command == SetIoptCmd)
            theField->SetIopt(SetIoptCmd->GetNewIntValue(newValues));   	    
  
  
  else if (command == SetPdynCmd) 
            theField->SetPdyn(SetPdynCmd->GetNewDoubleValue(newValues));   
  
  else if (command == SetTiltAngleCmd) 
            theField->SetTiltAngle(SetTiltAngleCmd->GetNewDoubleValue(newValues));  
  
  else if (command == SetDstCmd) 
            theField->SetDst(SetDstCmd->GetNewDoubleValue(newValues));
  
  else if (command == SetImfyCmd) 
            theField->SetImfy(SetImfyCmd->GetNewDoubleValue(newValues));   
  
  else if (command == SetImfzCmd) 
            theField->SetImfz(SetImfzCmd->GetNewDoubleValue(newValues));

  else if (command == SetG1Cmd) 
            theField->SetG1(SetG1Cmd->GetNewDoubleValue(newValues));   
  
  else if (command == SetG2Cmd) 
            theField->SetG2(SetG2Cmd->GetNewDoubleValue(newValues));
  
  else if( command == ReadTSY2001ParameterCmd)
            theField->ReadTSY2001Parameter(newValues);	       	  
  
  else if( command == PrintTSY2001ParameterCmd)
                                theField->PrintStormParameterTSY2001(); 
  
  else if (command == SetW1Cmd) 
            theField->SetW1(SetW1Cmd->GetNewDoubleValue(newValues));     

  else if (command == SetW2Cmd) 
            theField->SetW2(SetW2Cmd->GetNewDoubleValue(newValues)); 
  
  else if (command == SetW3Cmd) 
            theField->SetW3(SetW3Cmd->GetNewDoubleValue(newValues)); 

  else if (command == SetW4Cmd) 
            theField->SetW4(SetW4Cmd->GetNewDoubleValue(newValues)); 

  else if (command == SetW5Cmd) 
            theField->SetW5(SetW5Cmd->GetNewDoubleValue(newValues)); 

  else if (command == SetW6Cmd) 
            theField->SetW6(SetW6Cmd->GetNewDoubleValue(newValues));   
  
  else if( command == ReadTSY2004ParameterCmd)
            theField->ReadTSY2004Parameter(newValues);	       	  
  
  else if( command == PrintTSY2004ParameterCmd)
                                theField->PrintStormParameterTSY2004();    
}
