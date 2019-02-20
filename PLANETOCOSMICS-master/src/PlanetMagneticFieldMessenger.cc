#include"PlanetMagneticField.hh"
#include"PlanetMagneticFieldMessenger.hh"
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

PlanetMagneticFieldMessenger::PlanetMagneticFieldMessenger (PlanetMagneticField* aField )
{ 
  theMotherField = aField;

  G4String candidates;
  G4String guidance;
  G4String cmd_name;
 
//directories
///////////////
 
  IntegrationDir= new G4UIdirectory("/PLANETOCOS/INTEGRATION/");
  IntegrationDir->SetGuidance("Numerical integration control");
  
  MagnetoDir = new G4UIdirectory("/PLANETOCOS/BFIELD/");
  MagnetoDir->SetGuidance("Magnetic field control"); 
 
//parameters
/////////////
  
  candidates ="PLA PLAG PSO PSEQ ";
  if (theMotherField->GetHasAGlobalField())
  			candidates +="PMAG PSMAG PSM ";
  G4UIparameter* coord_sys_param = 
               new G4UIparameter("Space coordinate system",'s',false);
  coord_sys_param->SetParameterCandidates(candidates);
  
  G4UIparameter* length_unit_param = new G4UIparameter("Length Unit",'s',false);
  length_unit_param->SetParameterCandidates("Rplanet km m");
  
  G4UIparameter* angle_unit_param = new G4UIparameter("Angle Unit",'s',false);
  angle_unit_param->SetParameterCandidates("degree deg rad radian mrad milliradian");
  
  G4UIparameter* altitude_param = new G4UIparameter("Altitude",'d',false);
  G4UIparameter* dalt_param = new G4UIparameter("dalt",'d',false);
  G4UIparameter* nAlt_param = new G4UIparameter("number of altitude",'i',false);
  
  G4UIparameter* x_param = new G4UIparameter("x",'d',false);
  G4UIparameter* dx_param = new G4UIparameter("dx",'d',false);
  G4UIparameter* nx_param = new G4UIparameter("nx",'i',false);
  
  G4UIparameter* y_param = new G4UIparameter("y",'d',false);
  G4UIparameter* dy_param = new G4UIparameter("dy",'d',false);
  G4UIparameter* ny_param = new G4UIparameter("ny",'i',false);
  
  G4UIparameter* nz_param = new G4UIparameter("nz",'i',false);
  
  G4UIparameter* longitude_param = new G4UIparameter("Longitude",'d',false);
  G4UIparameter* dlong_param = new G4UIparameter("delta longitude",'d',false);
  G4UIparameter* nLong_param = new G4UIparameter("number of longitude",'i',false);
  
  G4UIparameter* latitude_param = new G4UIparameter("Latitude",'d',false);
  G4UIparameter* dlat_param = new G4UIparameter("delta latitude",'d',false);
  G4UIparameter* nLat_param = new G4UIparameter("number of altitude",'i',false);
  
  
  G4UIparameter* theta_param = new G4UIparameter("Theta",'d',false);
  theta_param->SetParameterRange("Theta >=0 ");
  G4UIparameter* phi_param = new G4UIparameter("Phi",'d',false);
  
  
  G4UIparameter* output_file_param = new G4UIparameter("OutputFile name",'s',false);
  
  G4UIparameter* int_model_parameter = new G4UIparameter("Internal Model",'s',false);
  candidates ="";
  std::vector< G4String > ListOfModels = theMotherField->GetListOfInternalFieldModels();
  for (unsigned int i=0;i<ListOfModels.size();i++){
  	candidates+=" "+ListOfModels[i];
  }
  int_model_parameter->SetParameterCandidates(candidates);
  
  G4UIparameter* ext_model_parameter = new G4UIparameter("External Model",'s',true);
  candidates ="";
  ListOfModels = theMotherField->GetListOfExternalFieldModels();
  for (unsigned int i=0;i<ListOfModels.size();i++){
  	candidates+=" "+ListOfModels[i];
  }
  ext_model_parameter->SetParameterCandidates(candidates);
  ext_model_parameter->SetDefaultValue("NOFIELD");
  
  
  G4UIparameter* year_param = new G4UIparameter("Year",'i',false);
 
  G4UIparameter* month_param  = new G4UIparameter("Month",'i',false);
  month_param->SetParameterRange("Month <13 && Month>0");
  
  G4UIparameter* day_param = new G4UIparameter("Day",'i',false);
  day_param->SetParameterRange("Day <32 && Day>0");
  
  G4UIparameter* hour_param = new G4UIparameter("Hour",'i',true);
  hour_param->SetDefaultValue("0");
  hour_param->SetParameterRange("Hour<24");
  
  G4UIparameter* minute_param = new G4UIparameter("Minute",'i',true);
  minute_param->SetDefaultValue("0");
  minute_param->SetParameterRange("Minute<60");
  
  G4UIparameter* second_param = new G4UIparameter("Second",'i',true);
  second_param->SetDefaultValue("0");
  second_param->SetParameterRange("Second<60");
    
// Integration control commands
/////////////////////////////////

  cmd_name = "/PLANETOCOS/INTEGRATION/SetPrecision";
  SetEpsilonCmd = new G4UIcmdWithADouble(cmd_name,this);
  SetEpsilonCmd->SetGuidance("Set the relative precision for integration");
  SetEpsilonCmd->SetParameterName("precision",false);
  SetEpsilonCmd->SetRange("precision <= 1.e-3  && precision >=1.e-8");
  SetEpsilonCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name = "/PLANETOCOS/INTEGRATION/SetG4MaxStep";
  SetG4DeltaChordCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetG4DeltaChordCmd->SetGuidance("Set the maximal integrating step allowed for G4 integration");
  SetG4DeltaChordCmd->SetParameterName("DeltaChord",false);
  SetG4DeltaChordCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetG4DeltaChordCmd->SetUnitCategory("Length");
  
  cmd_name = "/PLANETOCOS/INTEGRATION/SetDeltaIntersection"; 
  SetDeltaIntersectionCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetDeltaIntersectionCmd->SetGuidance("Set the precision for crossing boundary");
  SetDeltaIntersectionCmd->SetParameterName("DeltaIntersection",false);
  SetDeltaIntersectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetDeltaIntersectionCmd->SetUnitCategory("Length");
 
  
  cmd_name = "/PLANETOCOS/INTEGRATION/SetDefaultIntegrationParameters";
  ResetIntegrationParametersCmd=  new G4UIcmdWithoutParameter(cmd_name,this);
  ResetIntegrationParametersCmd->SetGuidance("Set the integration parameters to their default values");
  ResetIntegrationParametersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/INTEGRATION/SetStepperModel";
  SetStepperCmd = new G4UIcmdWithAString(cmd_name,this);
  SetStepperCmd->SetGuidance("Set the stepper model for the G4Integration method ");
  SetStepperCmd->SetParameterName("choice",true);
  SetStepperCmd->SetDefaultValue("CashKarpRKF45");
  SetStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  candidates = "ExplicitEuler ImplicitEuler SimpleRunge ClassicalRK4 ";
  candidates += "CashKarpRKF45 RKG3_Stepper";
  SetStepperCmd->SetCandidates(candidates);  
  
  // magnetic field model commands
  /////////////////////////////////////
  
  cmd_name = "/PLANETOCOS/BFIELD/SetTimeOfB";
  SetTimeOfBCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetTimeOfBCmd->SetGuidance("Set the time from the start date at which B is computed ");
  SetTimeOfBCmd->SetParameterName("time",false);
  SetTimeOfBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetTimeOfBCmd->SetUnitCategory("Time");
  
  cmd_name = "/PLANETOCOS/BFIELD/SetStartDate";
  SetStartDateCmd = new G4UIcommand(cmd_name,this);
  SetStartDateCmd->SetGuidance("Set the start reference date"); 
  SetStartDateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetStartDateCmd->SetParameter(year_param);
  SetStartDateCmd->SetParameter(month_param);
  SetStartDateCmd->SetParameter(day_param);
  SetStartDateCmd->SetParameter(hour_param);
  SetStartDateCmd->SetParameter(minute_param);
  SetStartDateCmd->SetParameter(second_param);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetNmaxForGauss";
  Setnmax_GaussCmd = new G4UIcmdWithAnInteger(cmd_name,this);
  guidance ="Set the maximum order of the spherical harmonics coefficients";
  guidance +=" used for computing the Gauss model"; 
  Setnmax_GaussCmd->SetGuidance(guidance);
  Setnmax_GaussCmd->SetParameterName("Nmax",false);
  Setnmax_GaussCmd->SetRange("Nmax <100 && Nmax >0");
  Setnmax_GaussCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetInternalFieldModel";
  SetInternalFieldCmd = new G4UIcmdWithAString(cmd_name,this);
  SetInternalFieldCmd->SetGuidance("Set the model of the internal magnetic field");
  SetInternalFieldCmd->SetParameterName("InternalModel",true);
  candidates ="";
  ListOfModels = theMotherField->GetListOfInternalFieldModels();
  for (unsigned int i=0;i<ListOfModels.size();i++){
  	candidates+=" "+ListOfModels[i];
  }
  SetInternalFieldCmd->SetCandidates(candidates);
  SetInternalFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name = "/PLANETOCOS/BFIELD/SetExternalFieldModel";
  SetExternalFieldCmd = new G4UIcmdWithAString(cmd_name,this);
  SetExternalFieldCmd->SetGuidance("Set the model of the internal magnetic field");
  SetExternalFieldCmd->SetParameterName("ExternalModel",true);
  candidates ="";
  std::vector< G4String > ListOfExternalFieldModels = theMotherField->GetListOfExternalFieldModels();
  for (unsigned int i=0;i<ListOfExternalFieldModels.size();i++){
  	candidates+=" "+ListOfExternalFieldModels[i];
  }
  SetExternalFieldCmd->SetCandidates(candidates);
  SetExternalFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 
  cmd_name = "/PLANETOCOS/BFIELD/SetMagnetopauseModel";
  SetMagnetopauseModelCmd = new G4UIcmdWithAString(cmd_name,this);
  SetMagnetopauseModelCmd->SetGuidance("Set the model of the internal magnetic field");
  SetMagnetopauseModelCmd->SetParameterName("ExternalModel",true);
  candidates ="";
  std::vector< G4String > ListOfMagnetopauseModels = theMotherField->GetListOfMagnetopauseModels();
  for (unsigned int i=0;i<ListOfMagnetopauseModels.size();i++){
  	candidates+=" "+ListOfMagnetopauseModels[i];
  }
  SetMagnetopauseModelCmd->SetCandidates(candidates);
  SetMagnetopauseModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  
  
  
  cmd_name = "/PLANETOCOS/BFIELD/SetDipoleB0";
  SetDipoleB0Cmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetDipoleB0Cmd->SetGuidance("Set the magnetic moment B0 of the magnetic dipole");
  SetDipoleB0Cmd->SetParameterName("B0",false);
  SetDipoleB0Cmd->SetRange("B0 > 0.");
  SetDipoleB0Cmd->SetUnitCategory("Magnetic flux density");
  SetDipoleB0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetDipoleCenter";
  SetDipoleShiftCmd = new G4UIcmdWith3VectorAndUnit(cmd_name,this);
  SetDipoleShiftCmd->SetGuidance("Define the center  of the plamagnetic dipole in PLA coordinates");
  SetDipoleShiftCmd->SetParameterName("X","Y","Z",false);
  SetDipoleShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetDipoleShiftCmd->SetUnitCategory("Length");
  
  SetDipoleAxisCmd = new G4UIcommand("/PLANETOCOS/BFIELD/SetDipoleAxis",this);
  SetDipoleAxisCmd->SetGuidance("Set the  axis of the global magnetic dipole");
  SetDipoleAxisCmd
    ->SetGuidance("[usage] /PLANETOCOS/BFIELD/SetDipoleAxis   Theta Phi Unit" ); 
  SetDipoleAxisCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetDipoleAxisCmd->SetParameter(theta_param);
  SetDipoleAxisCmd->SetParameter(phi_param);
  SetDipoleAxisCmd->SetParameter(angle_unit_param);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetRadiusMagnetosphere";
  SetRadiusMagnetosphereCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetRadiusMagnetosphereCmd->SetGuidance("Set the radius of the spherical magnetosphere");
  SetRadiusMagnetosphereCmd->SetParameterName("R",false);
  SetRadiusMagnetosphereCmd->SetRange("R > 0.");
  SetRadiusMagnetosphereCmd->SetUnitCategory("Length");
  SetRadiusMagnetosphereCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  cmd_name = "/PLANETOCOS/BFIELD/BfieldVsPositions";
  ComputeBfieldAtDifferentPositions = new G4UIcommand(cmd_name,this);
  guidance = "Compute the magnetic field at # positions";
  ComputeBfieldAtDifferentPositions->SetGuidance(guidance);
  ComputeBfieldAtDifferentPositions->SetParameter(coord_sys_param);
  ComputeBfieldAtDifferentPositions->SetParameter(altitude_param);
  ComputeBfieldAtDifferentPositions->SetParameter(dalt_param);
  ComputeBfieldAtDifferentPositions->SetParameter(nAlt_param);
  ComputeBfieldAtDifferentPositions->SetParameter(length_unit_param);
  ComputeBfieldAtDifferentPositions->SetParameter(latitude_param );
  ComputeBfieldAtDifferentPositions->SetParameter(dlat_param);
  ComputeBfieldAtDifferentPositions->SetParameter(nLat_param);
  ComputeBfieldAtDifferentPositions->SetParameter(longitude_param);
  ComputeBfieldAtDifferentPositions->SetParameter(dlong_param);
  ComputeBfieldAtDifferentPositions->SetParameter(nLong_param);
  ComputeBfieldAtDifferentPositions->SetParameter(angle_unit_param);
  ComputeBfieldAtDifferentPositions->SetParameter(output_file_param);
  ComputeBfieldAtDifferentPositions->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/BfieldVsXYZ";
  ComputeBfieldAtDifferentPositions1 = new G4UIcommand(cmd_name,this);
  guidance = "Compute the magnetic field at # positions";
  ComputeBfieldAtDifferentPositions1->SetGuidance(guidance);
  ComputeBfieldAtDifferentPositions1->SetParameter(altitude_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(dalt_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(nAlt_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(x_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(dx_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(nx_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(y_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(dy_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(ny_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(length_unit_param);
  ComputeBfieldAtDifferentPositions1->SetParameter(output_file_param);
  ComputeBfieldAtDifferentPositions1->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/DRAW/TraceBlineAtPositions";
  TraceBlineAtDifferentPositions = new G4UIcommand(cmd_name,this);
  guidance = "Trace blines at # positions";
  TraceBlineAtDifferentPositions->SetGuidance(guidance);
  TraceBlineAtDifferentPositions->SetParameter(coord_sys_param);
  TraceBlineAtDifferentPositions->SetParameter(altitude_param);
  TraceBlineAtDifferentPositions->SetParameter(length_unit_param);
  TraceBlineAtDifferentPositions->SetParameter(latitude_param );
  TraceBlineAtDifferentPositions->SetParameter(dlat_param);
  TraceBlineAtDifferentPositions->SetParameter(nLat_param);
  TraceBlineAtDifferentPositions->SetParameter(longitude_param);
  TraceBlineAtDifferentPositions->SetParameter(dlong_param);
  TraceBlineAtDifferentPositions->SetParameter(nLong_param);
  TraceBlineAtDifferentPositions->SetParameter(angle_unit_param);
  TraceBlineAtDifferentPositions->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/DRAW/TraceBlineAtXYAltPos";
  TraceBlineAtDifferentPositions1 = new G4UIcommand(cmd_name,this);
  guidance = "Trace blines at # Positions";
  TraceBlineAtDifferentPositions1->SetGuidance(guidance);
  TraceBlineAtDifferentPositions1->SetParameter(altitude_param);
  TraceBlineAtDifferentPositions1->SetParameter(x_param );
  TraceBlineAtDifferentPositions1->SetParameter(dx_param);
  TraceBlineAtDifferentPositions1->SetParameter(nx_param);
  TraceBlineAtDifferentPositions1->SetParameter(y_param);
  TraceBlineAtDifferentPositions1->SetParameter(dy_param);
  TraceBlineAtDifferentPositions1->SetParameter(ny_param);  
  TraceBlineAtDifferentPositions1->SetParameter(length_unit_param);
  TraceBlineAtDifferentPositions1->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SetHomogeneousField";
  SetHomogeneousFieldCmd = new G4UIcmdWith3VectorAndUnit(cmd_name,this); 
  guidance = "Set the homogeneous magnetic field";
  SetHomogeneousFieldCmd->SetGuidance(guidance);
  SetHomogeneousFieldCmd->SetParameterName("Bx","By","Bz",false);
  SetHomogeneousFieldCmd->SetUnitCategory("Magnetic flux density");
  SetHomogeneousFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  cmd_name = "/PLANETOCOS/BFIELD/SetHomogeneousFieldFromSelectedModel";
  SetHomogeneousFieldFromSelectedModelCmd = new G4UIcommand(cmd_name, this);
  guidance = "Set the homogeneous magnetic field to the field value at ";
  guidance = "the given position";
  SetHomogeneousFieldFromSelectedModelCmd->SetGuidance(guidance);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(altitude_param);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(length_unit_param);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(latitude_param);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(longitude_param);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(coord_sys_param);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(int_model_parameter);
  SetHomogeneousFieldFromSelectedModelCmd->SetParameter(ext_model_parameter);
  
  SetHomogeneousFieldFromSelectedModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  
  cmd_name = "/PLANETOCOS/BFIELD/SetWorldCenterPositionForFlatGeometry";
  SetWorldCenterPositionCmd = new G4UIcommand(cmd_name, this);
  guidance = "Set the geographical position of the center of the world box in ";
  guidance += "the case of the Flat geometry";
  SetWorldCenterPositionCmd->SetParameter(latitude_param);
  SetWorldCenterPositionCmd->SetParameter(longitude_param);
  SetWorldCenterPositionCmd->SetParameter(coord_sys_param);
  SetWorldCenterPositionCmd->AvailableForStates(G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SwitchOn";
  SwitchOnFieldCmd= new G4UIcmdWithoutParameter(cmd_name,this);
  SwitchOnFieldCmd->SetGuidance("Activate the effect of the magnetic field");
  SwitchOnFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/SwitchOff";
  SwitchOffFieldCmd= new G4UIcmdWithoutParameter(cmd_name,this);
  SwitchOffFieldCmd->SetGuidance("Desactivate the effect of the magnetic field");
  SwitchOffFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/ComputeInterpolMatrixForFlatGeometry";
  ComputeInterpolMatrixForFlatGeometryCmd = new G4UIcommand(cmd_name, this);
  guidance = "Compute the interpolation matrix for the flat grid model";
  ComputeInterpolMatrixForFlatGeometryCmd->SetGuidance(guidance);
  ComputeInterpolMatrixForFlatGeometryCmd->SetParameter(nx_param);
  ComputeInterpolMatrixForFlatGeometryCmd->SetParameter(ny_param);
  ComputeInterpolMatrixForFlatGeometryCmd->SetParameter(nz_param);
  ComputeInterpolMatrixForFlatGeometryCmd->SetParameter(output_file_param);
  ComputeInterpolMatrixForFlatGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/BFIELD/ReadInterpolMatrixForFlatGeometry";
  ReadInterpolMatrixForFlatGeometryCmd = new G4UIcmdWithAString(cmd_name, this);
  guidance = "Read the interpolation matrix for the flat grid model";
  ReadInterpolMatrixForFlatGeometryCmd->SetGuidance(guidance);
  ReadInterpolMatrixForFlatGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}
////////////////////////////////////////////////////////////////////////////////
//
PlanetMagneticFieldMessenger::~PlanetMagneticFieldMessenger()
{ delete SetEpsilonCmd;
  delete SetG4DeltaChordCmd; 
  delete SetDeltaIntersectionCmd;
  delete ResetIntegrationParametersCmd;
  delete SetStepperCmd;
  delete SetTimeOfBCmd;
  delete SetStartDateCmd;
  delete Setnmax_GaussCmd;
  delete SetInternalFieldCmd;
  delete SetExternalFieldCmd;
  delete SetMagnetopauseModelCmd;
  delete SetDipoleB0Cmd;
  delete SetDipoleShiftCmd; 
  delete SetDipoleAxisCmd; 
  delete SetConsiderDipoleShiftCmd;
  delete ComputeBfieldAtDifferentPositions;
  delete ComputeBfieldAtDifferentPositions1;
  delete TraceBlineAtDifferentPositions;
  delete TraceBlineAtDifferentPositions1;
  delete SetHomogeneousFieldCmd;
  delete SetHomogeneousFieldFromSelectedModelCmd;
  delete SetWorldCenterPositionCmd;
  delete SwitchOnFieldCmd; 
  delete SwitchOffFieldCmd;
  delete ComputeInterpolMatrixForFlatGeometryCmd;
  delete ReadInterpolMatrixForFlatGeometryCmd;
  delete mainDir;
  delete IntegrationDir;
  delete MagnetoDir;
}		  
////////////////////////////////////////////////////////////////////////////////
//
bool PlanetMagneticFieldMessenger::SetMotherNewValue(G4UIcommand * command,G4String newValues)
{  
  
   // integration parameters command 
  
  bool result=true; 
   
  if (command == SetEpsilonCmd) 
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
            theMotherField->SetTimeOfB(SetTimeOfBCmd->GetNewDoubleValue(newValues)/s);
   
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
	    	    
  else if ( command == SetMagnetopauseModelCmd )
            theMotherField->SetMagnetopauseModel(newValues);
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
   else if( command == SetRadiusMagnetosphereCmd ) 
            theMotherField->SetRadiusMagnetosphere(SetRadiusMagnetosphereCmd->GetNewDoubleValue(newValues));
       
  
  
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
	G4String coor_sys, length_unit, int_model, ext_model;
	std::istringstream is((char*)paramString);
	is >> altitude >> length_unit >> latitude >> longitude >>coor_sys
		>>int_model>>ext_model;
	altitude *= G4UnitDefinition::GetValueOf(length_unit);
	latitude *=degree;
	longitude *=degree;
	theMotherField->ComputeOrientationAndBfieldConst(altitude,
				  latitude, longitude,
	                          coor_sys,
			           int_model, ext_model);
	
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
  
  else if (command == ComputeInterpolMatrixForFlatGeometryCmd){
  	const char* paramString=newValues;
	G4int nx,ny,nz;
	G4String output_file_name;
	std::istringstream is((char*)paramString);
	is >> nx >> ny >> nz >>output_file_name;

	theMotherField->ComputeInterpolMatrixForFlatGeometry(nx,ny,nz,output_file_name);	
  }
  else if (command == ReadInterpolMatrixForFlatGeometryCmd){
  	theMotherField->ReadInterpolMatrixForFlatGeometry(newValues);
  }
 
  else result = false;
  
  return result;

}












