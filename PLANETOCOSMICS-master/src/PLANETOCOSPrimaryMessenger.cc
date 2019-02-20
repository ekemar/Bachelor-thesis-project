#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "PLANETOCOSPrimaryMessenger.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4VSensitiveDetector.hh" 
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "SpaceCoordinatePlanet.hh"
// #include "G4VViewer.hh"

PLANETOCOSPrimaryMessenger::PLANETOCOSPrimaryMessenger(PLANETOCOSPrimaryGeneratorAction
* aGenerator)
:myGenerator(aGenerator)
{ G4UIparameter* param;
  
  
  G4UIparameter* coorsys_param_with_PLAG = new G4UIparameter("CoordSys",'s',false);
  G4UIparameter* coorsys_param_without_PLAG = new G4UIparameter("CoordSys",'s',false);
  
  std::vector<G4String> theListOfCoorSystems = 
          SpaceCoordinatePlanet::GetInstance()->GetListOfCoordinateSystems();
  G4String coorsys_candidates =" ";
  G4String coorsys_candidates1 =" ";
  for (unsigned int i=0;i<theListOfCoorSystems.size();i++){
  	coorsys_candidates+= theListOfCoorSystems[i]+" ";
  	if ( theListOfCoorSystems[i] != "PLAG") {
		coorsys_candidates1+= theListOfCoorSystems[i]+" ";
	}
  }	  
  coorsys_param_with_PLAG->SetParameterCandidates(coorsys_candidates);
  coorsys_param_without_PLAG->SetParameterCandidates(coorsys_candidates1);	  
  
  
  
  myStartGeoDir = new G4UIdirectory("/PLANETOCOS/SOURCE/");
  myStartGeoDir
   ->SetGuidance("Definition of the start position and direction of a particle in space coordinate system");
   
 /* myRigidityFilterDir = new G4UIdirectory("/PLANETOCOS/MAGNETIC/");
  myRigidityFilterDir
   ->SetGuidance("Command to compute the rigidity filter");
  */ 

  SetPositionVectorCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/SetPositionVector",this);
  SetPositionVectorCmd
    ->SetGuidance("Set the vector start position in your selected sapce system");
  SetPositionVectorCmd
    ->SetGuidance("[usage] /PLANETOCOS/SOURCE/SetPositionVector  CoordSys X Y Z unit" ); 
  SetPositionVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPositionVectorCmd->SetParameter(coorsys_param_without_PLAG);
  param = new G4UIparameter("X",'d',false);
  SetPositionVectorCmd->SetParameter(param);
  param = new G4UIparameter("Y",'d',false);
  SetPositionVectorCmd->SetParameter(param);
  param = new G4UIparameter("Z",'d',false);
  SetPositionVectorCmd->SetParameter(param);
  param = new G4UIparameter("Unit",'s',false);
  SetPositionVectorCmd->SetParameter(param);
  
  
  SetPositionCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/SetPosition",this);
  SetPositionCmd
    ->SetGuidance("Set the altitude, latitude and longitude for the start position in your selected sapce system");
  SetPositionCmd
    ->SetGuidance("[usage] /PLANETOCOS/SOURCE/SetPosition  CoordSys altitude length_ unit latitude longitude angle_unit" ); 
  SetPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPositionCmd->SetParameter(coorsys_param_with_PLAG);
  param = new G4UIparameter("altitude",'d',false);
  SetPositionCmd->SetParameter(param);
  param = new G4UIparameter("length_unit",'s',false);
  SetPositionCmd->SetParameter(param);
  param = new G4UIparameter("latitude",'d',false);
  SetPositionCmd->SetParameter(param);
  param = new G4UIparameter("longitude",'d',false);
  SetPositionCmd->SetParameter(param);
  param = new G4UIparameter("angle_unit",'s',false);
  SetPositionCmd->SetParameter(param);
  
  
  
  /*SetPositionOnDipoleMagneticShellCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/SetPositionOnDipoleMagneticShell",this);
  SetPositionOnDipoleMagneticShellCmd
    ->SetGuidance("Set the L parameter in Re, latitude and longitude for the start PositionOnDipoleMagneticShell in your selected sapce system");
  SetPositionOnDipoleMagneticShellCmd
    ->SetGuidance("[usage] /PLANETOCOS/SOURCE/SetPositionOnDipoleMagneticShell CoordSys L latitude longitude angle_unit" ); 
  SetPositionOnDipoleMagneticShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  param = new G4UIparameter("CoordSys",'s',false);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(param);
  param = new G4UIparameter("L",'d',false);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(param);
  param = new G4UIparameter("latitude",'d',false);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(param);
  param = new G4UIparameter("longitude",'d',false);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(param);
  param = new G4UIparameter("angle_unit",'s',false);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(param);*/
  
  
  
  SetDirectionVectorCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/SetDirectionVector",this);
  SetDirectionVectorCmd
    ->SetGuidance("Set the vector start Direction in your selected sapce system");
  SetDirectionVectorCmd
    ->SetGuidance("[usage] /PLANETOCOS/SOURCE/SetDirectionVector  CoordSys X Y Z unit" ); 
  SetDirectionVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDirectionVectorCmd->SetParameter(coorsys_param_without_PLAG);
  param = new G4UIparameter("X",'d',false);
  SetDirectionVectorCmd->SetParameter(param);
  param = new G4UIparameter("Y",'d',false);
  SetDirectionVectorCmd->SetParameter(param);
  param = new G4UIparameter("Z",'d',false);
  SetDirectionVectorCmd->SetParameter(param);
 
  
  
  SetDirectionCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/SetDirection",this);
  SetDirectionCmd
    ->SetGuidance("Set the zenith and azimuth for the start Direction in your selected sapce system");
  SetDirectionCmd
    ->SetGuidance("[usage] /PLANETOCOS/SOURCE/SetDirection  CoordSys zenith azimuth angle_unit" ); 
  SetDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDirectionCmd->SetParameter(coorsys_param_without_PLAG);
  param = new G4UIparameter("zenith",'d',false);
  SetDirectionCmd->SetParameter(param);
  param = new G4UIparameter("azimuth",'d',false);
  SetDirectionCmd->SetParameter(param);
  param = new G4UIparameter("angle_unit",'s',false);
  SetDirectionCmd->SetParameter(param);
  
  G4String title_cmd, guidance;
  G4UIparameter* enuc_min = new G4UIparameter("enuc_min",'d',false);
  enuc_min->SetParameterRange("enuc_min>0");
  
  G4UIparameter* enuc_max = new G4UIparameter("enuc_max",'d',false);
  enuc_max->SetParameterRange("enuc_max>0");
  
  G4UIparameter* enuc_unit = new G4UIparameter("enuc_unit",'s',true);
  enuc_unit->SetDefaultValue("GeV/nuc");
  enuc_unit->SetParameterCandidates("GeV/nuc MeV/nuc TeV/nuc PeV/nuc keV/nuc");
  
  G4UIparameter* particle = new G4UIparameter("particle",'s',false);
  
  G4UIparameter* phi_mod = new G4UIparameter("phi_mod",'d',false);
  phi_mod->SetParameterRange("phi_mod >= 0 && phi_mod < 1500");
  
  G4UIparameter* altitude = new G4UIparameter("altitude",'d',false);
  altitude->SetParameterRange(" altitude >= 0");
  
  G4UIparameter* altitude_unit = new G4UIparameter("altitude_unit",'s',true);
  altitude_unit->SetDefaultValue("km");
  altitude_unit->SetParameterCandidates("km m");
  
  // position and direction command for flat geometry
  //--------------------------------------------------
  title_cmd ="/PLANETOCOS/SOURCE/SetPositionForFlatGeometry";
  SetPositionVectorFlatCaseCmd = new G4UIcmdWith3VectorAndUnit(title_cmd,this);
  guidance ="Set position in case of flat geometry. The third parameter represents the altitude.";
  SetPositionVectorFlatCaseCmd->SetGuidance(guidance);
  SetPositionVectorFlatCaseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetPositionVectorFlatCaseCmd->SetUnitCandidates("Rplanet rplanet km m");
  
  title_cmd ="/PLANETOCOS/SOURCE/SetDirectionVectorForFlatGeometry";
  SetDirectionVectorFlatCaseCmd = new G4UIcmdWith3Vector(title_cmd,this);
  guidance ="Set direction in case of flat geometry.";
  SetDirectionVectorFlatCaseCmd->SetGuidance(guidance);
  SetDirectionVectorFlatCaseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  title_cmd ="/PLANETOCOS/SOURCE/SetDirectionForFlatGeometry";
  SetDirectionFlatCaseCmd = new G4UIcommand(title_cmd,this);
  guidance ="Set the zenith and azimuth angle of the incoming direction. ";
  SetDirectionFlatCaseCmd->SetGuidance(guidance);
  SetDirectionFlatCaseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  param = new G4UIparameter("zenith",'d',false);
  SetDirectionFlatCaseCmd->SetParameter(param);
  param = new G4UIparameter("azimuth",'d',false);
  SetDirectionFlatCaseCmd->SetParameter(param);
  param = new G4UIparameter("angle_unit",'s',false);
  SetDirectionFlatCaseCmd->SetParameter(param);
 
  
  title_cmd ="/PLANETOCOS/SOURCE/SelectSolMaxGalacticFlux";
  SelectSolMaxGalacticFluxCmd =new G4UIcommand(title_cmd,this);
  guidance ="Select the flux of galactic comsic ray at solar maximum ";
  SelectSolMaxGalacticFluxCmd->SetGuidance(guidance);
  SelectSolMaxGalacticFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SelectSolMaxGalacticFluxCmd->SetParameter(particle);
  SelectSolMaxGalacticFluxCmd->SetParameter(enuc_min);
  SelectSolMaxGalacticFluxCmd->SetParameter(enuc_max);
  SelectSolMaxGalacticFluxCmd->SetParameter(enuc_unit);
 /* SelectSolMaxGalacticFluxCmd->SetParameter(altitude);
  SelectSolMaxGalacticFluxCmd->SetParameter(altitude_unit);*/
  
  
  title_cmd ="/PLANETOCOS/SOURCE/SelectSolMinGalacticFlux";
  SelectSolMinGalacticFluxCmd =new G4UIcommand(title_cmd,this);
  guidance ="Select the flux of galactic comsic ray at solar minimum ";
  SelectSolMinGalacticFluxCmd->SetGuidance(guidance);
  SelectSolMinGalacticFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SelectSolMinGalacticFluxCmd->SetParameter(particle);
  SelectSolMinGalacticFluxCmd->SetParameter(enuc_min);
  SelectSolMinGalacticFluxCmd->SetParameter(enuc_max);
  SelectSolMinGalacticFluxCmd->SetParameter(enuc_unit);
 /* SelectSolMinGalacticFluxCmd->SetParameter(altitude);
  SelectSolMinGalacticFluxCmd->SetParameter(altitude_unit);*/
  
  
  
  title_cmd ="/PLANETOCOS/SOURCE/SelectMeanGalacticFlux";
  SelectMeanGalacticFluxCmd =new G4UIcommand(title_cmd,this);
  guidance ="Select the mean flux of galactic comsic ray";
  SelectMeanGalacticFluxCmd->SetGuidance(guidance);
  SelectMeanGalacticFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SelectMeanGalacticFluxCmd->SetParameter(particle);
  SelectMeanGalacticFluxCmd->SetParameter(enuc_min);
  SelectMeanGalacticFluxCmd->SetParameter(enuc_max);
  SelectMeanGalacticFluxCmd->SetParameter(enuc_unit);
 /* SelectMeanGalacticFluxCmd->SetParameter(altitude);
  SelectMeanGalacticFluxCmd->SetParameter(altitude_unit);*/
  
  
  
  
  title_cmd ="/PLANETOCOS/SOURCE/SelectModulatedGalacticFlux";
  SelectModulatedGalacticFluxCmd =new G4UIcommand(title_cmd,this);
  guidance ="Select a modulated flux of galactic comsic ray";
  SelectModulatedGalacticFluxCmd->SetGuidance(guidance);
  SelectModulatedGalacticFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SelectModulatedGalacticFluxCmd->SetParameter(particle);
  SelectModulatedGalacticFluxCmd->SetParameter(phi_mod);
  SelectModulatedGalacticFluxCmd->SetParameter(enuc_min);
  SelectModulatedGalacticFluxCmd->SetParameter(enuc_max);
  SelectModulatedGalacticFluxCmd->SetParameter(enuc_unit);
 /* SelectModulatedGalacticFluxCmd->SetParameter(altitude);
  SelectModulatedGalacticFluxCmd->SetParameter(altitude_unit);*/
  
  
  
  title_cmd ="/PLANETOCOS/SOURCE/ReadPrimaryFlux";
  ReadPrimaryFluxCmd = new G4UIcmdWithAString(title_cmd,this);
  guidance ="Read a file describing the primary spectrum";
  ReadPrimaryFluxCmd->SetGuidance(guidance);
  ReadPrimaryFluxCmd->SetParameterName("file_name",false);
  ReadPrimaryFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  title_cmd ="/PLANETOCOS/SOURCE/RandomIsotropicDistribution";				//PVD
  RandomIsotropicDistributionCmd = new G4UIcommand(title_cmd,this);				//PVD
  RandomIsotropicDistributionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);		//PVD

  param = new G4UIparameter("particle_distribution",'s',false);
  RandomIsotropicDistributionCmd->SetParameter(param);
  RandomIsotropicDistributionCmd->SetParameter(coorsys_param_with_PLAG);
  param = new G4UIparameter("length_unit",'s',false);
  RandomIsotropicDistributionCmd->SetParameter(param);
  param = new G4UIparameter("angle_unit",'s',false);
  RandomIsotropicDistributionCmd->SetParameter(param);
  param = new G4UIparameter("planet_radius",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);
  param = new G4UIparameter("atmo_height",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);
  param = new G4UIparameter("det_height",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);  
  param = new G4UIparameter("low_theta",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);  
  param = new G4UIparameter("up_theta",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param); 
  param = new G4UIparameter("low_phi",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);  
  param = new G4UIparameter("up_phi",'d',false);
  RandomIsotropicDistributionCmd->SetParameter(param);        
  param = new G4UIparameter("eventNumber",'i',false);
  RandomIsotropicDistributionCmd->SetParameter(param);  
  
  
    //commands
  title_cmd="/PLANETOCOS/SOURCE/ConsiderCutoff";
  SetConsiderCutoffCmd = new G4UIcmdWithABool(title_cmd,this);
  guidance ="If true the cut off is taken into account";
  SetConsiderCutoffCmd->SetGuidance(guidance); 
  SetConsiderCutoffCmd->SetParameterName("Consider cutoff",true);
  SetConsiderCutoffCmd->SetDefaultValue(true);
  SetConsiderCutoffCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  
  title_cmd ="/PLANETOCOS/SOURCE/SetCutoffRigidityValue";
  SetCutoffRigidityCmd = new G4UIcmdWithADoubleAndUnit(title_cmd,this);
  guidance ="Set the value of the cutoff rigidity ";
  SetCutoffRigidityCmd->SetGuidance(guidance);
  SetCutoffRigidityCmd->SetParameterName("CutoffRigidity",false);
  SetCutoffRigidityCmd->SetUnitCategory("Electric potential");
  SetCutoffRigidityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  title_cmd ="/PLANETOCOS/SOURCE/ReadCutoffVsDirection";
  ReadCutoffVsDirectionCmd = new G4UIcmdWithAString(title_cmd,this);

 /* SetZenithCmd =
     new G4UIcmdWithADouble("/PLANETOCOS/STARTPOSITION/SetZenith",this);
  SetZenithCmd->SetGuidance("Define the zenith for the incoming direction");
  SetZenithCmd->SetParameterName("Zenith",false);
  SetZenithCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetAzimuthCmd =
    new G4UIcmdWithADouble("/PLANETOCOS/STARTPOSITION/SetAzimuth",this);
  SetAzimuthCmd->SetGuidance("Define the azimuth for the incoming direction");
  SetAzimuthCmd->SetParameterName("Azimuth",false);
  SetAzimuthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
  SetRigidityCmd =
    new G4UIcmdWithADoubleAndUnit("/PLANETOCOS/SOURCE/SetRigidity",this);
  SetRigidityCmd->SetGuidance("Define the Rigidity of the incoming particle");
  SetRigidityCmd->SetParameterName("Rigidity",false);
  SetRigidityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  /*SetPlaPositionCmd =
    new G4UIcmdWith3Vector("/PLANETOCOS/STARTPOSITION/SetPlaPosition",this);
  SetPlaPositionCmd->SetGuidance("Define the position  for the particle in PLA");
  SetPlaPositionCmd->SetParameterName("Altitude km","Latitude °","Longitude °",
                                                       false);
  SetPlaPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ComputeStartPositionCmd= 
   new G4UIcmdWithoutParameter("/PLANETOCOS/STARTPOSITION/ComputeStartPosition",this);
       
  ComputeStartPositionCmd->SetGuidance("Compute the start position");
   ComputeStartPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
   
  AddValuesToRigidityVectorCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/AddValuesToRigidityVector",this);
  AddValuesToRigidityVectorCmd
    ->SetGuidance("Increase the rigidity vector with new values");
  AddValuesToRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  param = new G4UIparameter("val0",'d',false);
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  param = new G4UIparameter("dval",'d',false);
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  param = new G4UIparameter("nvalues",'i',true);
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  
  SetDefaultRigidityVectorCmd= 
   new G4UIcmdWithoutParameter("/PLANETOCOS/MAGNETIC/SetDefaultRigidityVector",this);
  SetDefaultRigidityVectorCmd->SetGuidance("Set default rigidity vector");
  SetDefaultRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  ResetRigidityVectorCmd= 
   new G4UIcmdWithoutParameter("/PLANETOCOS/MAGNETIC/ResetRigidityVector",this);
  ResetRigidityVectorCmd->SetGuidance("Reset rigidity vector");
  ResetRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  G4String cmd_name;
  
  cmd_name = "/PLANETOCOS/SOURCE/verbose";
  SetVerbosityCmd = new G4UIcmdWithAnInteger(cmd_name,this);
  SetVerbosityCmd->SetGuidance("If >0 primary information are printed");
  SetVerbosityCmd->SetParameterName("verbosity", false);
  SetVerbosityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/SOURCE/BfieldAtPrimaryPosition";
  PrintBfieldAtPrimaryCmd = new G4UIcmdWithoutParameter(cmd_name,this);
  PrintBfieldAtPrimaryCmd->SetGuidance("Print the Bfield at primary position");
  PrintBfieldAtPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


#ifdef USE_ROOT_SOURCE
  
  
  
  title_cmd="/PLANETOCOS/SOURCE/UseRootSource";
  SetUseRootSourceCmd = new G4UIcmdWithABool(title_cmd,this);
  guidance ="If true the root source is used";
  SetUseRootSourceCmd->SetGuidance(guidance); 
  SetUseRootSourceCmd->SetParameterName("Use Root Source",true);
  SetUseRootSourceCmd->SetDefaultValue(true);
  SetUseRootSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  title_cmd ="/PLANETOCOS/SOURCE/SetVar1NameForRootSource";
  SetVar1NameForRootSourceCmd = new G4UIcmdWithAString(title_cmd,this);
  guidance ="Set the type of the var1 of the Root Source";
  SetVar1NameForRootSourceCmd->SetGuidance(guidance);
  SetVar1NameForRootSourceCmd->SetParameterName("file_name",false);
  SetVar1NameForRootSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  title_cmd ="/PLANETOCOS/SOURCE/SetVar2NameForRootSource";
  SetVar2NameForRootSourceCmd = new G4UIcmdWithAString(title_cmd,this);
  guidance ="Set the type of the var2 of the Root Source";
  SetVar2NameForRootSourceCmd->SetGuidance(guidance);
  SetVar2NameForRootSourceCmd->SetParameterName("file_name",false);
  SetVar2NameForRootSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  ReadFirstHistoForRootSourceCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/ReadFirstHistoForRootSource",this);
  ReadFirstHistoForRootSourceCmd
    ->SetGuidance("Read the first histogram of the Root Source");
  ReadFirstHistoForRootSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  ReadFirstHistoForRootSourceCmd->SetParameter(new G4UIparameter("FileName",'s',false));
  ReadFirstHistoForRootSourceCmd->SetParameter(new G4UIparameter("Path",'s',false));
  
  ReadSecondHistoForRootSourceCmd = new
             G4UIcommand("/PLANETOCOS/SOURCE/ReadSecondHistoForRootSource",this);
  ReadSecondHistoForRootSourceCmd
    ->SetGuidance("Read the second histogram of the Root Source");
  ReadSecondHistoForRootSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  ReadSecondHistoForRootSourceCmd->SetParameter(new G4UIparameter("FileName",'s',false));
  ReadSecondHistoForRootSourceCmd->SetParameter(new G4UIparameter("Path",'s',false));

  

    
#endif  
   
 
  
  
}

PLANETOCOSPrimaryMessenger::~PLANETOCOSPrimaryMessenger()
{ delete   myStartGeoDir;

}

void PLANETOCOSPrimaryMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{/* if( command == SetZenithCmd)
           myGenerator->SetZenith(SetZenithCmd->GetNewDoubleValue(newValues));
  
  else if( command == SetAzimuthCmd)
           myGenerator->SetAzimuth(SetAzimuthCmd->GetNewDoubleValue(newValues));
  */
  
   
   
   
    if ( command == SetPositionVectorCmd )
       { const char* paramString;
         paramString=newValues;
         G4double x,y,z;
	 G4String  coord_sys,unit_str;
	 std::istringstream is((char*)paramString);
         is >> coord_sys >> x >> y >> z >> unit_str ;
         G4ThreeVector PLAposition =
	     G4ThreeVector(x,y,z)*G4UnitDefinition::GetValueOf(unit_str);
	 G4bool test;
	 test= myGenerator->SetPosition(coord_sys,PLAposition);
       }
       
    if ( command == SetPositionCmd )
       { const char* paramString;
         paramString=newValues;
         G4double altitude,latitude,longitude;
	 G4String  coord_sys,length_unit,angle_unit;
	 std::istringstream is((char*)paramString);
         is >> coord_sys >> altitude >> length_unit >> latitude >> 
	       longitude >> angle_unit ;
         altitude *=G4UnitDefinition::GetValueOf(length_unit);
	 longitude *=G4UnitDefinition::GetValueOf(angle_unit);
	 latitude *=G4UnitDefinition::GetValueOf(angle_unit);
	 G4bool test;
	 test= myGenerator->SetPosition(coord_sys,altitude,longitude,latitude);
       }  
  /*   if ( command == SetPositionOnDipoleMagneticShellCmd )
       { const char* paramString;
         paramString=newValues;
         G4double L,latitude,longitude;
	 G4String  coord_sys,angle_unit;
	 std::istringstream is((char*)paramString);
         is >> coord_sys >> L >> latitude >> 
	       longitude >> angle_unit ;
	 longitude *=G4UnitDefinition::GetValueOf(angle_unit);
	 latitude *=G4UnitDefinition::GetValueOf(angle_unit);
	 G4bool test;
	 test=
myGenerator->SetPositionOnDipoleMagneticShell(coord_sys,L,latitude,longitude);
       }  */   
        
    if ( command == SetDirectionVectorCmd )
       { const char* paramString;
         paramString=newValues;
         G4double x,y,z;
	 G4String  coord_sys;
	 std::istringstream is((char*)paramString);
         is >> coord_sys >> x >> y >> z ;
         G4ThreeVector PLADirection = G4ThreeVector(x,y,z);
	 PLADirection/=PLADirection.mag();
	 G4bool test;
	 test= myGenerator->SetDirection(coord_sys,PLADirection);
       }
       
     if ( command == SetDirectionCmd )
       { const char* paramString;
         paramString=newValues;
         G4double azimuth,zenith;
	 G4String  coord_sys,angle_unit;
	 std::istringstream is((char*)paramString);
         is >> coord_sys >> zenith >> azimuth >> angle_unit ;
	 zenith *=G4UnitDefinition::GetValueOf(angle_unit);
	 azimuth *=G4UnitDefinition::GetValueOf(angle_unit);
	 G4bool test;
	 test= myGenerator->SetDirection(coord_sys,zenith,azimuth);
       }     
  //flat geoemtry position and direction command
   else if( command == SetPositionVectorFlatCaseCmd){
   	const PLANETOCOSGeometryConstruction* theGeometry =
			dynamic_cast<const PLANETOCOSGeometryConstruction*>
					(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
   	if (theGeometry->GetGeometryType() != "SPHERICAL"){
		myGenerator->SetPosition(
			SetPositionVectorFlatCaseCmd->GetNew3VectorValue(newValues));
	} 
	else {
		G4cout<<"This command should be used for flat geometry"<<std::endl;
	}
   }
   
   else if( command == SetDirectionVectorFlatCaseCmd){
        const PLANETOCOSGeometryConstruction* theGeometry =
			dynamic_cast<const PLANETOCOSGeometryConstruction*>
					(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
   	
   	if (theGeometry->GetGeometryType() != "SPHERICAL"){
		
		myGenerator->SetDirection(
			SetDirectionVectorFlatCaseCmd->GetNew3VectorValue(newValues));
	} 
	else {
		G4cout<<"This command should be used for flat geometry"<<std::endl;
	}
   }
   
   else if( command == SetDirectionFlatCaseCmd){
   	const PLANETOCOSGeometryConstruction* theGeometry =
			dynamic_cast<const PLANETOCOSGeometryConstruction*>
					(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
   	
   	if (theGeometry->GetGeometryType() != "SPHERICAL"){
		const char* paramString;
         	paramString=newValues;
         	G4double azimuth,zenith;
	 	G4String  angle_unit;
	 	std::istringstream is((char*)paramString);
         	is >> zenith >> azimuth >> angle_unit ;
	 	zenith *=G4UnitDefinition::GetValueOf(angle_unit);
	 	azimuth *=G4UnitDefinition::GetValueOf(angle_unit);
		myGenerator->SetDirection(zenith,azimuth);
		
	} 
	else {
		G4cout<<"This command should be used for flat geometry"<<std::endl;
	}
   }
   
   else if( command == SetRigidityCmd)
           myGenerator->SetRigidity(SetRigidityCmd->GetNewDoubleValue(newValues));
  
  
  
   else if ( command == AddValuesToRigidityVectorCmd )
       { const char* paramString;
         paramString=newValues;
         G4double val0,dval;
	 G4int   nvalues;
         std::istringstream is((char*)paramString);
         is >> val0 >> dval >> nvalues;
         myGenerator->AddValuesToRigidityVector(nvalues,val0,dval);
       }
   
   else if( command == SetDefaultRigidityVectorCmd) 
                       myGenerator->SetDefaultRigidityVector();    
       
   else if( command == ResetRigidityVectorCmd) 
                       myGenerator->ResetRigidityVector();  
		       
   else if( command == SetVerbosityCmd)
        myGenerator->SetVerbosity(SetVerbosityCmd->GetNewIntValue(newValues)); 
   else if( command == PrintBfieldAtPrimaryCmd)
   	                                 myGenerator->PrintBfieldAtPrimary();      
   
   else if (command == SelectSolMaxGalacticFluxCmd)
    {const char* paramString=newValues;
     std::istringstream is((char*)paramString);
     G4String particle, enuc_unit;
     G4double enuc_min,enuc_max;
     is >> particle>>enuc_min>>enuc_max>>enuc_unit;
     enuc_min*=G4UnitDefinition::GetValueOf(enuc_unit);
     enuc_max*=G4UnitDefinition::GetValueOf(enuc_unit);
     
     G4ParticleDefinition* aParticle = 
              G4ParticleTable::GetParticleTable()->FindParticle(particle);
     myGenerator->SelectModulatedGalacticFluxAtSolMax(aParticle,
                                                      enuc_min,
						      enuc_max); 
    }
  else if (command == SelectSolMinGalacticFluxCmd)
    {const char* paramString=newValues;
     std::istringstream is((char*)paramString);
     G4String particle, enuc_unit;
     G4double enuc_min,enuc_max;
     is >> particle>>enuc_min>>enuc_max>>enuc_unit;
     enuc_min*=G4UnitDefinition::GetValueOf(enuc_unit);
     enuc_max*=G4UnitDefinition::GetValueOf(enuc_unit);
     G4ParticleDefinition* aParticle = 
              G4ParticleTable::GetParticleTable()->FindParticle(particle);
     myGenerator->SelectModulatedGalacticFluxAtSolMin(aParticle,
                                                      enuc_min,
						      enuc_max); 
    }
  else if (command == SelectMeanGalacticFluxCmd)
    {const char* paramString=newValues;
     std::istringstream is((char*)paramString);
     G4String particle, enuc_unit;
     G4double enuc_min,enuc_max;
     is >> particle>>enuc_min>>enuc_max>>enuc_unit;
     enuc_min*=G4UnitDefinition::GetValueOf(enuc_unit);
     enuc_max*=G4UnitDefinition::GetValueOf(enuc_unit);
     G4ParticleDefinition* aParticle = 
              G4ParticleTable::GetParticleTable()->FindParticle(particle);
     myGenerator->SelectMeanModulatedGalacticFlux(aParticle,
                                                  enuc_min,
						  enuc_max); 
    } 
   else if (command == SelectModulatedGalacticFluxCmd)
    {const char* paramString=newValues;
     std::istringstream is((char*)paramString);
     G4String particle, enuc_unit;
     G4double enuc_min,enuc_max,phi_mod;
     is >> particle>>phi_mod>>enuc_min>>enuc_max>>enuc_unit;
     enuc_min*=G4UnitDefinition::GetValueOf(enuc_unit);
     enuc_max*=G4UnitDefinition::GetValueOf(enuc_unit);
     G4ParticleDefinition* aParticle = 
              G4ParticleTable::GetParticleTable()->FindParticle(particle);
     myGenerator->SelectModulatedGalacticFlux(aParticle, enuc_min,
                                                  enuc_max,
						  phi_mod); 
    }
    else if (command ==  ReadPrimaryFluxCmd)
     myGenerator->ReadUserSpectrum(newValues);
     
    else if (command ==  RandomIsotropicDistributionCmd)//PVD
     {const char* paramString;
      paramString=newValues;
      
      G4String particle_distribution, coord_sys, length_unit, angle_unit;	
      G4double planet_radius, atmo_height, det_height, low_theta, up_theta, low_phi, up_phi;
      G4int eventNumber;
	
      std::istringstream is((char*)paramString);
      is >> particle_distribution >> coord_sys >> length_unit >> angle_unit >> planet_radius >> atmo_height >> det_height >> low_theta >> up_theta >> low_phi >> up_phi >> eventNumber;
            
      myGenerator->RandomIsotropicDistribution(particle_distribution, coord_sys, length_unit, angle_unit, planet_radius, atmo_height, det_height, low_theta,  up_theta,  low_phi, up_phi, eventNumber);
     }     
    
    else if (command == SetConsiderCutoffCmd)
      myGenerator->SetConsiderCutoff(SetConsiderCutoffCmd->GetNewBoolValue(newValues));
    
    else if (command == SetCutoffRigidityCmd)
      myGenerator->SetCutoffRigidity(SetCutoffRigidityCmd->GetNewDoubleValue(newValues));
  
    else if (command == ReadCutoffVsDirectionCmd)   
      myGenerator->ReadCutOffVsDirection(newValues);

#ifdef USE_ROOT_SOURCE
    else if (command == SetUseRootSourceCmd)
     				myGenerator->SetUseRootSource(SetUseRootSourceCmd->GetNewBoolValue(newValues));
    
    else if (command == SetVar1NameForRootSourceCmd)
     				myGenerator->SetVar1Type(newValues);
    else if (command == SetVar2NameForRootSourceCmd)
     				myGenerator->SetVar2Type(newValues);				
    
    else if (command ==  ReadFirstHistoForRootSourceCmd) {
    	const char* paramString=newValues;
     	std::istringstream is((char*)paramString);
     	G4String file_name, path;
    	is >> file_name>>path;
    	myGenerator->ReadROOTFirstSourceHisto(file_name,path); 
    }
    else if (command ==  ReadSecondHistoForRootSourceCmd) {
    	const char* paramString=newValues;
     	std::istringstream is((char*)paramString);
     	G4String file_name, path;
    	is >> file_name>>path;
    	myGenerator->ReadROOTSecondSourceHisto(file_name,path); 
    }
    

   
    
    
    
     	
#endif      

}
