#include"MagneticShieldingToolMessenger.hh"
#include"MagneticShieldingTool.hh"
#include"MagneticShieldingToolEventAction.hh"

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
#include "SpaceCoordinatePlanet.hh"

//units
#include "G4UnitsTable.hh"

MagneticShieldingToolMessenger::MagneticShieldingToolMessenger
            (MagneticShieldingTool* aMagneticShieldingTool )
{ theMagneticShieldingTool = aMagneticShieldingTool;
 
  //directory
  //----------
  
  MagneticShieldingToolDir=new G4UIdirectory("/PLANETOCOS/MAGNETIC/");
  MagneticShieldingToolDir->SetGuidance("Interactive commands to compute cutoff rigidities, visualiuse trajectories and magnetic field lines");
 
  //parameters
  //----------

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
  
  
  
  G4UIparameter* length_unit_param = new G4UIparameter("Length Unit",'s',false);
  length_unit_param->SetParameterCandidates("Re km m");
  
  G4UIparameter* angle_unit_param = new G4UIparameter("Angle Unit",'s',false);
  angle_unit_param->SetParameterCandidates("degree deg rad radian mrad milliradian");
  
  G4UIparameter* altitude_param = new G4UIparameter("Altitude",'d',false);
  
  G4UIparameter* longitude_param = new G4UIparameter("Longitude",'d',false);
  G4UIparameter* dlong_param = new G4UIparameter("delta longitude",'d',false);
  G4UIparameter* nLong_param = new G4UIparameter("number of longitudes",'i',false);
  
  G4UIparameter* latitude_param = new G4UIparameter("Latitude",'d',false);
  G4UIparameter* dlat_param = new G4UIparameter("delta latitude",'d',false);
  G4UIparameter* nLat_param = new G4UIparameter("number of latitudes",'i',false);
  
  G4UIparameter* zenith_param = new G4UIparameter("Zenith",'d',false);
  G4UIparameter* dZenith_param = new G4UIparameter("delta zenith",'d',false);
  G4UIparameter* nZenith_param = new G4UIparameter("number of zenith angles",'i',false);
  
  G4UIparameter* azimuth_param = new G4UIparameter("Azimuth",'d',false);
  G4UIparameter* dAzimuth_param = new G4UIparameter("delta azimuth",'d',false);
  G4UIparameter* nAzimuth_param = new G4UIparameter("number of azimuth angles",'i',false);
  
  G4UIparameter* L_param = new G4UIparameter("L",'d',false);
  L_param->SetParameterRange(" L > 0.");
  L_param->SetGuidance("L shell parameter in Re");
  
  G4UIparameter* output_file_param = new G4UIparameter("OutputFile name",'s',false);
  //G4UIparameter* input_file_param = new G4UIparameter("InputFile name",'s',false);
  
  G4UIparameter* start_time_param = new G4UIparameter("start time",'d',false);
  G4UIparameter* dTime_param = new G4UIparameter("delta time",'d',false);
  G4UIparameter* nTime_param = new G4UIparameter("number of different times",'i',false);
  G4UIparameter*  time_unit = new G4UIparameter("Time Unit",'s',false);
  time_unit->SetParameterCandidates("s second hour  minute day");
  
  // Scenario commands
  //-----------------

  SetStopAltitudeCmd = new G4UIcmdWithADoubleAndUnit("/PLANETOCOS/MAGNETIC/SetStopAltitude",this);
  SetStopAltitudeCmd->SetGuidance("Set the altitude below which a trajectory is stopped");
  SetStopAltitudeCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  //forward trajectory
  BlineCmd = new G4UIcmdWithoutParameter("/PLANETOCOS/DRAW/TraceBline",this);
  BlineCmd->SetGuidance("Trace a Magnetic field line");
  BlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  
  //forward trajectory
  ParticleTrajectoryCmd = new G4UIcmdWithAnInteger("/PLANETOCOS/MAGNETIC/ComputeTrajectories",this);
  ParticleTrajectoryCmd->SetParameterName("nb_particles",true);
  ParticleTrajectoryCmd->SetDefaultValue(1);
  ParticleTrajectoryCmd->SetGuidance("Compute  particle trajectories in the magnetic field.");
  ParticleTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  //bacward trajectory
  ReverseParticleTrajectoryCmd = new G4UIcmdWithAnInteger("/PLANETOCOS/MAGNETIC/ComputeReverseTrajectories",this);
  ReverseParticleTrajectoryCmd->SetParameterName("nb_particles",true);
  ReverseParticleTrajectoryCmd->SetDefaultValue(1);
  G4String guidance = "Compute particle trajectories bacward in time in the magnetic field.";
  ReverseParticleTrajectoryCmd->SetGuidance(guidance);
  ReverseParticleTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  //asymptotic directions
  ComputeRigidityFilterCmd =
    new G4UIcmdWithAString("/PLANETOCOS/MAGNETIC/ComputeAsymptoticDirections",this);
  guidance ="Compute the asymptotic directions and cutoff rigidities ";
  guidance +="for the primary position and direction that you should have ";
  guidance +="previously defined by using the /PLANETOCOS/SOURCE commands";
  ComputeRigidityFilterCmd->SetGuidance(guidance);
  ComputeRigidityFilterCmd->SetParameterName("Output_file name",false);
  ComputeRigidityFilterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  /*// asymptotic directions
  //----------------------
  
  ComputeDirectionFilterCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/AsymptoticDirVsDirection",this);
  ComputeDirectionFilterCmd
    ->SetGuidance("Computing of asymptotic direction and filter value for different zenith and azimuth");
  ComputeDirectionFilterCmd
    ->SetGuidance(
   "[usage] /PLANETOCOS/MAGNETIC/AsymptoticDirVsDirection Rigidity unit CoordSys CosZen0 dCosZen nZen  Azim0 dAzim nAzim  OutputFile" ); 
  
  ComputeDirectionFilterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   param = new G4UIparameter("Rigidity",'d',false);
  param->SetDefaultValue("5.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("unit",'s',false);
  param->SetDefaultValue("GV");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  
  
  param = new G4UIparameter("CoordSys",'s',false);
  param->SetDefaultValue("GEOID");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  
  
  
  param = new G4UIparameter("CosZen0",'d',false);
  param->SetDefaultValue("1.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
 
  param = new G4UIparameter("dCosZen",'d',false);
  param->SetDefaultValue("-.1");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("nZen",'i',false);
  param->SetDefaultValue("10.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("Azim0",'d',false);
  param->SetDefaultValue("0");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("dAzim",'d',false);
  param->SetDefaultValue("40");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("nAzim",'i',false);
  param->SetDefaultValue("9");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("OutputFile",'s',true);
  param->SetDefaultValue("Filter.txt");
  ComputeDirectionFilterCmd->SetParameter(param);*/
  
  
  
  // rigidity cutoff for different position
  
  RCutoffVsPositionCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsPosition",this);
  RCutoffVsPositionCmd
    ->SetGuidance("Computing cutoff rigidities  at different longitudes and latitudes");
  RCutoffVsPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsPositionCmd->SetParameter(coorsys_param_with_PLAG);
  RCutoffVsPositionCmd->SetParameter(altitude_param);
  RCutoffVsPositionCmd->SetParameter(length_unit_param);
  RCutoffVsPositionCmd->SetParameter(latitude_param);
  RCutoffVsPositionCmd->SetParameter(dlat_param);
  RCutoffVsPositionCmd->SetParameter(nLat_param);
  RCutoffVsPositionCmd->SetParameter(longitude_param);
  RCutoffVsPositionCmd->SetParameter(dlong_param);
  RCutoffVsPositionCmd->SetParameter(nLong_param);
  RCutoffVsPositionCmd->SetParameter(zenith_param);
  RCutoffVsPositionCmd->SetParameter(azimuth_param);
  RCutoffVsPositionCmd->SetParameter(angle_unit_param);
  RCutoffVsPositionCmd->SetParameter(output_file_param);
  
  // rigidity cutoff for different position on dipole-magnetic shell
  
 /* RCutoffVsPositionOnLShellCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsPositionOnLShell",this);
  guidance ="Computing of cutoff rigidities at different longitudes and latitudes ";	     
  guidance +="on a given dipole magnetic shell";	     
  RCutoffVsPositionOnLShellCmd->SetGuidance(guidance);
  RCutoffVsPositionOnLShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsPositionOnLShellCmd->SetParameter(mag_shell_coord_sys_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(L_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(latitude_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(dlat_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(nLat_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(longitude_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(dlong_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(nLong_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(zenith_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(azimuth_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(angle_unit_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(output_file_param);
  
  // rigidity cutoff for spenvis trajectory
  
  RCutoffVsSpenvisTrajectoryCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsSpenvisTrajectory",this);
  RCutoffVsSpenvisTrajectoryCmd
    ->SetGuidance("Computing cutoff rigidities along a trajectory");
  RCutoffVsSpenvisTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(input_file_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(output_file_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(zenith_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(azimuth_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(angle_unit_param);
  
  // rigidity cutoff for spenvis position grid
  
  RCutoffVsSpenvisPositionGridCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsSpenvisPositionGrid",this);
  RCutoffVsSpenvisPositionGridCmd
    ->SetGuidance("Computing cutoff rigidities along a trajectory");
  RCutoffVsSpenvisPositionGridCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(input_file_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(output_file_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(zenith_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(azimuth_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(angle_unit_param);
 

  */
  // rigidity cutoff for different directions
 
  RCutoffVsDirectionCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsDirection",this);
  guidance ="Computing of cutoff rigidities for different directions ";
  guidance +="defined by a grid of different zenith and azimuth angles";
  RCutoffVsDirectionCmd->SetGuidance(guidance);
  RCutoffVsDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsDirectionCmd->SetParameter(coorsys_param_with_PLAG);
  RCutoffVsDirectionCmd->SetParameter(zenith_param);
  RCutoffVsDirectionCmd->SetParameter(dZenith_param);
  RCutoffVsDirectionCmd->SetParameter(nZenith_param);
  RCutoffVsDirectionCmd->SetParameter(azimuth_param);
  RCutoffVsDirectionCmd->SetParameter(dAzimuth_param);
  RCutoffVsDirectionCmd->SetParameter(nAzimuth_param);
  RCutoffVsDirectionCmd->SetParameter(output_file_param);
  

  
 // Cutoff for different times
 
  RCutoffVsTimeCmd = new
             G4UIcommand("/PLANETOCOS/MAGNETIC/RCutoffVsTime",this);
  RCutoffVsTimeCmd->SetGuidance("Computing of cutoff rigidities  for different times");
  RCutoffVsTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsTimeCmd->SetParameter(start_time_param);
  RCutoffVsTimeCmd->SetParameter(dTime_param);
  RCutoffVsTimeCmd->SetParameter(nTime_param);
  RCutoffVsTimeCmd->SetParameter(time_unit);
  RCutoffVsTimeCmd->SetParameter(output_file_param);
  
  
  //AutoDetectionOfPenumbra for the Earth only 
  SetAutoDetectionOfPenumbra = new 
      G4UIcmdWithABool("/PLANETOCOS/MAGNETIC/AutomaticDetectionOfPenumbra",this);
  guidance =" If the parameter Auto is true the penumbra is detected ";
  guidance += " automatically for the Earth (only) for the case where cutoff rigidities ";
  guidance += "are computed in function of position, direction of ";
  guidance += "incidence and time ";
  SetAutoDetectionOfPenumbra->SetGuidance(guidance);  
  SetAutoDetectionOfPenumbra->SetParameterName("Auto",true);
  SetAutoDetectionOfPenumbra->SetDefaultValue(true);    
  SetAutoDetectionOfPenumbra->AvailableForStates(G4State_PreInit,G4State_Idle);
 
 
  //Command for registering results in SpenvisCSV file
  
 /* SetRegisterResultsInSpenvisCSVFileCmd = new 
      G4UIcmdWithABool("/PLANETOCOS/MAGNETIC/RegisterResultsInSpenvisCSVFile",this);
  guidance =" If true (fase)l the results of cutoff rigidity ";
  guidance += "and asymptotic direction computation are (not) registered in a SpenvisCSV file";
  SetRegisterResultsInSpenvisCSVFileCmd->SetGuidance(guidance);  
  SetRegisterResultsInSpenvisCSVFileCmd->SetParameterName("Register",false);
  SetRegisterResultsInSpenvisCSVFileCmd->SetDefaultValue(true);    
  SetRegisterResultsInSpenvisCSVFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  SetSpenvisCSVFileNameCmd =
    new G4UIcmdWithAString("/PLANETOCOS/MAGNETIC/SetSpenvisCSVFileName",this);
  guidance ="Set the the name of the Spenvis csv file";
  SetSpenvisCSVFileNameCmd->SetGuidance(guidance);
  SetSpenvisCSVFileNameCmd->SetParameterName("SpenvisCSV file  name",false);
  SetSpenvisCSVFileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
 */ 
 
  InitialiseCmd = new G4UIcmdWithoutParameter("/PLANETOCOS/InitialiseForMagneticShieldingStudy",this);
  InitialiseCmd->SetGuidance("Initialise for magnetic shielding study. Atmosphere and soil are not considered. And  Hadronic and EM physics are neglected.");
  InitialiseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

MagneticShieldingToolMessenger::~MagneticShieldingToolMessenger()
{
}		  

void MagneticShieldingToolMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ if (command == SetStopAltitudeCmd) {
  		G4double stop_alt = SetStopAltitudeCmd->GetNewDoubleValue(newValues);
                theMagneticShieldingTool->SetStopAltitude(stop_alt);
  }
  else if (command == BlineCmd) {
                theMagneticShieldingTool->Bline();
  }
  
  else if (command == ParticleTrajectoryCmd) {
  		G4int nb_part = ParticleTrajectoryCmd->GetNewIntValue(newValues);
                theMagneticShieldingTool->ParticleTrajectoryInBfield(nb_part);
  }	       
  else if (command == ReverseParticleTrajectoryCmd){
  		G4int nb_part = ReverseParticleTrajectoryCmd->GetNewIntValue(newValues);
                theMagneticShieldingTool->ReverseParticleTrajectoryInBfield(nb_part);
  }	       	       
  else if ( command == ComputeRigidityFilterCmd) 
                  theMagneticShieldingTool->ComputeRigidityFilter(newValues);
/*  else if ( command == ComputeDirectionFilterCmd ){
  	const char* paramString=newValues;
        G4double  cos_zen0,dcos_zen,azim0,dazim,rigidity;
	G4int   nzen,nazim;
	G4String OutputFile,CoordSys,unit;
        std::istringstream is((char*)paramString);
        is >> rigidity >> unit >> CoordSys>>cos_zen0 >> dcos_zen >> nzen >>
	azim0 >> dazim >> nazim>>OutputFile;
	rigidity *=G4UnitDefinition::GetValueOf(unit);
        theMagneticShieldingTool->ComputeDirectionFilter
	          		(CoordSys,rigidity,cos_zen0,dcos_zen,nzen,
		   		 azim0*degree,dazim*degree,nazim,OutputFile);
  }
  */		  			      
  else if ( command == RCutoffVsPositionCmd ){
  	const char* paramString=newValues;
        G4double  Alt,lat0,dlat,long0,dlong,zenith,azimuth;
	G4int   nlat,nlong;
	G4String OutputFile,CoordSys,LengthUnit,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>> Alt >> LengthUnit>>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	zenith >> azimuth >> AngleUnit >> OutputFile;
	Alt *=G4UnitDefinition::GetValueOf(LengthUnit);
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMagneticShieldingTool->RCutoffVsPosition
	          	(CoordSys,Alt,lat0,dlat,nlat,long0,dlong,nlong,zenith,
		   	 azimuth,OutputFile);
  }
  /*
  else if ( command == RCutoffVsPositionOnLShellCmd ){
  	const char* paramString=newValues;
       	G4double  L,lat0,dlat,long0,dlong,zenith,azimuth;
	G4int   nlat,nlong;
	G4String OutputFile,CoordSys,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>> L >>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	zenith >> azimuth >> AngleUnit >> OutputFile;
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMagneticShieldingTool->RCutoffVsPositionOnDipoleMagneticShell
	                 (CoordSys,L,lat0,dlat,nlat,long0,dlong,nlong,zenith,
		         azimuth,OutputFile);
  }
  
  else if ( command == RCutoffVsSpenvisTrajectoryCmd ){
  	const char* paramString=newValues;
        G4double  zenith,azimuth;
	G4String OutputFile,InputFile,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> InputFile>> OutputFile >>zenith >> azimuth >> AngleUnit;
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMagneticShieldingTool->RCutoffVsSpenvisTrajectory
	          	(InputFile, OutputFile, zenith, azimuth);
  }
  else if ( command == RCutoffVsSpenvisPositionGridCmd ){
  	const char* paramString=newValues;
        G4double  zenith,azimuth;
	G4String OutputFile,InputFile,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> InputFile>> OutputFile >>zenith >> azimuth >> AngleUnit;
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theMagneticShieldingTool->RCutoffVsSpenvisPositionGrid
	          	(InputFile, OutputFile, zenith, azimuth);
  }
 */ 
  else if ( command == RCutoffVsDirectionCmd ){
  	const char* paramString=newValues;
        G4double  zen0,dzen,azim0,dazim;
	G4int   nzen,nazim;
	G4String OutputFile,CoordSys;
        std::istringstream is((char*)paramString);
        is >> CoordSys>>zen0 >> dzen >> nzen>>azim0 >> dazim >> nazim>>OutputFile;
        theMagneticShieldingTool->RCutoffVsDirection
	          (CoordSys,zen0*degree,dzen*degree,nzen,azim0*degree,dazim*degree,nazim,OutputFile);
  }		  
  
  else if ( command == RCutoffVsTimeCmd ){
  	const char* paramString=newValues;
       	G4double  time0,dtime;
	G4int   ntime;
	G4String TimeUnit,OutputFile;
        std::istringstream is((char*)paramString);
        is >> time0 >> dtime >> ntime;
	is >>TimeUnit;
	is>>OutputFile;
	G4double t_unit = G4UnitDefinition::GetValueOf(TimeUnit);
        time0 *=t_unit;
	dtime *=t_unit;
        theMagneticShieldingTool->RCutoffVsTime(time0,dtime,ntime,OutputFile);
  }
  else if (command == SetAutoDetectionOfPenumbra){
  	G4bool aBool = SetAutoDetectionOfPenumbra->GetNewBoolValue(newValues);
	theMagneticShieldingTool->SetAutomaticDetectionPenumbra(aBool);
  }
  /*
  else if (command == SetRegisterResultsInSpenvisCSVFileCmd){
  	G4bool aBool = SetRegisterResultsInSpenvisCSVFileCmd->GetNewBoolValue(newValues);
	theMagneticShieldingTool->SetRegisterResultsInSpenvisFile(aBool);
  }
  else if ( command == SetSpenvisCSVFileNameCmd) 
                  theMagneticShieldingTool->SetSpenvisFileName(newValues);
	*/
  else if (command == InitialiseCmd) {
                theMagneticShieldingTool->InitialiseForMagneticShieldingStudy();
  }    				      
}












