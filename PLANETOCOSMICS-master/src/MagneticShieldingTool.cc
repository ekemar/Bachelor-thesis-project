#include"MagneticShieldingTool.hh"
#include"MagneticShieldingToolMessenger.hh"
#include"MagneticShieldingToolAnalysisManager.hh"
#include"MagneticShieldingToolEventAction.hh"
#include"MagneticShieldingToolSteppingAction.hh"
#include"PLANETOCOSSpenvisManager.hh"
#include"G4UImanager.hh"
#include"G4ParticleDefinition.hh"
#include"G4Event.hh"
#include"G4UnitsTable.hh"
#include"G4ios.hh"
#include"G4RunManager.hh"
#include"G4LogicalVolumeStore.hh"
#include"PLANETOCOSTrackingAction.hh"
#include"PLANETOCOSSteppingAction.hh"
#include"PLANETOCOSStackingAction.hh"
#include"PLANETOCOSApplicationScenario.hh"
#include"PLANETOCOSEventAction.hh"
#include"PLANETOCOSPhysicsList.hh"
#include"PLANETOCOSGeometryConstruction.hh"
#include"PLANETOCOSPrimaryGeneratorAction.hh"
#include"DurationManager.hh"
#include"PlanetManager.hh"
#include"PlanetSoil.hh"
#include"PlanetEquationOfMotion.hh"
#include"PlanetMagneticField.hh"
#include"SpaceCoordinatePlanet.hh"
#include"PLANETOCOSAnalysisManager.hh"
#include"G4ProcessTable.hh"
#include"G4ParticleDefinition.hh"
#include"MYGeneralParticleSource.hh"
#include"G4ChargedGeantino.hh"

#include <vector>
 
MagneticShieldingTool* MagneticShieldingTool::instance = 0;

//////////////////////////////////////////////////////////////////
MagneticShieldingTool::MagneticShieldingTool()
{
 fMessenger = new MagneticShieldingToolMessenger(this);
 fSteppingAction = new MagneticShieldingToolSteppingAction(this) ;
 fEventAction = new  MagneticShieldingToolEventAction(this);
 fSteppingAction->SetBlineMode(false);
 
 theSpenvisManager = PLANETOCOSSpenvisManager::GetInstance();
  
 theBfield = PlanetManager::GetInstance()->GetMagneticField(); 
 TotalDurationReached = false;
 AutomaticDetectionOfThePenumbra = false;
 stop_altitude = 0.;
 
 if (PlanetManager::GetInstance()->GetPlanetName() == "Earth"){
 	AutomaticDetectionOfThePenumbra = true;
	stop_altitude = 19.99*km;
 }
 
 RegisterAsymptoticDirection =false;
 
}
///////////////////////////////////////////////////////////////////////
MagneticShieldingTool::~MagneticShieldingTool()
{delete fMessenger;
 delete fSteppingAction;
 delete fEventAction; 
}
////////////////////////////////////////////////////////////////////////////////
//
MagneticShieldingTool* MagneticShieldingTool::GetInstance()
{ if (instance == 0 ) instance = new MagneticShieldingTool;
  return instance;
}		  
////////////////////////////////////////////////////////////////////

void MagneticShieldingTool::BeginOfRunAction(const G4Run*)
{ 
  theGeometry = 
   	const_cast<PLANETOCOSGeometryConstruction*>
		(dynamic_cast<const PLANETOCOSGeometryConstruction*>(
				G4RunManager::GetRunManager()
					->GetUserDetectorConstruction()));

  geometry_type=theGeometry->GetGeometryType();

  ZpositionAtZeroAltitude = theGeometry->GetZpositionAtZeroAltitude(); 
  MagnetosphereMaxStepLength = theGeometry->GetMagnetosphereMaxStepLength();
  AtmosphereMaxStepLength = theGeometry->GetAtmosphereMaxStepLength();
  RplanetAtEquator=PlanetManager::GetInstance()->GetRplanetAtEquator();

  theGeometry->SetAtmosphereMaxStepLength(MagnetosphereMaxStepLength);
 

  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
  
  DurationManager* theDurationManager = DurationManager::GetInstance();
  TotalDurationReached = !theDurationManager->CheckDurationAtBeginOfRun();
  if (TotalDurationReached) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  }
  

 
}  
///////////////////////////////////////////////////////////////////////
void MagneticShieldingTool::EndOfRunAction(const G4Run* )
{ ResetMaxStepLength();
  
}
//////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::Bline()
{ 
   // Switch on the field
  theBfield->SwitchOn();
  
  //change user actions
  //fSteppingAction->SetBlineMode(true);
  G4double old_stop_altitude = stop_altitude;
  stop_altitude =0.;
  SetActions();
  
  
  //PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(false);

  //get a pointer to the particle source 
  MYGeneralParticleSource* theParticleSource =
                 fUserPrimaryAction->GetParticleSource();

  G4ParticleDefinition* theOldParticleDefinition
                   = theParticleSource->GetParticleDefinition();
  theParticleSource->SetParticleDefinition(G4ChargedGeantino::ChargedGeantino());

 //get a pointer to the equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();

  //select the bline equation
  theEquation->SetEquationType("BFIELD_LINE");

  //SwitchOff all physics process

  SwitchOffAllProcessesExceptTransportation();

  //southward integration
  
  theEquation->SetReverseTimeMode(true);
  G4RunManager::GetRunManager()->BeamOn(1);
  
 //northward integration
  theEquation->SetReverseTimeMode(false);
  G4RunManager::GetRunManager()->BeamOn(1); 
  
  //Back to Lorentz Equation
  theEquation->SetEquationType("LORENTZ_MOTION");
  
  //Back to old particle definition
   
  if (theOldParticleDefinition)  
      theParticleSource->SetParticleDefinition(theOldParticleDefinition);
   //go back to user actions
  ResetUserActions(); 
  fSteppingAction->SetBlineMode(false);
  stop_altitude = old_stop_altitude; 
  
  //SwitchOn all physics process
  SwitchOnAllProcesses();  
}

////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ParticleTrajectoryInBfield(G4int nparticles)
{ 
   // Switch on the field
  theBfield->SwitchOn();
  
  //change user actions
  
 // SetActions();
  
  //PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(false); 
  
  //get a pointer to the equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();
 
  //select the Lorentz equation
  theEquation->SetEquationType("LORENTZ_MOTION");
 
 	
  //SwitchOff all physics process
  
  SwitchOffAllProcessesExceptTransportation();
	
	
  //integration
  theEquation->SetReverseTimeMode(false);
  G4RunManager::GetRunManager()->BeamOn(nparticles);
  	
  //SwitchOn all physics process
  SwitchOnAllProcesses();
  
  //PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(true);
  
  //go back to user actions
 // ResetUserActions(); 
}
///////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ReverseParticleTrajectoryInBfield(G4int nparticles)
{ // Switch on the field
  theBfield->SwitchOn();
  
  //change user actions
  
   SetActions();
  //PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(false); 
  
  //get a pointer to the equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();
 
  //select the Lorentz equation
  theEquation->SetEquationType("LORENTZ_MOTION");
 
 	
  //SwitchOff all physics process
  
  SwitchOffAllProcessesExceptTransportation();
	
	
  //integration
  theEquation->SetReverseTimeMode(true);
  G4RunManager::GetRunManager()->BeamOn(nparticles);
  theEquation->SetReverseTimeMode(false);
  	
  //SwitchOn all physics process
  SwitchOnAllProcesses();
  
  //PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(true); 
  ResetUserActions();
   
}
////////////////////////////////////////////////////////////////
bool MagneticShieldingTool::CheckIfParticleShouldBeKilledAndSetMaxStepLength(G4ThreeVector pos)
{
 //G4cout<<geometry_type<<std::endl;
 if (geometry_type == "SPHERICAL"){
	
	
	G4double r=pos.mag();
	if (theBfield->OutsideMagnetosphere(pos)) return true;
	if (r> (RplanetAtEquator +stop_altitude)){
		G4double max_step_length= std::min(r-(RplanetAtEquator +stop_altitude),MagnetosphereMaxStepLength);
		max_step_length = std::max(max_step_length,1.*km);
		theGeometry->SetAtmosphereMaxStepLength(max_step_length);
		theGeometry->SetMagnetosphereMaxStepLength(max_step_length);
		return false;
	}	
	
	
	
	
	G4double altitude,latitude,longitude;
	altitude =10.;
	
	SpaceCoordinatePlanet::GetInstance()
				->ComputePLAGCoordinatesFromPLAPosition
	              			(pos,altitude,longitude,latitude);
    	
	if (altitude < stop_altitude ) {
		
		/*G4cout<<"STOP"<<std::endl;
		G4cout<<"Rplanet "<<RplanetAtEquator/km<<std::endl;
		G4cout<<stop_altitude/km<<std::endl;
		G4cout<<altitude/km<<std::endl;*/
		return true;
	}	
	
	//change step_length
	if (PlanetManager::GetInstance()->GetPlanetName() != "Jupiter"){
		if (r< (RplanetAtEquator +stop_altitude)){
			theGeometry->SetAtmosphereMaxStepLength(1.*km);
			theGeometry->SetMagnetosphereMaxStepLength(1.*km);
			return false;
		}
		else if (r< (RplanetAtEquator +stop_altitude+10.*km)){
			theGeometry->SetAtmosphereMaxStepLength(2.*km);
			theGeometry->SetMagnetosphereMaxStepLength(2.*km);
			return false;
		}
		else if (r< (RplanetAtEquator +stop_altitude+100.*km)){
		
			theGeometry->SetAtmosphereMaxStepLength(10.*km);
			theGeometry->SetMagnetosphereMaxStepLength(10.*km);
			return false;
		}
	}
	if (PlanetManager::GetInstance()->GetPlanetName() == "Jupiter"){
		G4double max_step_length= std::max(1.2*(altitude-stop_altitude),.5*km);
		max_step_length=std::min(max_step_length,MagnetosphereMaxStepLength);
		theGeometry->SetAtmosphereMaxStepLength(max_step_length);
		theGeometry->SetMagnetosphereMaxStepLength(max_step_length);
		
		return false;
	}	
	
		
	
  }
  else {
  	G4double altitude = pos.z()-ZpositionAtZeroAltitude;
	if (altitude< stop_altitude) return true;
	else if (altitude < stop_altitude+10.*km){
		theGeometry->SetAtmosphereMaxStepLength(2.*km);
		theGeometry->SetMagnetosphereMaxStepLength(2.*km);
		return false;
	}
	else if (altitude < stop_altitude+100.*km){
		theGeometry->SetAtmosphereMaxStepLength(10.*km);
		theGeometry->SetMagnetosphereMaxStepLength(10.*km);
		return false;
		
	}
	else {  
		G4double max_step_length= std::min(altitude-stop_altitude,MagnetosphereMaxStepLength);
		if (max_step_length <5.*m) max_step_length = 100.*m;
		theGeometry->SetAtmosphereMaxStepLength(max_step_length);
		theGeometry->SetMagnetosphereMaxStepLength(max_step_length);
		return false;
	}
	
  }
  return false;
 
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::SetAutomaticDetectionPenumbra(G4bool aBool)
{ if (PlanetManager::GetInstance()->GetPlanetName() == "Earth"){
 	AutomaticDetectionOfThePenumbra = aBool; 
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ComputeRigidityFilter(G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  RegisterAsymptoticDirection=true;
  
  
  //change user actions
 SetActions();
  
  //csv file
  if (RegisterResultsInSpenvisFile)
  	theSpenvisManager->InitialiseAsymptoticDirectionCSVFile(); 
 
  MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  theAnalysisManager->OpenAsymptoticDirectionFile(outputfile_name);
  ComputeRigidityCutoff();

  theAnalysisManager->CloseAsymptoticDirectionFile(Rc,Rm,Rs); 
			  
  if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->SaveCSVFile(SpenvisFileName);
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
  RegisterAsymptoticDirection = false; 
  ResetUserActions();
 
}
///////////////////////////////////////////////////////////////////////////////
//
/*void MagneticShieldingTool::ComputeDirectionFilter
                          (G4String CoordSys,G4double rigidity,G4double cos_zen0, 
			   G4double delta_cos_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim,
			   G4String outputfile_name)
{ 
  if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  //Open the ascii file
  MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  theAnalysisManager->OpenAsymptoticDirVsDirFile(CoordSys,rigidity,outputfile_name);
  
 
  
  //get a pointer to the equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();
 
 
 
  //select the bline equation
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);
 
 
  //pointer on PLANETOCOSTrackingAction
  PLANETOCOSTrackingAction* theTrackingAction =
    (PLANETOCOSTrackingAction*)
          G4RunManager::GetRunManager()->GetUserTrackingAction();
  theTrackingAction->SetRegisterLastPoint(true);
 
  //get a pointer to the primary generator action
  PLANETOCOSPrimaryGeneratorAction* fUserPrimaryAction =
   (PLANETOCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  fUserPrimaryAction->SetRigidity(rigidity);
  
  for (int i=0; i<nzen;i++){
  	G4double zenith=std::acos( cos_zen0+double(i)*delta_cos_zen);
   	for (int j=0; j<nazim;j++){
      		G4double azimuth=azim0+double(j)*delta_azim;
       		if (fUserPrimaryAction->SetDirection(CoordSys,zenith,azimuth)){ 
          		G4RunManager::GetRunManager()->BeamOn(1);
          		theAnalysisManager
	     			->RegisterAsymptoticDirVsDir
	                            (zenith,azimuth,LastPositions[0],LastMomentaOnCharge[0],
				                    FilterValues[0]);
						    
	  		LastPositions.clear();
	  		LastMomentaOnCharge.clear();
	  		FilterValues.clear();   
	  	}
		else { 
	 		j=nazim;
	  		i=nzen;
	 	}          
       
    	}
  }
  theEquation->SetReverseTimeMode(false);
  theTrackingAction->SetRegisterLastPoint(false); 
  theAnalysisManager->CloseAsciiFile(); 
}
*/ 
///////////////////////////////////////////////////////////////////////////////
// 
void MagneticShieldingTool::RCutoffVsPosition
                          (G4String CoordSys, G4double Altitude,
                           G4double lat0,G4double delta_lat, G4int nlat,
                           G4double long0, G4double delta_long, G4int nlong,
			   G4double zenith,G4double azimuth,
			   G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  //Open the ascii file
  MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  theAnalysisManager->OpenCutoffVsPositionFile(outputfile_name,CoordSys,
	                                       Altitude, zenith, azimuth);
  //test cvs
  if (RegisterResultsInSpenvisFile)
  	theSpenvisManager->InitialiseCutoffVsPosCSVFile(CoordSys,zenith, azimuth);
 
  //change user actions
  SetActions();	  
  
  //avoid drawing particles
  G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
  fUserEventAction->SetDrawTrajectory(false);	  
  
   
	
  for (int i=0; i<nlat;i++){
  	G4double latitude=lat0+double(i)*delta_lat;
   	for (int j=0; j<nlong;j++){
      		G4double longitude=long0+double(j)*delta_long;
       		G4double test = 
          		fUserPrimaryAction
	    			->SetPositionAndDirection(CoordSys,Altitude,longitude,latitude, 
	                              zenith,azimuth);
		G4ThreeVector PLAPosition = 
				fUserPrimaryAction->GetPLAPosition();
		
		std::vector<double> ILat;	
		ILat.push_back(0);
		ILat.push_back(0);
		ILat.push_back(0);	
/*		ILat.push_back(SpaceCoordinatePlanet::GetInstance()
					->ComputePLADipoleInvariantLatitude("PLA",PLAPosition ));		      
		if (getenv("INVARIANT_LATITUDE")){
			std::vector<double> LMc = theBfield->ComputeMcIlwainLParameter(PLAPosition,89.999*degree);
			for (int i=0;i<2;i++){
				
				if (LMc[i] >= 1. ) ILat.push_back(std::acos(std::sqrt(1/LMc[i])));
				else ILat.push_back(-999999.*degree);
				
			}	
		}	
	
*/					      
			
				      
       		if (test){
			//G4cout<<"MAG1"<<std::endl; 
         		if (AutomaticDetectionOfThePenumbra ) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
				//G4cout<<"MAG2"<<std::endl;
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude,longitude,ILat);
				//G4cout<<"MAG3"<<std::endl;
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(theBfield->GetReferenceDate().SpenvisModifiedJulianDay());
					values.push_back(Altitude/km);
					values.push_back(latitude/degree);
					values.push_back(longitude/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
       		else{
         		j=nlong;
	  		i=nlat;
	 	}
       }
  }
  if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);
  theAnalysisManager->CloseAsciiFile();   

  fUserEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;
  
  ResetUserActions();		    
}


///////////////////////////////////////////////////////////////////////////////////////////
// 
void MagneticShieldingTool::RCutoffVsSpenvisPositionGrid(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    			     G4double zenith , 
				    			     G4double azimuth)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  
  //Read the position grid
  
  std::vector<double> altitude;
  std::vector<double>  latitude;
  std::vector<double>  longitude;
  
  
  G4bool test = theSpenvisManager
  			->ReadSpenvisPositionGridCSVFile(SpenvisCSVFileName,
							 altitude,
							 latitude,
							 longitude,
    							 1);
 
  //change user actions
  SetActions();	
  if (test) {
  	
	
	
	//Open the ascii file
  	MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  	theAnalysisManager->OpenCutoffVsSpenvisPositionGridFile(outputfile_name,
								"PLAID",
	                                        		zenith, 
								azimuth);
 	//test cvs
  	if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->InitialiseCutoffVsPosCSVFile("PLAID",zenith, azimuth);
 
   
  	//avoid drawing particles
	//----------------------- 
	G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
  	fUserEventAction->SetDrawTrajectory(false);	  
  
  	for (unsigned int i=0; i<altitude.size();i++){
	 	G4double test1 = fUserPrimaryAction
					->SetPositionAndDirection("PLAID",
							     	  altitude[i],
							          longitude[i],
							          latitude[i], 
	                              			          zenith,
							          azimuth);
		G4ThreeVector PLAPosition = 
				fUserPrimaryAction->GetPLAPosition();
	        std::vector<double> ILat;	
		ILat.push_back(0);
		ILat.push_back(0);
		ILat.push_back(0);
		if (test1){ 
         		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
          		if (!TotalDurationReached){
				theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude[i],longitude[i],ILat);
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(theBfield->GetReferenceDate().SpenvisModifiedJulianDay());
					values.push_back(altitude[i]/km);
					values.push_back(latitude[i]/degree);
					values.push_back(longitude[i]/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
		
	}
	if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);						     
  	theAnalysisManager->CloseAsciiFile();
	fUserEventAction->SetDrawTrajectory(draw_trajectory);
	clock2=clock();
    	double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  	G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
	
  }
  else {
  	G4cout<<"A problem occured when reading the CSV file."<<std::endl; 
	G4cout<<"The rigidity cutoffs will not be computed."<<std::endl;
	
  }
  ResetUserActions();
 
  
		    
}
///////////////////////////////////////////////////////////////////////////////////////////
// 
void MagneticShieldingTool::RCutoffVsSpenvisTrajectory(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    			     G4double zenith, 
				    			     G4double azimuth)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  
  //Read the position grid
  
  std::vector<double> altitude;
  std::vector<double>  latitude;
  std::vector<double>  longitude;
  std::vector<double>  smjd;
  
  
  G4bool test = theSpenvisManager
  			->ReadSpenvisTrajectoryCSVFile(SpenvisCSVFileName,
							 altitude,
							 latitude,
							 longitude,
							 smjd,
							 1);
  //change user actions
  SetActions();	   
  if (test) {
  	
	
	DateAndTime theDate = DateAndTime(2000,1,1,0,0,0);
	DateAndTime StartDate = DateAndTime(2000,1,1,0,0,0);
	StartDate.ConvertToSpenvisModifiedJulianDay(smjd[0]);
	

	//Open the ascii file
  	MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  	theAnalysisManager->OpenCutoffVsSpenvisPositionGridFile(outputfile_name,
								"PLAID",
	                                        		zenith, 
								azimuth);
 	//cvs file
  	if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->InitialiseCutoffVsPosCSVFile("PLAID",
								zenith, 
								azimuth);
 
   
  	//avoid drwaing particles
	G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
  	fUserEventAction->SetDrawTrajectory(false);	  
  
  	theBfield->SetStartDate(StartDate);
	G4cout<<altitude.size()<<std::endl;
	for (unsigned int i=0; i<altitude.size();i++){
		theDate.ConvertToSpenvisModifiedJulianDay(smjd[i]);
		G4double t =theDate.DifferenceInSeconds(StartDate);
		theBfield->SetTimeOfB(t);
	 	G4double test1 = fUserPrimaryAction
					->SetPositionAndDirection("PLAID",
							     	  altitude[i],
							          longitude[i],
							          latitude[i], 
	                              			          zenith,
							          azimuth);
		G4ThreeVector PLAPosition = 
				fUserPrimaryAction->GetPLAPosition();
	        std::vector<double> ILat;	
		ILat.push_back(0);
		ILat.push_back(0);
		ILat.push_back(0);
		if (test1){ 
         		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude[i],longitude[i],ILat);
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(smjd[i]);
					values.push_back(altitude[i]/km);
					values.push_back(latitude[i]/degree);
					values.push_back(longitude[i]/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
		
	}
	if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);						     
  	theAnalysisManager->CloseAsciiFile();
	fUserEventAction->SetDrawTrajectory(draw_trajectory);
	clock2=clock();
    	double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  	G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
	
  }
  else {
  	G4cout<<"A problem occured when reading the CSV file."<<std::endl; 
	G4cout<<"The rigidity cutoffs will not be computed."<<std::endl;
	
  }
  ResetUserActions();
 
}

////////////////////////////////////////////////////////////////////////////////
//
/*void MagneticShieldingTool::RCutoffVsPositionOnDipoleMagneticShell
                          (G4String CoordSys, G4double L,
                           G4double lat0,G4double delta_lat, G4int nlat,
                           G4double long0, G4double delta_long, G4int nlong,
			   G4double zenith,G4double azimuth,
			   G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  //Open the ascii file
  MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  theAnalysisManager->OpenCutoffVsPositionOnLShellFile
                       (outputfile_name,CoordSys, L, zenith, azimuth);


  // avoid drawing particle trajectory
  G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
  fUserEventAction->SetDrawTrajectory(false);
  
  
  // compute grid
  for (int i=0; i<nlat;i++){
  	G4double latitude=lat0+double(i)*delta_lat;
   	for (int j=0; j<nlong;j++){
      		G4double longitude=long0+double(j)*delta_long;
       		G4double test = fUserPrimaryAction
	    				->SetPositionOnDipoleMagneticShell
	                        		(CoordSys,L,latitude,longitude);
		G4ThreeVector PLAPosition = 
				fUserPrimaryAction->GetPLAPosition();
		std::vector<double> ILat;	
		ILat.push_back(0);
		ILat.push_back(0);
		ILat.push_back(0);			      
       		if (test){ 
         		fUserPrimaryAction
	                 	->SetDirection(CoordSys, zenith, azimuth);
	  		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude,
						longitude,ILat);
			}		
	  	}
       		else {
         		j=nlong;
	  		i=nlat;
	 	}
       }
  }
  theAnalysisManager->CloseAsciiFile();   
  
  fUserEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;
}
*/
/////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::RCutoffVsDirection
                          (G4String CoordSys,G4double zen0, G4double delta_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim,
			   G4String outputfile_name)
{ 
if (TotalDurationReached) 
	{
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
	}
clock_t clock1,clock2;
clock1 = clock();
G4cout<<"Dir1 "<<std::endl;
//Open the ascii file
MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  
//PVD
  
//----------------------
  
ifstream in(outputfile_name.c_str());

std::vector<double> vector_zenith;
std::vector<double> vector_azimuth;
std::vector<double> vector_Ru;
std::vector<double> vector_Rc;
std::vector<double> vector_Rl;	

G4double Oldx = -1000, Oldy = -1000, Oldz = -1000;

for(double zenith = zen0; zenith < zen0 + delta_zen * (nzen-1); zenith  = zenith + delta_zen)
	{
	for(double azimuth = azim0; azimuth < azim0 + delta_azim * (nazim-1); azimuth  = azimuth + delta_azim)
		{
		vector_zenith.push_back(zenith);
		vector_azimuth.push_back(azimuth);
		vector_Ru.push_back(1e10);
		vector_Rc.push_back(1e10);
		vector_Rl.push_back(1e10);
		}
	}
	
int ctr_missing = 0;

if(in)
	{
	//read header of ascii files
		
	while(!in.eof())
		{
		G4String Text;
		in>>Text;
		if(!Text.compare("Position:")) break;
		}

	if(!in.eof())
		{
		in>>Oldx;
		in>>Oldy;
		in>>Oldz;
		}
		
	while(!in.eof())
		{
		G4String Text;
		in>>Text;
		if(!Text.compare("Rl")) break;
		}		
	
	//read cutoffs
			
	double zenith, azimuth, Ru, Rc, Rl;
	
	while(1)
		{	
		if(!in.good()) break;
		else
			{
			in>>zenith>>azimuth>>Ru>>Rc>>Rl;
			
			for(unsigned int i = 0; i < vector_zenith.size(); i++)
				{	
				if(fabs(vector_zenith.at(i)-zenith*degree) < 1e-4 && fabs(vector_azimuth.at(i)-azimuth*degree) < 1e-4)
					{
					vector_Ru.at(i) = Ru;
					vector_Rc.at(i) = Rc;
					vector_Rl.at(i) = Rl;
				
					break;
					}
				}
			}
		}
	}
in.close();  
  
//----------------------

G4double Currentx = 1000, Currenty = 1000, Currentz = 1000;
if (CoordSys!="PLAG")
	{
	PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator = (PLANETOCOSPrimaryGeneratorAction*) G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();  
	G4ThreeVector PLAposition = thePrimaryGenerator->GetPLAPosition();
	G4ThreeVector position = SpaceCoordinatePlanet::GetInstance()->Transform(PLAposition,"PLA",CoordSys)/PlanetManager::GetInstance()->GetRplanet();
		
	Currentx = position.x();
	Currenty = position.y();
	Currentz = position.z();
	}

if(fabs(Oldx-Currentx) > 1e-3 || fabs(Oldy-Currenty) > 1e-3 || fabs(Oldz-Currentz) > 1e-3)
	{
	for(unsigned int i = 0; i < vector_zenith.size(); i++) 
		{
		vector_Ru.at(i) = 1e10;
		vector_Rc.at(i) = 1e10;
		vector_Rl.at(i) = 1e10;
		}
	}
	
theAnalysisManager->OpenCutoffVsDirectionFile(outputfile_name,CoordSys, false);

G4cout<<"Dir2 "<<std::endl;
//change user actions
SetActions();	

//initialise the cvs file
if (RegisterResultsInSpenvisFile) 
	{
	G4ThreeVector thePLAPosition = fUserPrimaryAction->GetPLAPosition();
	SpaceCoordinatePlanet* theCoordinateConvertor = SpaceCoordinatePlanet::GetInstance();  
	G4double altitude,longitude,latitude;
	if (CoordSys == "PLAG" || CoordSys == "GEODETIC") theCoordinateConvertor->ComputePLAGCoordinatesFromPLAPosition(thePLAPosition, altitude, longitude, latitude);
	else 
		{
		G4ThreeVector thePosition = theCoordinateConvertor->Transform(thePLAPosition, "PLA", CoordSys);
		altitude = thePosition.mag()-Re;
		latitude = 90.*degree -thePosition.theta();

		longitude = thePosition.phi();
		if (longitude > 180. *degree) longitude = longitude -360.*degree;		   
		}
	theSpenvisManager->InitialiseCutoffVsDirCSVFile(CoordSys, altitude, latitude, longitude);
	}  
G4cout<<"Dir3 "<<std::endl;
//avoid particels drawing
G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
fUserEventAction->SetDrawTrajectory(false);
G4cout<<"Dir4 "<<std::endl;

for(unsigned int i = 0; i < vector_zenith.size(); i++)
	{
	G4double zenith = vector_zenith.at(i);
	G4double azimuth = vector_azimuth.at(i);
	
	if(vector_Ru.at(i) == 1e10 || vector_Rc.at(i) == 1e10 || vector_Rl.at(i) == 1e10)
		{
		G4double test = fUserPrimaryAction->SetDirection(CoordSys,zenith,azimuth);
		G4cout<<"Dir5 "<<std::endl;
		if(test)
			{
			if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
			else ComputeRigidityCutoff();
			if (!TotalDurationReached)
				{
				theAnalysisManager->RegisterCutoffVsDirection(Rc,Rm,Rs,zenith,azimuth);
				
				if (RegisterResultsInSpenvisFile)
					{
					std::vector<double > values;
					values.clear();
					values.push_back(zenith/degree);
					values.push_back(azimuth/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
					}
				}	
			}         
		}
	else
		{
		Rc = vector_Rc.at(i);	
		Rm = vector_Ru.at(i);
		Rs = vector_Rl.at(i);
			
		theAnalysisManager->RegisterCutoffVsDirection(Rc,Rm,Rs,zenith,azimuth);
				
		if (RegisterResultsInSpenvisFile)
			{
			std::vector<double > values;
			values.clear();
			values.push_back(zenith/degree);
			values.push_back(azimuth/degree);
			values.push_back(Rm);
			values.push_back(Rc);
			values.push_back(Rs);
			theSpenvisManager->WriteData(values);
			}
		}
	}
theAnalysisManager->CloseAsciiFile();
if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);

fUserEventAction->SetDrawTrajectory(draw_trajectory);
ResetUserActions();
clock2=clock();
double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
}  
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::RCutoffVsTime
                         (G4double time0, G4double delta_t, G4int ntime, 
			  G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  clock1=clock();
  
  //change user actions
  SetActions();	  
  
  //Open the ascii file
  MagneticShieldingToolAnalysisManager* theAnalysisManager= MagneticShieldingToolAnalysisManager::GetInstance();
  theAnalysisManager->OpenCutoffVsTimeFile(outputfile_name);

   
   //initialise the cvs file
  
  if (RegisterResultsInSpenvisFile) { 
  
       G4ThreeVector thePLAPosition = fUserPrimaryAction->GetPLAPosition();
	G4ThreeVector thePLADirection =  fUserPrimaryAction->GetPLADirection();
  	G4double altitude,longitude,latitude, zenith,azimuth;
	altitude = thePLAPosition.mag()-Re;
	latitude = 90.*degree - thePLAPosition.theta();
	longitude = thePLAPosition.phi();
	if (longitude > 180.*degree) longitude += -360.*degree; 
	thePLADirection.rotateZ(-thePLAPosition.phi());
	thePLADirection.rotateY(-thePLAPosition.theta());
	thePLADirection=-thePLADirection;
	zenith = thePLADirection.theta();
	azimuth = thePLADirection.phi();
	
	theSpenvisManager->InitialiseCutoffVsTimeCSVFile("PLA", 
   				     			altitude,
				     			latitude,
				     			longitude,
							zenith,
							azimuth);
  }  
   
  
  //avoid drawing particles
  G4bool draw_trajectory = fUserEventAction->GetDrawTrajectory();
  fUserEventAction->SetDrawTrajectory(false);
  
  for (int i=0; i<ntime;i++){
  	G4double time=time0+double(i)*delta_t;
   	theBfield->SetTimeOfB(time/s);
   	//theBfield->PrintStormParameter();
   	if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
   	else ComputeRigidityCutoff();
	if (!TotalDurationReached){
   		theAnalysisManager
	     		->RegisterCutoffVsTime(Rc,Rm,Rs,time);
		if (RegisterResultsInSpenvisFile){
				std::vector<double > values;
				values.clear();
				values.push_back(theBfield->GetReferenceDate().SpenvisModifiedJulianDay());
				values.push_back(Rm);
				values.push_back(Rc);
				values.push_back(Rs);
				theSpenvisManager->WriteData(values);
		}
	}	
  }

  theAnalysisManager->CloseAsciiFile();
  if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);
  fUserEventAction->SetDrawTrajectory(draw_trajectory);
  
  ResetUserActions();
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;   

}
///////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::RegisterTrackLastPoint(G4ThreeVector position, 
                             G4ThreeVector momentum,G4double filter_value)
{ LastPositions.push_back(position);
  LastMomentaOnCharge.push_back(momentum);
  FilterValues.push_back(G4int(filter_value));
 
  if (RegisterAsymptoticDirection){
    	MagneticShieldingToolAnalysisManager::GetInstance()
                ->RegisterAsymptoticDirection(position,momentum,
		                                  G4int(filter_value));
  	if (RegisterResultsInSpenvisFile) 
	    theSpenvisManager->RegisterAsymptoticDirection(position,momentum,
		                                           G4int(filter_value));
  
  } 
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ComputeRigidityCutoff()
{ if (TotalDurationReached) return;
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();

 // Switch on the field
  theBfield->SwitchOn();

    //change user actions
 // SetActions();

  //No detection of secondary particles
 // PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(false); 

  //Set the  equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);

  //get a pointer to the primary generator action
  fUserPrimaryAction->SelectTypeOfPrimaries("RigidityFilter");
  fUserPrimaryAction->ResetRigidityIndex();
 
  //pointer on PLANETOCOSTrackingAction
  fUserTrackingAction->SetRegisterLastPoint(true);

  //Compute trajectory and rigidity cutoff

  G4int n_particles = fUserPrimaryAction->GetNumberOfRigidity();  
  if (n_particles> 0){
    	G4RunManager::GetRunManager()->BeamOn(n_particles);
    
    	// RM RS and RC computation
    	G4double Rss=0,Rmm=0;
    	G4double Rcc=0.;
    	G4double last_rigidity=100.,fac=0.,last_filter_value=1.;
    	for (int j=0; j<n_particles;j++){
      		G4double rigidity = LastMomentaOnCharge[j].mag()/1000.;
       		G4int filter_value = FilterValues[j];
		if (filter_value <0) filter_value =0;
       		if ( (fac==0.) && filter_value==0){
        		fac=1.;
         		Rmm=last_rigidity;
         	}
       		Rcc += fac*filter_value*(rigidity-last_rigidity);
       		last_rigidity = rigidity;
      		if ( filter_value==1) Rss=rigidity;
      		last_filter_value=filter_value;
      	}
     
     	Rcc+=Rmm;
     	Rc=Rcc;
     	Rs=Rss;
     	Rm=Rmm;
      
  }
  else G4cout<<"The rigidity vector is empty"<<G4endl;
  if (Rm==Rc) Rs=Rc;
  
 
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
 
 
  //go back to normal conditions
 
  theEquation->SetReverseTimeMode(false);
  fUserPrimaryAction->SelectTypeOfPrimaries("Standard");
  fUserTrackingAction->SetRegisterLastPoint(false); 
  PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(true); 
 // ResetUserActions();
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ComputeCutoff()
{ 
 //In this method the cutoff rigidity is computed in the following way:
 // a) A start rigidity Sc is computed from the Störmer formula by considering the 
 //   geomagnetic  dipole
 // b) It is looked in the following rigidity series [Sc,Sc+1.GV,Sc+2.GV,....]
 // 	for the first rigidity with an allowed trajectory. This rigidity is called R1st .
 // c)It is looked in the serie [R1st-0.01GV, R1st-0.02GV,...,R1st-n*0.01GV,.... ]
 //    for the first rigidity with a forbidden trajectory. We call this rigidity R1stForb.
 // d)We compute trajectories for [R1st +0.01GV, R1st+0.02GV,...,R1st+n*0.01GV,.... ]
 //   till we have find a 1GV  band of rigidity with only allowed trajectories. Note that
 // the band [R1stForb+0.01GV, R1st] is considered as allowed and if a serie of rigidity
 // is found directly above [R1stForb+0.01GV, R1st] such that together they form a 1GV
 // band of allowed trajectories the 1GV band is considered as being found.
 //  We define Rm as the lowest rigidity in the 1GV allowed band.
 // e)It is looked in the serie [R1stForb, R1stForb-0.01GV,...,R1stForb-n*0.01GV,.... ]
  //   to the first 1GV band with all trajectories forbidden. Rs is defined by the 
  // highest rigidity of this band +0.01 GV.
  // f) Nallowed the number of allowed trajectory betwenn
  //   Rm and Rs-0.01 GV is computed during a) to e). 
  //   Rc =Rm- Nallowed*0.01*GV
  	

  if (TotalDurationReached) return;
  
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
  
  // Switch on the field
  theBfield->SwitchOn();
  
  
  
  //No detection of secondary particels
 // PLANETOCOSAnalysisManager::GetInstance()->SetCallEventActions(false); 
  
 
  //Set the  equation of motion
  PlanetEquationOfMotion* theEquation =  theBfield->GetEquationOfMotion();
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);
  
 
	  
	  
  //PLAPosition PLAdirection
  G4ThreeVector PLAPosition = fUserPrimaryAction->GetPLAPosition();
  G4ThreeVector PLADirection = fUserPrimaryAction->GetPLADirection();

  //PLADipole
  G4double DipoleB0=theBfield->GetDipoleB0();
  G4double DipolePhi=theBfield->GetDipolePhi();
  G4double DipoleTheta=theBfield->GetDipoleTheta();
  G4ThreeVector DipoleSchift=theBfield->GetDipoleShift();

  // compute the Störmer cutoff in a geodipole field 	  
  G4ThreeVector PositionInDipole=PLAPosition-DipoleSchift;
  PositionInDipole.rotateZ(-DipolePhi).rotateY(-DipoleTheta);
  G4double phi=PositionInDipole.phi();
  G4double r=PositionInDipole.r()/re;
  G4double cos_geomagnetic_lat=sin(PositionInDipole.theta());
  G4double cos3=cos_geomagnetic_lat*cos_geomagnetic_lat*cos_geomagnetic_lat;
  G4double cos4=cos3*cos_geomagnetic_lat;
  G4ThreeVector EastDirection=G4ThreeVector(-sin(phi),cos(phi),0); 
  G4ThreeVector DirectionInDipole=PLADirection;
  DirectionInDipole.rotateZ(-DipolePhi).rotateY(-DipoleTheta);
  G4double cos_gam=-EastDirection.dot(DirectionInDipole);
 
  // G4cout<<cos_gam<<" cos gam"<<endl;
  G4double c=2.99792458e8*m/s;
  G4double term=1.+std::sqrt(1.-cos3*cos_gam);
  G4double SCutoff=DipoleB0*re*c*cos4/(r*r*term*term);
  SCutoff=G4int(SCutoff/(0.01*GV))*(0.01*GV);
  SCutoff = std::max(SCutoff, 0.01*GV);
  //G4cout<<SCutoff/GV<<G4endl;

  //pointer on PLANETOCOSTrackingAction
  PLANETOCOSTrackingAction* theTrackingAction =
       (PLANETOCOSTrackingAction*)
               G4RunManager::GetRunManager()->GetUserTrackingAction();
  theTrackingAction->SetRegisterLastPoint(true);
  

  G4int nmax=100;
  G4int n=0;
  G4double rigidity=SCutoff;
  G4double first_rigidity=SCutoff;
  G4double rig_min = .01*GV;
  G4double rig_max = std::max(100. * GV, SCutoff+10.*GV);
  G4int n_forbidden=0;
  G4int n_allowed=0;
  //G4cout<<SCutoff/GV<<" GV SCutoff"<<std::endl;
 //find a first allowed trajectory
 //---------------------------

  first_rigidity-=1.*GV;
  G4int filter_value= 0;
while (filter_value ==0){
first_rigidity+= 1.*GV;
fUserPrimaryAction->SetRigidity(first_rigidity);
G4RunManager::GetRunManager()->BeamOn(1);

if (TotalDurationReached) return;
filter_value = FilterValues[0];
if (filter_value <0) filter_value =0;

}
//G4cout<<"first "<<first_rigidity/GV<<std::endl;
Rm=first_rigidity;  
  
  //Find the first forbidden trajectory below first_rigidity 
  G4double first_forbid_rigidity=first_rigidity;
  while (filter_value ==1 &&  first_forbid_rigidity>= rig_min+0.01*GV){
  	first_forbid_rigidity-=0.01*GV;
	fUserPrimaryAction->SetRigidity(first_forbid_rigidity);
	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;
	
  }
  if (filter_value ==0) {
  	Rm=first_forbid_rigidity+0.01*GV;
  }
  else {
  	Rm=first_forbid_rigidity;
  	Rm/=GV;
  	Rc=Rm;
  	Rs=Rm;
  	LastPositions.clear();
 	 LastMomentaOnCharge.clear();
  	FilterValues.clear();
 
 
  	//go back to normal 
  	theEquation->SetReverseTimeMode(false);
  	fUserPrimaryAction->SelectTypeOfPrimaries("Standard");
  	theTrackingAction->SetRegisterLastPoint(false);
	
	if (Rm<0.015 *GV){
		Rm=0.;
		Rc=0.;
		Rs=0.;
  	}
	return;
  }
  Rs=Rm;
  
  
  // find first 1GV band of allowed trajectories above first_forbid_rigidity 
  //------------------------------------------------------------------------ 
  
  n=int((first_rigidity-Rm)/(0.01*GV))+1;
  //G4cout<<"Rm "<<Rm/GV<<std::endl;
  //G4cout<<"n "<<n<<std::endl;
  rigidity= first_rigidity;
  
  while (n< nmax && rigidity <= rig_max-.01*GV) {
     	rigidity += .01*GV;
     	//G4cout<<"rigidity "<<rigidity<<std::endl;
      	fUserPrimaryAction->SetRigidity(rigidity);
      	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;

      	filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;
      	n++;
     
    	//G4cout<<filter_value<<" filter"<<std::endl;
      
      	if (filter_value == 0 ){ 
        	n_forbidden+=1;
		Rm=rigidity + 0.01*GV;
	    	n=0;
	}      		 	 
  }
  
  //G4cout<<"Rm "<<Rm/GV<<std::endl;
  Rs=first_forbid_rigidity+0.01*GV;
  
  
  //find first 1GV band of forbidden trajectories below first_forbid_rigidity
  //---------------------------------
  
  n_allowed=int ((Rm- first_forbid_rigidity)/(0.01*GV));
   //we add 1 at n_forbidden because the lowest computed rigidity is forbidden
  n_forbidden+=1;
  n_allowed-=n_forbidden;
  //G4cout<<"n_allowed "<<n_allowed<<std::endl;
 
  n=1;
  rigidity=first_forbid_rigidity;
  while (n< 100 && rigidity >= rig_min + 0.01*GV ){ 
     	rigidity -= .01*GV;
      	fUserPrimaryAction->SetRigidity(rigidity);
      	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	G4int filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;    
      	n++;
      	if (filter_value == 1 ){ 
        	n=0;
	  	n_allowed+=1; 
          	Rs=rigidity;
		//G4cout<<"new Rs :"<<Rs<<std::endl;
	}	 	 
  }
  //G4cout<<"n_allowed "<<n_allowed<<std::endl;
  //G4cout<<Rm/GV<<std::endl;
  if (Rm<0.015 *GV){
	Rm=0.;
	Rc=0.;
	Rs=0.;
  }
  
  else {
  	Rc=Rm-n_allowed*.01*GV; 
  }
  Rc/=GV;
  Rs/=GV;
  Rm/=GV;
  //G4cout<<Rc<<" "<<Rs<< " "<<Rm<<std::endl;
  
  if (Rm==Rs) Rc=Rs;
  
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
 
 
  //go back to normal 
  theEquation->SetReverseTimeMode(false);
  fUserPrimaryAction->SelectTypeOfPrimaries("Standard");
  theTrackingAction->SetRegisterLastPoint(false);
  
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::ResetMaxStepLength()
{ theGeometry->SetAtmosphereMaxStepLength(AtmosphereMaxStepLength);
  theGeometry->SetMagnetosphereMaxStepLength(MagnetosphereMaxStepLength);
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::SwitchOffAllProcessesExceptTransportation()
{ G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  std::vector<G4String>* theProcessNames = theProcessTable->GetNameList();
  for (unsigned int i=0; i<theProcessNames->size();i++){
	theProcessTable->SetProcessActivation((*theProcessNames)[i],false);
  } 
  theProcessTable->SetProcessActivation("Transportation",true);
  theProcessTable->SetProcessActivation("StepLimiter",true);
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::SwitchOnAllProcesses()
{ G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  std::vector<G4String>* theProcessNames = theProcessTable->GetNameList();
  for (unsigned int i=0; i<theProcessNames->size();i++){
	theProcessTable->SetProcessActivation((*theProcessNames)[i],true);
  } 
}
////////////////////////////////////////////////////////////////////////////////
//
void  MagneticShieldingTool::SetActions(){
 G4RunManager* theRunManager =  G4RunManager::GetRunManager();
 fUserPrimaryAction=dynamic_cast<PLANETOCOSPrimaryGeneratorAction*>(
 		    	const_cast<G4VUserPrimaryGeneratorAction* >(
				theRunManager->GetUserPrimaryGeneratorAction()
			)
		    );
 
 fUserTrackingAction=dynamic_cast<PLANETOCOSTrackingAction*>(
 		    	const_cast<G4UserTrackingAction* >(
				theRunManager->GetUserTrackingAction()
			)
		    );
 fUserEventAction=dynamic_cast<PLANETOCOSEventAction*>(
 		    	const_cast<G4UserEventAction* >(
				theRunManager->GetUserEventAction()
			)
		    );
 fEventAction->SetUserEventAction(fUserEventAction);
 fUserSteppingAction=dynamic_cast<PLANETOCOSSteppingAction*>(
 		    	const_cast<G4UserSteppingAction* >(
				theRunManager->GetUserSteppingAction()
			)
		    );
 fSteppingAction->SetUserSteppingAction(fUserSteppingAction);
 fUserRunAction=dynamic_cast<PLANETOCOSApplicationScenario*>(
 		    	const_cast<G4UserRunAction*>(
				theRunManager->GetUserRunAction()
			)
		    ); 
 //Replace the user action by the BlineTool actions
  
 
  theRunManager->SetUserAction(this);
  theRunManager->SetUserAction(fSteppingAction);
  theRunManager->SetUserAction(fEventAction);
}
///////////////////////////////////////////////////////////////////////////////
//
void  MagneticShieldingTool::ResetUserActions(){
  //Reset the user actions
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  theRunManager->SetUserAction(fUserRunAction);
  theRunManager->SetUserAction(fUserSteppingAction);
  theRunManager->SetUserAction(fUserEventAction);
}
//////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingTool::InitialiseForMagneticShieldingStudy()
{ theGeometry = 
   	const_cast<PLANETOCOSGeometryConstruction*>
		(dynamic_cast<const PLANETOCOSGeometryConstruction*>(
				G4RunManager::GetRunManager()
					->GetUserDetectorConstruction()));
 theGeometry->SetConsiderAtmosphere(false);
 PlanetManager::GetInstance()->GetSoil()->ResetLayers();
 PLANETOCOSPhysicsList* pPhysicsList  = 
   	const_cast<PLANETOCOSPhysicsList*>
		(dynamic_cast<const PLANETOCOSPhysicsList*>(
				G4RunManager::GetRunManager()
					->GetUserPhysicsList()));
 
 pPhysicsList->SelectEMPhysicsType("NONE");
 pPhysicsList->SelectHadronicPhysicsType("NOHADRONIC"); 
 pPhysicsList->BuildList();
 G4RunManager::GetRunManager()->Initialize();					
  					
}
