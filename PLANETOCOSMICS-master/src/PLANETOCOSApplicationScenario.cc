#include"PLANETOCOSApplicationScenario.hh"
#include"PLANETOCOSScenarioMessenger.hh"
#include"PLANETOCOSPrimaryGeneratorAction.hh"
#include"PLANETOCOSEventAction.hh"
#include"PLANETOCOSTrackingAction.hh"
#include"PLANETOCOSSteppingAction.hh"
#include"PlanetEquationOfMotion.hh"
#include"PlanetMagneticField.hh"
#include"PLANETOCOSAnalysisManager.hh"
#include"PLANETOCOSMagneticShieldingAnalyser.hh"
#include"G4TransportationManager.hh"
#include"G4FieldManager.hh"
#include"G4UImanager.hh"
#include"G4ParticleDefinition.hh"
#include"G4Circle.hh"
#include"G4Event.hh"
#include"G4UnitsTable.hh"
#include"G4ios.hh"
#include"G4RunManager.hh"
#include"G4ProcessTable.hh"
#include"G4GeneralParticleSource.hh"
#include"G4Proton.hh"
#include"G4ChargedGeantino.hh"
#include"PLANETOCOSGeometryConstruction.hh"
#include"DurationManager.hh"
#include<time.h>

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSApplicationScenario::PLANETOCOSApplicationScenario()
{ theMessenger = new PLANETOCOSScenarioMessenger(this);
  TotalDurationReached = false;
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSApplicationScenario::~PLANETOCOSApplicationScenario()
{ delete theMessenger;
}		  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSApplicationScenario::BeginOfRunAction(const G4Run*) // aRun
{ DurationManager* theDurationManager = DurationManager::GetInstance();
  TotalDurationReached = !theDurationManager->CheckDurationAtBeginOfRun();
  if (TotalDurationReached) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  }
  
  if (RandomSeedNeeded)  SetRandomSeed();
  
#ifdef ANALYSIS_SOIL
  const PLANETOCOSSteppingAction* 
			    theSteppingAction = 
		               dynamic_cast<const PLANETOCOSSteppingAction*>
  				  (G4RunManager::GetRunManager()->GetUserSteppingAction());  
  const_cast<PLANETOCOSSteppingAction*>(theSteppingAction)->ResetLowestZorRReached();	
#endif  
}  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSApplicationScenario::EndOfRunAction(const G4Run*)  //  aRun
{ /*const PLANETOCOSSteppingAction* 
			    theSteppingAction = 
		               dynamic_cast<const PLANETOCOSSteppingAction*>
				  (G4RunManager::GetRunManager()->GetUserSteppingAction());
			
  G4double lowest_altitude_needed;
  lowest_altitude_needed = theSteppingAction->GetLowestAltitudeNeeded();
  //G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
  */
#ifdef ANALYSIS_SOIL
  G4double LowestZorRReached= theSteppingAction->GetLowestZorRReached();
  const PLANETOCOSGeometryConstruction* 
			    theGeometry = 
		               dynamic_cast<const PLANETOCOSGeometryConstruction*>
				  (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double LowestDepthSoilReached = 
    (theGeometry->GetZpositionForAtmosphereBottom() - LowestZorRReached	)/m;			
  G4cout<<"Lowest soild depth reached "<<LowestDepthSoilReached<<" m"<<std::endl;
#endif 
			
}


////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSApplicationScenario::SetRandomSeed()
{ 
  //-------original code ---------
  //time_t time1=time(0);
  //int second= (*(gmtime(&time1))).tm_sec;
  //int minute= (*(gmtime(&time1))).tm_min;
  //int hour= (*(gmtime(&time1))).tm_hour;
  //int day_year = (*(gmtime(&time1))).tm_yday;
 
  //long nextSeed = (long) (day_year + second + minute + hour);
  //HepRandom::setTheSeed(nextSeed );
  //nextSeed = (long) (100000000L * HepRandom::getTheGenerator()->flat());
  //nextSeed +=  (long) (day_year + second + minute + hour);
  //HepRandom::setTheSeed(nextSeed );
  //G4cout<<" seed value :"<< nextSeed <<G4endl;
  //----------------------------

  //Set the seed as the time in nanoseconds since the epoch
  //AM
  struct timespec tsp;
  clock_gettime(CLOCK_REALTIME, &tsp); 
  long ts=tsp.tv_nsec + 1000000000*tsp.tv_sec;
  long nextSeed = (long) (ts);
  HepRandom::setTheSeed(nextSeed );
  G4cout<<" seed value :"<< nextSeed <<G4endl;
}



