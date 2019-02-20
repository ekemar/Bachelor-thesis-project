#ifndef MagneticShieldingTool_h
#define MagneticShieldingTool_h 1

#include "G4ThreeVector.hh"
#include "G4UserRunAction.hh"
#include<vector>
class PLANETOCOSSpenvisManager;
class MagneticShieldingToolMessenger;
class MagneticShieldingToolSteppingAction;
class MagneticShieldingToolEventAction;
class PlanetMagneticField;
class PLANETOCOSPrimaryGeneratorAction;
class PLANETOCOSTrackingAction;
class PLANETOCOSSteppingAction;
class PLANETOCOSStackingAction;
class PLANETOCOSApplicationScenario;
class PLANETOCOSEventAction;
class PLANETOCOSGeometryConstruction;

class MagneticShieldingTool : public G4UserRunAction 
{ public:
  
   MagneticShieldingTool();
   virtual ~MagneticShieldingTool();
   
   static MagneticShieldingTool* GetInstance();
  
  public:
  
   virtual void BeginOfRunAction(const G4Run* aRun);
   virtual void EndOfRunAction(const G4Run* aRun);
  
// possible applications
   void Bline();
   void ParticleTrajectoryInBfield(G4int n);
   void ReverseParticleTrajectoryInBfield(G4int n);
   void ComputeRigidityFilter(G4String outputfile_name);
   void RCutoffVsPosition(G4String CoordSys, 
                          G4double Altitude,
                          G4double lat0,G4double delta_lat, G4int nlat,
                          G4double long0, G4double delta_long, G4int nlong,
			  G4double zenith, G4double azimuth, 
			  G4String outputfile_name);
   void RCutoffVsSpenvisPositionGrid(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    G4double zenith =0., 
				    G4double azimuth = 0.);			  
   void RCutoffVsSpenvisTrajectory(G4String SpenvisCSVFileName,G4String outputfile_name,
  				  G4double zenith =0., 
				  G4double azimuth = 0.);			  
   void RCutoffVsDirection(G4String CoordSys, G4double zen0, G4double delta_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim, 
			   G4String outputfile_name);
   
   void RCutoffVsTime(G4double time0, G4double delta_t, G4int ntime, 
		     G4String outputfile_name);			  
 
    
    
    
   void RegisterTrackLastPoint(G4ThreeVector position, 
                             G4ThreeVector momentum,G4double filter_value);
       
   void SetAutomaticDetectionPenumbra(G4bool aBool);

   inline void SetRegisterResultsInSpenvisFile(bool aBool) 
  		{RegisterResultsInSpenvisFile = aBool;}		 
    
   inline void SetSpenvisFileName(G4String aName){SpenvisFileName = aName;}
   
  
  
  
  
   bool CheckIfParticleShouldBeKilledAndSetMaxStepLength(G4ThreeVector pos);
//set methods
   void SetStopAltitude(G4double aVal){stop_altitude=aVal;}
   void ResetMaxStepLength();
 //get methods
   inline   MagneticShieldingToolEventAction* GetEventAction(){return fEventAction;};
   
   void InitialiseForMagneticShieldingStudy();
   
   		       
 private:
 //private methods
  // bool CheckMagneticFields();
   void SwitchOffAllProcessesExceptTransportation();
   void SwitchOnAllProcesses();
   void ComputeRigidityCutoff();
   void ComputeCutoff();
   void SetActions();
   void ResetUserActions(); 
 
 private:
 
   static MagneticShieldingTool* instance;
   PLANETOCOSSpenvisManager* theSpenvisManager;
   
   MagneticShieldingToolMessenger* fMessenger;
   MagneticShieldingToolSteppingAction* fSteppingAction;
   MagneticShieldingToolEventAction* fEventAction;
   
   //user defined primary generator action
   PLANETOCOSPrimaryGeneratorAction* fUserPrimaryAction;
   PLANETOCOSApplicationScenario* fUserRunAction;
   PLANETOCOSStackingAction* fUserStackingAction;
   PLANETOCOSSteppingAction* fUserSteppingAction;
   PLANETOCOSTrackingAction* fUserTrackingAction;
   PLANETOCOSEventAction* fUserEventAction;

    //Last position, momentum and filter value for a particle trajectory 
   std::vector<G4ThreeVector>  LastPositions;
   std::vector<G4ThreeVector>  LastMomentaOnCharge;
   std::vector<G4int>  FilterValues;
   
   
   //Rigidity cutoff
    G4double  Rs;
    G4double  Rc;
    G4double  Rm;
    
   // run duration 
    bool TotalDurationReached;  
    
    bool RegisterAsymptoticDirection;
    bool AutomaticDetectionOfThePenumbra;
    G4double MagnetosphereMaxStepLength ;
    G4double AtmosphereMaxStepLength;
    G4String geometry_type;
    G4double ZpositionAtZeroAltitude;
    G4double RplanetAtEquator;
    G4double stop_altitude;
    PLANETOCOSGeometryConstruction* theGeometry;
    PlanetMagneticField* theBfield;
    G4bool BlineMode;
    
    
    //Spenvis
    bool RegisterResultsInSpenvisFile;
    G4String SpenvisFileName;
    

    
  
};












#endif
