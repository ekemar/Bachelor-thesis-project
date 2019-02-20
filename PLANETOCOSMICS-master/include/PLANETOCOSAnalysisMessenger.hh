#ifndef PLANETOCOSAnalysisMessenger_h
#define PLANETOCOSAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PLANETOCOSAnalysisManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;


class PLANETOCOSAnalysisMessenger: public G4UImessenger
{
  public:
    PLANETOCOSAnalysisMessenger(PLANETOCOSAnalysisManager * theManager);
    ~PLANETOCOSAnalysisMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    PLANETOCOSAnalysisManager*  pAnalysisManager;
    G4UIdirectory*             analysisDir;
    G4UIdirectory*             primfluxDir;
    G4UIdirectory*             detfluxDir;
    G4UIdirectory*             edepDir;
    G4UIdirectory*             pseudotrappingDir;
    
    G4UIdirectory*             completeEventDir;   
   
    G4UIcmdWithAnInteger* SelectDetectorCmd;
    G4UIcmdWithAnInteger* DeselectDetectorCmd;
    G4UIcmdWithoutParameter* SelectAllDetectorsCmd;
    G4UIcmdWithoutParameter* DeselectAllDetectorsCmd;
    
    //create edep histogram for atmosphere command
   
    G4UIcommand* CreateEdepVsDepthHistoCmd;
    G4UIcommand* CreateEdepVsAltitudeHistoCmd;
    
    //create edep histogram for soil command
   
    G4UIcommand* CreateSoilEdepVsDepthHistoCmd;
    G4UIcommand* CreateSoilEdepVsThicknessHistoCmd;
    
   
    //create flux histogram command
   
    G4UIcommand* CreateDownwardFluxVsEkinHistoCmd;
    G4UIcommand* CreateUpwardFluxVsEkinHistoCmd;
    G4UIcommand* CreateEnergyVsCosZenithHistoCmd;
    G4UIcommand* CreateCosZenithHistoCmd;
    G4UIcommand* CreateAzimuthHistoCmd;
    G4UIcommand* CreateDownwardFluxVsPosHistoCmd;
    G4UIcommand* CreateUpwardFluxVsPosHistoCmd;
    G4UIcommand* CreateDownwardFluxVsPosFlatHistoCmd;
    G4UIcommand* CreateUpwardFluxVsPosFlatHistoCmd;
    G4UIcommand* SetEnergyLimitsCmd;
    G4UIcommand* SetLatLongLimitsCmd;
    G4UIcmdWithAString* SetDivideByCosThCmd;
    
    
    //create primary histogram command
    G4UIcommand* CreateEnergyVsZenithPrimaryCmd;
    G4UIcommand* CreateZenithVsAzimuthPrimaryCmd;
    G4UIcommand* CreateOmniFluxPrimaryCmd;
    G4UIcommand* CreateZenithPrimaryCmd;
    G4UIcommand* CreateAzimuthPrimaryCmd; 
    
    // PVD Create root files with trees for primary and 
    G4UIcommand* CreateFluxRootFileCmd;
    G4UIcmdWithoutParameter* WriteFluxRootFileCmd;    
    
   //create pseudo tracking histogram command
    G4UIcommand*  CreateNbOfPlanetTurnVsStartEkinHistoCmd;
    G4UIcommand*  CreateNbOfPlanetTurnVsLifeTimeHistoCmd;
    G4UIcommand*  CreateNbOfEquatorCrossingVsStartEkinHistoCmd;
    G4UIcommand*  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd;
    G4UIcommand*  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd;
    G4UIcommand*  CreateLifeTimeVsStartEkinHistoCmd;
    
    //tree coommands
    G4UIcmdWithoutParameter* ResetHistoCmd;
    G4UIcmdWithoutParameter* ResetTreeCmd;
    G4UIcommand* SaveTreeCmd;
    G4UIcommand* ListTreeCmd;
    G4UIcommand* PrintTreeCmd;
    G4UIcommand* AddTreeFilesCmd;
    G4UIcommand* AddFileToTreeCmd;
    G4UIcommand* AddTreeToFileCmd;
    
    //security save
    G4UIcmdWithABool* SetSecuritySaveCmd;
    G4UIcommand* SetSecurityFileNameCmd;
    G4UIcmdWithAnInteger* SetSecurityNbOfEventsCmd;
    
    //scaling command
    G4UIcmdWithADoubleAndUnit* SetPrimaryIncidentFluxCmd;
    
    //verbose
    
    G4UIcmdWithAnInteger* SetBugVerboseCmd; 
     	    
    
    
    
       
  
    
    
    
 };

#endif

