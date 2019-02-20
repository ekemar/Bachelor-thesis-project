#ifndef PLANETOCOSSteppingActionMessenger_h
#define PLANETOCOSSteppingActionMessenger_h 1

class PLANETOCOSSteppingAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class PLANETOCOSSteppingActionMessenger: public G4UImessenger
{
  public:
    PLANETOCOSSteppingActionMessenger(PLANETOCOSSteppingAction* msa);
    ~PLANETOCOSSteppingActionMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    PLANETOCOSSteppingAction * myAction;
    G4UIdirectory*             myStackingDir;
    G4UIdirectory*             myStopBoundaryDir;
  
  private: //commands
    
    //Commands for stopping particle below user defined energy
    G4UIcommand* AddUntrackedParticleCmd;
    G4UIcmdWithAString* RemoveUntrackedParticleCmd;
     
    //Commands for controlling stopping at user defined boundary
    G4UIcmdWithABool* SetStopUpFluxAtSelectedBoundaryCmd;
    G4UIcmdWithABool* SetStopDownFluxAtSelectedBoundaryCmd;
    G4UIcmdWithAString* SetNameOfStopBoundaryForUpFluxCmd;
    G4UIcmdWithAString* SetNameOfStopBoundaryForDownFluxCmd;
    
    //Commands for controlling stopping at magnetopause
    G4UIcmdWithABool* SetStopAtMagnetopauseCmd;
    G4UIcmdWithADouble* SetMagnetopauseOutFactorCmd;
    
    
    
    //Commands for controlling stopping at user defined altitude
    G4UIcmdWithADoubleAndUnit* SetStopAltitudeForUpwardCmd;
    G4UIcmdWithADoubleAndUnit* SetStopAltitudeForDownwardCmd;
    G4UIcmdWithABool*  SetStopUpwardPrimaryCmd;
    G4UIcmdWithABool*  SetStopDownwardPrimaryCmd;
    
    //Commands for stopping particle after some turns around the planet 
    G4UIcmdWithADouble* SetMaxNumberOfTurnAroundThePlanetCmd;
    
    
};

#endif


