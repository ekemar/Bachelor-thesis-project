
#ifndef MagneticShieldingToolMessenger_h
#define MagneticShieldingToolMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MagneticShieldingTool;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

class MagneticShieldingToolMessenger : public G4UImessenger
{public:
   MagneticShieldingToolMessenger(MagneticShieldingTool* aMagneticShieldingTool);
  ~MagneticShieldingToolMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
   MagneticShieldingTool* theMagneticShieldingTool;
   G4UIdirectory*             MagneticShieldingToolDir;
   
 
   // Scenario command
   G4UIcmdWithoutParameter* BlineCmd;
   G4UIcmdWithAnInteger* ParticleTrajectoryCmd;
   G4UIcmdWithAnInteger* ReverseParticleTrajectoryCmd;
   G4UIcmdWithADoubleAndUnit* SetStopAltitudeCmd;
   G4UIcmdWithAString* ComputeRigidityFilterCmd;
   G4UIcommand* ComputeDirectionFilterCmd;
   G4UIcommand* RCutoffVsPositionCmd;
   G4UIcommand* RCutoffVsPositionOnLShellCmd;
   G4UIcommand* RCutoffVsSpenvisPositionGridCmd;
   G4UIcommand* RCutoffVsSpenvisTrajectoryCmd;
   G4UIcommand* RCutoffVsDirectionCmd;
   G4UIcommand* RCutoffVsTimeCmd;
   G4UIcmdWithABool* SetAutoDetectionOfPenumbra;
   G4UIcmdWithoutParameter* InitialiseCmd; 
   
   //SpenvisCSV command
   G4UIcmdWithABool* SetRegisterResultsInSpenvisCSVFileCmd;
   G4UIcmdWithAString* SetSpenvisCSVFileNameCmd;
 
  
   
    

};












#endif
