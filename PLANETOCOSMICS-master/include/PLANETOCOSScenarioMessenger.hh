#ifndef PLANETOCOSScenarioMessenger_h
#define PLANETOCOSScenarioMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PLANETOCOSApplicationScenario;
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

class PLANETOCOSScenarioMessenger : public G4UImessenger
{public:
   PLANETOCOSScenarioMessenger(PLANETOCOSApplicationScenario* ApplicationScenario);
  ~PLANETOCOSScenarioMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
   PLANETOCOSApplicationScenario* theApplicationScenario;
   G4UIdirectory*             ScenarioDir;
   
 
  
   //General run command
   
   G4UIcmdWithoutParameter* SetRandomSeedCmd;
   G4UIcmdWithABool* SetRandomSeedAtRunStartCmd;
   G4UIcmdWithADoubleAndUnit* SetMaxDurationForARunCmd;
   
   
   
   
    

};












#endif
