#ifndef PlanetSoilMessenger_h
#define PlanetSoilMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PlanetSoil;
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

class PlanetSoilMessenger : public G4UImessenger
{public:
   PlanetSoilMessenger(PlanetSoil* aField);
  ~PlanetSoilMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
     
     PlanetSoil* theSoil;
    
     G4UIdirectory* SoilDir;
     
     G4UIcommand*  AddMonoElementLayerCmd;
     G4UIcommand*  AddLayerCmd;
     G4UIcommand*  AddElementToLayerCmd;
     G4UIcmdWithoutParameter* ResetLayersCmd;
   
     
     
  



    
};












#endif
