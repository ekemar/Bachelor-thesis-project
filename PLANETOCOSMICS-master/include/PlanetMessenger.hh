#ifndef PlanetMessenger_h
#define PlanetMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PlanetManager;
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


class PlanetMessenger: public G4UImessenger
{
  public:
    PlanetMessenger(PlanetManager * theManager);
    ~PlanetMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    PlanetManager*  pPlanetManager;
    G4UIdirectory*  planetDir;
    G4UIdirectory*  testDir;
    
    G4UIcmdWithAString* SelectPlanetCmd;
    
    //Space coordinate 
    G4UIcommand* TestMarinerOrbitCmd;
    G4UIcommand* TestMarsOdysseyCmd;
    G4UIcommand* TestImageOrbitCmd;
    G4UIcommand* TestWindOrbitCmd;

//#ifdef USE_SPICE
    G4UIcmdWithABool* SetUseSpiceCmd;
//#endif   
       
  
    
    
    
 };

#endif

