#ifndef EarthAtmosphereMessenger_h
#define EarthAtmosphereMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class EarthAtmosphereModel;
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


class EarthAtmosphereMessenger: public G4UImessenger
{
  public:
    EarthAtmosphereMessenger(EarthAtmosphereModel * myMod);
    ~EarthAtmosphereMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    EarthAtmosphereModel* myModel;
    
    
    //atmosphere model  command specific for the Earth
    G4UIcommand*  SetAtmosphereReferenceDateCmd;
    G4UIcmdWithADouble*  SetApCmd;
    G4UIcmdWithADouble*  SetF107Cmd;
    G4UIcmdWithADouble*  SetF107ACmd;
    G4UIcommand*   SetPositionCmd;
   

};

#endif

