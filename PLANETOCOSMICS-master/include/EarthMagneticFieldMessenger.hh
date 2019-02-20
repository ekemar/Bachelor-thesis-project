#ifndef EarthMagneticFieldMessenger_h
#define EarthMagneticFieldMessenger_h 1

#include "globals.hh"
#include"PlanetMagneticFieldMessenger.hh"

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
class EarthMagneticField;


class EarthMagneticFieldMessenger : public PlanetMagneticFieldMessenger
{
public:
   EarthMagneticFieldMessenger(EarthMagneticField* aField);
   ~EarthMagneticFieldMessenger();
   
   
   void SetNewValue(G4UIcommand * command,G4String newValues);
   
   
   
   private:
   EarthMagneticField* theField;  
   
   
   G4UIcmdWithoutParameter* SetTiltedDipoleParameterFromIGRFCmd;
   G4UIcmdWithoutParameter* SetEccentricDipoleParameterFromIGRFCmd;
 
   
   // Magnetic Activity Command
     
   G4UIcmdWithAnInteger* SetIoptCmd;
   G4UIcmdWithADouble* SetPdynCmd;
   G4UIcmdWithADouble* SetTiltAngleCmd;
   G4UIcmdWithADoubleAndUnit* SetDstCmd;
   G4UIcmdWithADoubleAndUnit* SetImfyCmd;
   G4UIcmdWithADoubleAndUnit* SetImfzCmd;
   G4UIcmdWithADouble* SetG1Cmd;
   G4UIcmdWithADouble* SetG2Cmd;
   G4UIcmdWithADouble* SetW1Cmd;
   G4UIcmdWithADouble* SetW2Cmd;
   G4UIcmdWithADouble* SetW3Cmd;
   G4UIcmdWithADouble* SetW4Cmd;
   G4UIcmdWithADouble* SetW5Cmd;
   G4UIcmdWithADouble* SetW6Cmd;
   G4UIcmdWithAString* ReadTSY2001ParameterCmd;
   G4UIcmdWithoutParameter* PrintTSY2001ParameterCmd;
   G4UIcmdWithAString* ReadTSY2004ParameterCmd;
   G4UIcmdWithoutParameter* PrintTSY2004ParameterCmd;   
   
};   
#endif
