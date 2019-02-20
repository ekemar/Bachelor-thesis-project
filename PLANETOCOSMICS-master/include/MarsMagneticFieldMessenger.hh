#ifndef MarsMagneticFieldMessenger_h
#define MarsMagneticFieldMessenger_h 1

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
class MarsMagneticField;


class MarsMagneticFieldMessenger : public PlanetMagneticFieldMessenger
{
public:
   MarsMagneticFieldMessenger(MarsMagneticField* aField);
   ~MarsMagneticFieldMessenger();
   
   
   void SetNewValue(G4UIcommand * command,G4String newValues);
   
   
   
   private:
   MarsMagneticField* theField;
   
};   
#endif
