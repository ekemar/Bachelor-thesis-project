#ifndef MercuryMagneticFieldMessenger_h
#define MercuryMagneticFieldMessenger_h 1

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
class MercuryMagneticField;


class MercuryMagneticFieldMessenger : public PlanetMagneticFieldMessenger
{
public:
   MercuryMagneticFieldMessenger(MercuryMagneticField* aField);
   ~MercuryMagneticFieldMessenger();
   
   
   void SetNewValue(G4UIcommand * command,G4String newValues);
   
   
   
   private:
   MercuryMagneticField* theField;
   
};   
#endif