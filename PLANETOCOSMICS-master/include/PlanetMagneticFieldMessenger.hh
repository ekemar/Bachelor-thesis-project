#ifndef PlanetMagneticFieldMessenger_h
#define PlanetMagneticFieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PlanetMagneticField;
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

class PlanetMagneticFieldMessenger : public G4UImessenger
{public:
   PlanetMagneticFieldMessenger(PlanetMagneticField* aField);
  ~PlanetMagneticFieldMessenger();
  
   virtual void SetNewValue(G4UIcommand * command,G4String newValues)=0;
 
 protected:
     
     G4UIdirectory*             mainDir;
     G4UIdirectory*             IntegrationDir;
     G4UIdirectory*             MagnetoDir;
     
     PlanetMagneticField* theMotherField;

   //integration command
  
     G4UIcmdWithADouble* SetEpsilonCmd;
     G4UIcmdWithADoubleAndUnit* SetG4DeltaChordCmd; 
     G4UIcmdWithADoubleAndUnit* SetDeltaIntersectionCmd;
     G4UIcmdWithoutParameter*  ResetIntegrationParametersCmd;
     G4UIcmdWithAString* SetStepperCmd;
     
   // magnetic field model command
     
     G4UIcmdWithADoubleAndUnit* SetTimeOfBCmd;
     G4UIcommand* SetStartDateCmd;
     G4UIcmdWithAnInteger* Setnmax_GaussCmd;
     G4UIcmdWithAString* SetInternalFieldCmd;
     G4UIcmdWithAString* SetExternalFieldCmd;
     G4UIcmdWithAString* SetMagnetopauseModelCmd;
     G4UIcmdWithADoubleAndUnit* SetDipoleB0Cmd;
     G4UIcmdWith3VectorAndUnit* SetDipoleShiftCmd; 
     G4UIcommand* SetDipoleAxisCmd; 
     G4UIcmdWithADoubleAndUnit* SetRadiusMagnetosphereCmd;
     G4UIcmdWithABool* SetConsiderDipoleShiftCmd; 
     G4UIcommand* ComputeBfieldAtDifferentPositions;
     G4UIcommand* ComputeBfieldAtDifferentPositions1;
     G4UIcommand* TraceBlineAtDifferentPositions;
     G4UIcommand* TraceBlineAtDifferentPositions1;
     G4UIcmdWith3VectorAndUnit* SetHomogeneousFieldCmd;
     G4UIcommand* SetHomogeneousFieldFromSelectedModelCmd;
     G4UIcommand* SetWorldCenterPositionCmd;
     G4UIcmdWithoutParameter* SwitchOnFieldCmd; 
     G4UIcmdWithoutParameter* SwitchOffFieldCmd;
     G4UIcommand* ComputeInterpolMatrixForFlatGeometryCmd;
     G4UIcmdWithAString*  ReadInterpolMatrixForFlatGeometryCmd;
protected :     
     
     bool SetMotherNewValue(G4UIcommand * command,G4String newValues);

    
};












#endif
