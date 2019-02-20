#ifndef PLANETOCOSGeometryMessenger_h
#define PLANETOCOSGeometryMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PLANETOCOSGeometryConstruction;
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


class PLANETOCOSGeometryMessenger: public G4UImessenger
{
  public:
    PLANETOCOSGeometryMessenger(PLANETOCOSGeometryConstruction * myDet);
    ~PLANETOCOSGeometryMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    PLANETOCOSGeometryConstruction* myGeometry;
    G4UIdirectory*             atmocosmicsDir;
    G4UIdirectory*             myGeometryDir;
    G4UIdirectory* 	       myUserlimitDir;
    
    //geometry command
    
    G4UIcmdWithADoubleAndUnit* AtmosphereMaxStepLengthCmd;
    G4UIcmdWithADoubleAndUnit* MagnetosphereMaxStepLengthCmd;
    G4UIcmdWithADoubleAndUnit* MagnetosphereMaxTrackLengthCmd;
    G4UIcmdWithADoubleAndUnit* MagnetosphereMaxTrackDurationCmd;
    
    G4UIcmdWithADoubleAndUnit* SetAtmosphereHmaxCmd;
    G4UIcmdWithADoubleAndUnit* SetAtmosphereHminCmd;
    G4UIcmdWithADoubleAndUnit* SetMagnetosphereHCmd;
    G4UIcmdWithADoubleAndUnit* SetPlanetHCmd;
    G4UIcmdWithAnInteger* SetVerbosityCmd;
    
    
    G4UIcmdWithADoubleAndUnit* SetMaxThicknessCmd;
    G4UIcmdWithADoubleAndUnit* SetMinThicknessCmd;
    G4UIcmdWithADoubleAndUnit* SetHalfLengthCmd;
    G4UIcmdWithADouble* SetDepthPercentCmd;
    G4UIcmdWithAString* SetGeometryTypeCmd;
    G4UIcmdWithoutParameter* UpdateGeometryCmd;
    G4UIcmdWithoutParameter* RemoveAllDetectorsCmd;
    G4UIcmdWithADoubleAndUnit* AddDetectorAtAltitudeCmd;
    G4UIcmdWithADoubleAndUnit* AddDetectorAtDepthCmd;
    
    
    
    //atmosphere model command
    
    G4UIcmdWithAString* SetAtmosphereModelCmd;
    G4UIcmdWithAString* ReadAtmosphereModelCmd;
    G4UIcmdWithABool* SetConsiderAtmosphereCmd;
    G4UIcmdWithABool* SetDetectInMiddleOfAtmosphereModeCmd;
    G4UIcmdWithABool* SetDetectionBelowSoilLayersCmd;
    

};

#endif

