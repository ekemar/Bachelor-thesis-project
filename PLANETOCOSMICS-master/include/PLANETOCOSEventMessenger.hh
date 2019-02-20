#ifndef PLANETOCOSEventMessenger_h
#define PLANETOCOSEventMessenger_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              PLANETOCOSGeometryConstruction.hh
//
// Version:		VERSION_NUMBER
// Date:		LAST_DATE
// Author:		L Desorgher
// Organisation:	ORGANISATION_NAME
// Project:		PROJECT_NAME
// Customer:		CUSTOMER_NAME
// Contract:		CONTRACT_NAME
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This class defines the geometry messenger for MAGNETOCOSMICS. 
// It defines UI commond that allows to control interactively 
// the PLANETOCOSGeometryContruction object.
//  
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
// 

#include "globals.hh"
#include "G4UImessenger.hh"

class PLANETOCOSEventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class PLANETOCOSEventMessenger: public G4UImessenger
{
  public:
    PLANETOCOSEventMessenger(PLANETOCOSEventAction * myAct);
    ~PLANETOCOSEventMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    PLANETOCOSEventAction* myAction;
    G4UIdirectory*     myDrawDir;
   
   /* G4UIcmdWithADouble* SetDrawLineWidthCmd;
    G4UIcmdWithAnInteger* SetDrawLineStyleCmd;*/
    
    G4UIcmdWith3Vector* SetDrawColourForBlineCmd;
    G4UIcmdWithABool*   SetDrawTrajectoryCmd;    
    G4UIcmdWithAString* SetDrawingCoordinateSystemCmd;
    G4UIcmdWithABool*  SetDrawPointsCmd;
    G4UIcmdWithADouble* SetPointSizeCmd;
    G4UIcmdWithoutParameter* DrawCmd;
    G4UIcmdWithoutParameter* ResetCmd;
    
    G4UIcmdWithADoubleAndUnit* TraceMagnetopauseLineCmd; 
    G4UIcommand* AddParticleToBeDrawnCmd;
    G4UIcmdWithAString* RemoveParticleToBeDrawnCmd;
    
 };

#endif

