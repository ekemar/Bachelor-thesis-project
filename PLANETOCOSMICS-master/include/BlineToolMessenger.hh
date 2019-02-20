//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolMessenger.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines interactive commands that allows to pilot 
//   the BlineTool object                          
// 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef BlineToolMessenger_h
#define BlineToolMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class BlineTool;
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

class BlineToolMessenger : public G4UImessenger
{public:
   BlineToolMessenger(BlineTool* aBlineTool);
  ~BlineToolMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
   BlineTool* theBlineTool;
   G4UIdirectory*             BlineToolDir;
   
 
   //  commands
   
   G4UIcmdWithAnInteger* BlineCmd;
   G4UIcmdWithADoubleAndUnit* SetMaxTrackingStepCmd;
   G4UIcmdWith3Vector* SetDrawColourCmd;
   G4UIcmdWithABool*   SetDrawBlineCmd;    
   G4UIcmdWithABool*  SetDrawPointsCmd;
   G4UIcmdWithADouble* SetPointSizeCmd;
   G4UIcmdWithoutParameter* DrawCmd;
   G4UIcmdWithoutParameter* ResetCmd;
  
   
    

};












#endif
