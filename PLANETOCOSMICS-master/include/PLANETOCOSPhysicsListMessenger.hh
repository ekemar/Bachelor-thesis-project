#ifndef PLANETOCOSPhysicsListMessenger_h
#define PLANETOCOSPhysicsListMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "PLANETOCOSPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSPhysicsListMessenger: public G4UImessenger
{
public:
  
  PLANETOCOSPhysicsListMessenger (PLANETOCOSPhysicsList* );
  ~PLANETOCOSPhysicsListMessenger ();
  
  void SetNewValue (G4UIcommand*, G4String);
  
private:
  
  PLANETOCOSPhysicsList* pPhysicsList;
  
  G4UIdirectory*              PhysDir;
  G4UIdirectory*              CutInRangeDir;
  
  G4UIcmdWithAString*         SelectTypeOfEMPhysicsCmd;
  G4UIcmdWithAString*         SelectTypeOfHadronicPhysicsCmd;
  G4UIcmdWithAString*         SelectTypeOfIonHadronicPhysicsCmd;
  G4UIcmdWithABool*	      SetConsiderElectroNuclearPhysicsCmd;
  G4UIcmdWithABool*	      SetConsiderMuonNuclearPhysicsCmd;
  G4UIcmdWithABool*	      SetConsiderSyncPhysicsCmd;
  G4UIcmdWithoutParameter*    ShowTypeOfPhysicsCmd;
  G4UIcmdWithoutParameter*    BuildListCmd;
  G4UIcmdWithADoubleAndUnit*  SetCutInDepthForAllLayersCmd;
  G4UIcmdWithoutParameter*    DeleteAllRegionsCmd;
  G4UIcommand* 		      SetCutInDepthForASpecificLayerCmd;	  
  G4UIcommand*                SetCutInLengthForASpecificLayerCmd;
  G4UIcmdWithADoubleAndUnit*  SetMinEForMultiFragCmd;
  G4UIcommand*  	      SetMaxAandZForFermiBreakUpCmd;
  G4UIcommand*  	      SetSteppingAlgorithmParametersForMscCmd;
 
};
////////////////////////////////////////////////////////////////////////////////
#endif
