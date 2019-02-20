#ifndef PLANETOCOSPrimaryMessenger_h
#define PLANETOCOSPrimaryMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PLANETOCOSPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class PLANETOCOSPrimaryMessenger: public G4UImessenger
{
  public:
    PLANETOCOSPrimaryMessenger(PLANETOCOSPrimaryGeneratorAction * aGenerator);
    ~PLANETOCOSPrimaryMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    PLANETOCOSPrimaryGeneratorAction* myGenerator;
    G4UIdirectory*             myStartGeoDir;
    G4UIdirectory*             myRigidityFilterDir;
    G4UIdirectory*             myBlineTracingDir;
    
  /*  G4UIcmdWithADouble* SetZenithCmd;
    G4UIcmdWithADouble* SetAzimuthCmd;*/
    G4UIcmdWithADoubleAndUnit* SetRigidityCmd;
   /* G4UIcmdWith3Vector* SetGeoPositionCmd;
    G4UIcmdWithoutParameter* ComputeStartPositionCmd; */
    G4UIcommand* SetPositionVectorCmd;
    G4UIcommand* SetPositionCmd;
    G4UIcommand* SetPositionOnDipoleMagneticShellCmd;
    G4UIcommand* SetDirectionCmd;
    G4UIcommand* SetDirectionVectorCmd;
    G4UIcommand* SetPositionAndDirectionVectorCmd;
    G4UIcommand* SetPositionAndDirectionCmd;
    
    G4UIcmdWith3VectorAndUnit* SetPositionVectorFlatCaseCmd;
    G4UIcmdWith3Vector* SetDirectionVectorFlatCaseCmd;
    G4UIcommand* SetDirectionFlatCaseCmd;
    
    
    
    
    G4UIcommand* AddValuesToRigidityVectorCmd;
    G4UIcmdWithoutParameter* SetDefaultRigidityVectorCmd;
    G4UIcmdWithoutParameter* ResetRigidityVectorCmd;
    G4UIcmdWithAString* ComputeRigidityFilterCmd;
    G4UIcmdWithAString* DefineBlineFootPointsCmd;
    G4UIcmdWithoutParameter* ComputeBlinesCmd;
    G4UIcmdWithoutParameter* TraceBlinesAndParticlesCmd;
    
    G4UIcmdWithAnInteger*    SetVerbosityCmd;
    G4UIcmdWithoutParameter* PrintBfieldAtPrimaryCmd;
    
    G4UIcommand* SelectModulatedGalacticFluxCmd;
    G4UIcommand* SelectSolMaxGalacticFluxCmd;
    G4UIcommand* SelectSolMinGalacticFluxCmd;
    G4UIcommand* SelectMeanGalacticFluxCmd; 
    G4UIcmdWithAString* ReadPrimaryFluxCmd;
    
    G4UIcommand* RandomIsotropicDistributionCmd;
    
    //Cutoff command
    
    G4UIcmdWithABool* SetConsiderCutoffCmd;
    G4UIcmdWithADoubleAndUnit* SetCutoffRigidityCmd;
    G4UIcmdWithAString* ReadCutoffVsDirectionCmd;

#ifdef USE_ROOT_SOURCE
    G4UIcmdWithABool* SetUseRootSourceCmd;
    G4UIcmdWithAString* SetVar1NameForRootSourceCmd;
    G4UIcmdWithAString* SetVar2NameForRootSourceCmd;
    G4UIcommand* ReadFirstHistoForRootSourceCmd;
    G4UIcommand* ReadSecondHistoForRootSourceCmd;
    
#endif
    
 };

#endif

