////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSMuonPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSMuonPhysics::PLANETOCOSMuonPhysics(const G4String& name)
                   :  G4VPhysicsConstructor(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSMuonPhysics::~PLANETOCOSMuonPhysics()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

void PLANETOCOSMuonPhysics::ConstructParticle()
{
  // Mu
 
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Tau
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();

}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ProcessManager.hh"

void PLANETOCOSMuonPhysics::ConstructProcess()
{
  G4ProcessManager * pmanager = 0;

  // Muon Plus Physics
  pmanager = G4MuonPlus::MuonPlus()->GetProcessManager();
   // add processes 
  pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
  pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
  pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
  pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);
 

  // Muon Minus Physics
  pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
  pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
  pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
  pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);
  pmanager->AddRestProcess(new G4MuonMinusCaptureAtRest);
  
}
////////////////////////////////////////////////////////////////////////////////
