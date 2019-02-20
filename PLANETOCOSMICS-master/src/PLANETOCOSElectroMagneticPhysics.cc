////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSElectroMagneticPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

      			
//Extra EM physics

#include "G4ElectroNuclearBuilder.hh"
#include "G4MuNuclearInteraction.hh"
#include "G4SynchrotronRadiation.hh"
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSElectroMagneticPhysics::PLANETOCOSElectroMagneticPhysics (const G4String& name,
							      const bool consider_emnuc= false,
							      const bool consider_munuc= false,
							       const bool consider_sync= false)
  : G4VPhysicsConstructor(name), mode(name), ConsiderEMNucPhysics(consider_emnuc),ConsiderMuonNucPhysics(consider_munuc),ConsiderSyncPhysics(consider_sync)
{
   theEMNucPhysics =0;
   theMuMinusNucInt =0;
   theMuPlusNucInt =0;
   theElectronSyncRad =0;
   thePositronSyncRad=0;
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSElectroMagneticPhysics::~PLANETOCOSElectroMagneticPhysics ()
{ if (theEMNucPhysics) delete theEMNucPhysics;
  if ( theMuMinusNucInt) delete  theMuMinusNucInt;
  if ( theMuPlusNucInt) delete  theMuPlusNucInt;
  if ( theElectronSyncRad) delete  theElectronSyncRad;
  if ( thePositronSyncRad) delete  thePositronSyncRad;	 

}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

//gamma, e-, positron
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4GenericIon.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4He3.hh"

//hadrons
#include "G4Proton.hh"
#include "G4AntiProton.hh"




#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiOmegaMinus.hh"






////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSElectroMagneticPhysics::ConstructParticle ()
{ // gamma
  G4Gamma::GammaDefinition();
 
  // electron positron 
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  
  //light ions
   G4IonConstructor pConstructor;
   pConstructor.ConstructParticle();  
  
  //hadrons
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

  G4SigmaPlus::SigmaPlusDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4XiMinus::XiMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();
  
  
  
  // for test of stepping bounddary algorithm for multiple scattering
   SteppingAlgorithmMsc=true;
   facrange = -1.;
}
////////////////////////////////////////////////////////////////////////////////

#include "G4ProcessManager.hh"




// gamma e. positron

//standard
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
//low energy
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 
#include "G4ProductionCutsTable.hh"


//hadrons
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4hLowEnergyIonisation.hh" 


//Muon
//-------
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"


void PLANETOCOSElectroMagneticPhysics::ConstructProcess()
{  
 
  

  if (mode == "LowEnergyEM"){
  	G4LowEnergyPhotoElectric* lowePhot = new G4LowEnergyPhotoElectric();
  	G4LowEnergyIonisation* loweIon  = new G4LowEnergyIonisation();
  	G4LowEnergyBremsstrahlung* loweBrem = new G4LowEnergyBremsstrahlung();	
  	G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
	//aMultipleScattering->MscStepLimitation(SteppingAlgorithmMsc, facrange);
  	
	//gamma
	G4ProcessManager * pmanager = G4Gamma::Gamma()->GetProcessManager();
	pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());
        pmanager->AddDiscreteProcess(lowePhot);
        pmanager->AddDiscreteProcess(new G4LowEnergyCompton());
        pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
	
	//electron
	pmanager = G4Electron::Electron()->GetProcessManager();
	pmanager->AddProcess(aMultipleScattering,    -1, 1,1);
        pmanager->AddProcess(loweIon,                -1, 2,2);
        pmanager->AddProcess(loweBrem,               -1,-1,3);
	
	//positron
	pmanager = G4Positron::Positron()->GetProcessManager();
  	pmanager->AddProcess(aMultipleScattering, -1, 1,1);
      	pmanager->AddProcess(new G4eIonisation(),        -1, 2,2);
      	pmanager->AddProcess(new G4eBremsstrahlung(),    -1,-1,3);
      	pmanager->AddProcess(new G4eplusAnnihilation(),   0,-1,4);
	
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV,1e5*MeV);
        //
        //fluorescence apply specific cut for flourescence from photons, electrons
        //and bremsstrahlung photons:
        G4double cut = 250*eV;
        lowePhot->SetCutForLowEnSecPhotons(cut);
        loweIon->SetCutForLowEnSecPhotons(cut);
        loweBrem->SetCutForLowEnSecPhotons(cut);
	
	
	//Muon
	//----------
	pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();
	pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
      	pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
      	pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
      	pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);
	
	pmanager = G4MuonPlus::MuonPlus()->GetProcessManager();
	pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
      	pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
      	pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
      	pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);
	
	
	
	
	
	//hadron and nuclei
	//-----------------
	std::vector<G4ParticleDefinition *> partVec;
 
  	partVec.push_back(G4GenericIon::GenericIon());
  	partVec.push_back(G4Deuteron::Deuteron());
  	partVec.push_back(G4Triton::Triton());
  	partVec.push_back(G4Alpha::Alpha());
 	partVec.push_back(G4He3::He3());
        partVec.push_back(G4PionPlus::PionPlus());
  	partVec.push_back(G4PionMinus::PionMinus());
 	partVec.push_back(G4KaonPlus::KaonPlus());
  	partVec.push_back(G4KaonMinus::KaonMinus());
 	partVec.push_back(G4Proton::Proton());
  	partVec.push_back(G4AntiProton::AntiProton());
  	partVec.push_back(G4SigmaMinus::SigmaMinus());
  	partVec.push_back(G4AntiSigmaMinus::AntiSigmaMinus());
  	partVec.push_back(G4SigmaPlus::SigmaPlus());
  	partVec.push_back(G4AntiSigmaPlus::AntiSigmaPlus());
  	partVec.push_back(G4XiMinus::XiMinus());
  	partVec.push_back(G4AntiXiMinus::AntiXiMinus());
	partVec.push_back(G4TauMinus::TauMinus());
	partVec.push_back(G4TauPlus::TauPlus());
	partVec.push_back(G4OmegaMinus::OmegaMinus());
 	partVec.push_back(G4AntiOmegaMinus::AntiOmegaMinus());
  
  	for (unsigned int i=0;i< partVec.size();i++){
  		pmanager=partVec[i]->GetProcessManager();
  		pmanager->AddProcess( new G4hLowEnergyIonisation(),ordInActive,2, 2);
		G4MultipleScattering* theMulScattering = new G4MultipleScattering();
		pmanager->AddProcess(theMulScattering);
		pmanager->SetProcessOrdering(theMulScattering, idxAlongStep, 1);
  		pmanager->SetProcessOrdering(theMulScattering, idxPostStep, 1);	
 	}
  }
  else {//standard mode// Add standard EM Processes
	
	
	
	
	theParticleIterator->reset();
  	while( (*theParticleIterator)() ){
    		G4ParticleDefinition* particle = theParticleIterator->value();
    		G4ProcessManager* pmanager = particle->GetProcessManager();
    		G4String particleName = particle->GetParticleName();
		//G4MultipleScattering* theMsc = new G4MultipleScattering();
		//G4cout<<"Stepping algorithm "<<SteppingAlgorithmMsc<<'\t'<<facrange<<std::endl;
		//theMsc->MscStepLimitation(SteppingAlgorithmMsc, facrange);
    		if (particleName == "gamma") {

      			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      			pmanager->AddDiscreteProcess(new G4ComptonScattering);
      			pmanager->AddDiscreteProcess(new G4GammaConversion);

    		}
			 
		else if (particleName == "e-") {
			pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      			pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      			pmanager->AddProcess(new G4eBremsstrahlung(),  -1, 3, 3);

    		} 
		
	

		else if (particleName == "e+") {
     			pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      			pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      			pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);
      			pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);
		
        	}
	
		
		
		else if (particleName == "mu+" ||
               				particleName == "mu-"    ) {

      			pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
      			pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
      			pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
      			pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);

    		}
		 
    		else if (particleName == "alpha" ||
               		particleName == "He3" ||
               		particleName == "GenericIon") {
	       		pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
               		pmanager->AddProcess(new G4ionIonisation,      -1, 2, 2);
	       
	       
	       
		} 
		else if (particleName == "anti_omega-" ||
               		particleName == "anti_proton" ||
               		particleName == "anti_sigma+" ||
               		particleName == "anti_sigma-" ||
               		particleName == "anti_xi-" ||
               		particleName == "deuteron" ||
               		particleName == "kaon+" ||
               		particleName == "kaon-" ||
               		particleName == "omega-" ||
               		particleName == "pi+" ||
               		particleName == "pi-" ||
              		particleName == "proton" ||
               		particleName == "sigma+" ||
               		particleName == "sigma-" ||
               		particleName == "tau+" ||
              		particleName == "tau-" ||
               		particleName == "triton" ||
               		particleName == "xi-" ) {
				pmanager->AddProcess(new G4MultipleScattering,-1, 1, 1);
      				pmanager->AddProcess(new G4hIonisation,       -1, 2, 2);
    			}
			
		
		
		
	}
		
  }
 /* G4EmProcessOptions opt;
  opt.SetVerbose(verbose);*/
 
       
  	
	
  //}
  
  
  
   theEMNucPhysics =0;
   theMuMinusNucInt =0;
   theMuPlusNucInt =0;
   theElectronSyncRad =0;
   thePositronSyncRad=0;
  
  
  
  if(ConsiderSyncPhysics) {
    G4ProcessManager* pManager = G4Electron::Electron()->GetProcessManager();
    theElectronSyncRad = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(theElectronSyncRad);

    pManager = G4Positron::Positron()->GetProcessManager();
    thePositronSyncRad = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(thePositronSyncRad);
  }
  
  if (ConsiderEMNucPhysics) {
    G4ElectroNuclearBuilder* theEMNucPhysics = new G4ElectroNuclearBuilder();
    theEMNucPhysics->Build();
  }
  if(ConsiderMuonNucPhysics) {
    G4ProcessManager* pManager  = G4MuonPlus::MuonPlus()->GetProcessManager();
    theMuPlusNucInt = new G4MuNuclearInteraction("muNucl");
    pManager->AddDiscreteProcess(theMuPlusNucInt);

    pManager  = G4MuonMinus::MuonMinus()->GetProcessManager();
    theMuMinusNucInt = new G4MuNuclearInteraction("muNucl");
    pManager->AddDiscreteProcess(theMuMinusNucInt);
  }
 
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
