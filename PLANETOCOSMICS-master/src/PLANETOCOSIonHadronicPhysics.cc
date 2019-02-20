////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSIonHadronicPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSIonHadronicPhysics::PLANETOCOSIonHadronicPhysics (const G4String& name)
  :G4VPhysicsConstructor(name), mode(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSIonHadronicPhysics::~PLANETOCOSIonHadronicPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4IonConstructor.hh"

void PLANETOCOSIonHadronicPhysics::ConstructParticle ()
{
}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ProcessManager.hh"
#include "G4hLowEnergyIonisation.hh"

#include "G4BinaryCascade.hh"


#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonsSihverCrossSection.hh"
#include "G4GeneralSpaceNNCrossSection.hh"


#include "G4CascadeInterface.hh"
#include "G4IonInelasticProcess.hh"
#include "G4TheoFSGenerator.hh"
#include "G4PreCompoundModel.hh"
#include "G4VPartonStringModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4WilsonAbrasionModel.hh"

#include "G4ElementTable.hh"
#include "G4Element.hh"



void PLANETOCOSIonHadronicPhysics::ConstructProcess ()
{ G4ProcessManager * pManager = 0;

//   // Elastic Process
//   //------------------
//   theElasticModel = new G4LElastic();
//   theIonElasticProcess.RegisterMe(theElasticModel);
// //   theDElasticProcess.RegisterMe(theElasticModel);
// //   theTElasticProcess.RegisterMe(theElasticModel);
// //   theAElasticProcess.RegisterMe(theElasticModel);
//   theHe3ElasticProcess.RegisterMe(theElasticModel);


 //InelasticProcess
  //----------------
  fDeuteronLEInelasticModel = new G4LEDeuteronInelastic();
  fTritonLEInelasticModel = new G4LETritonInelastic();
  fAlphaLEInelasticModel = new G4LEAlphaInelastic();
  
  
  //Cross sections
  //--------------
  
  G4IonsSihverCrossSection *  theSihverCS = new G4IonsSihverCrossSection;
  G4GeneralSpaceNNCrossSection * theGeneralSpaceNNCS = new  G4GeneralSpaceNNCrossSection;
  
  //Cross sections for alpha , deuterium and tritium
  //-------------------------------------------------
  
  theDeuteronInelasticProcess.AddDataSet(theSihverCS);
  theDeuteronInelasticProcess.AddDataSet(theGeneralSpaceNNCS); 
  
  theAlphaInelasticProcess.AddDataSet(theSihverCS);
  theAlphaInelasticProcess.AddDataSet(theGeneralSpaceNNCS); 
  
  theTritonInelasticProcess.AddDataSet(theSihverCS);
  theTritonInelasticProcess.AddDataSet(theGeneralSpaceNNCS); 
  
  //deuteron
  theDeuteronInelasticProcess.RegisterMe(fDeuteronLEInelasticModel);
  fDeuteronLEInelasticModel->SetMinEnergy(19.999*GeV);	
  //triton
  theTritonInelasticProcess.RegisterMe(fTritonLEInelasticModel);
  fTritonLEInelasticModel->SetMinEnergy(19.999*GeV);	
  //alpha
  theAlphaInelasticProcess.RegisterMe(fAlphaLEInelasticModel);
  fAlphaLEInelasticModel->SetMinEnergy(19.999*GeV);
  
  
  
  
  if (mode =="BIC"  || mode == "ABRASION"){ 
 
  	//Binary cascade
  	
// 	// this will be the model class for high energies not available now
//   	G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator; 
//   	// all models for treatment of thermal nucleus 
//   	G4Evaporation * theEvaporation = new G4Evaporation;
//   	G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
//   	G4StatMF * theMF = new G4StatMF;
// 
//   	// Evaporation logic
//   	G4ExcitationHandler * theHandler = new G4ExcitationHandler;
//   	theHandler->SetEvaporation(theEvaporation);
//   	theHandler->SetFermiModel(theFermiBreakUp);
//   	theHandler->SetMultiFragmentation(theMF);
//   	theHandler->SetMaxAandZForFermiBreakUp(12, 6);
//   	theHandler->SetMinEForMultiFrag(5*MeV);
// 	
//   	// Pre equilibrium stage 
//   	G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);
//     
//   	// a no-cascade generator-precompound interaface
//   	G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
//   	theCascade->SetDeExcitation(thePreEquilib);  
// 	
//   	G4VPartonStringModel * theStringModel;
//   	theStringModel = new G4QGSModel<G4QGSParticipants>;
//   	theTheoModel->SetTransport(theCascade);
//   	theTheoModel->SetHighEnergyGenerator(theStringModel);
//   	theTheoModel->SetMinEnergy(19.999*GeV);  // 15 GeV may be the right limit
//   	theTheoModel->SetMaxEnergy(100*TeV);
	
	
	G4BinaryLightIonReaction * theLightIonBC= new G4BinaryLightIonReaction;
	G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
	
	//Cross section for generic ions and He3
	
	theIonInelasticProcess.AddDataSet(theSihverCS);
	theIonInelasticProcess.AddDataSet(theGeneralSpaceNNCS);

	theHe3InelasticProcess.AddDataSet(theSihverCS);
	theHe3InelasticProcess.AddDataSet(theGeneralSpaceNNCS);
	
	
	if (mode == "BIC"){
  
  		// light Ion BC
  		theLightIonBC->SetMinEnergy(0.*MeV);
  		theLightIonBC->SetMaxEnergy(20*GeV);
   
		//Deuteron
	        theDeuteronInelasticProcess.RegisterMe(theLightIonBC);
  		//theDeuteronInelasticProcess.RegisterMe(theTheoModel);
	
		//Triton
  		theTritonInelasticProcess.RegisterMe(theLightIonBC);
  		//theTritonInelasticProcess.RegisterMe(theTheoModel);
	
		//Alpha
		theAlphaInelasticProcess.RegisterMe(theLightIonBC);
  		//theAlphaInelasticProcess.RegisterMe(theTheoModel);
	
		//Generic ion
		//------------
		//theGenIonBC->SetMinEnergy(0*MeV);
        	//theGenIonBC->SetMaxEnergy(10*GeV);
        	theIonInelasticProcess.RegisterMe(theLightIonBC);
       		// theIonInelasticProcess.RegisterMe(theTheoModel);
		
		//He3
		theHe3InelasticProcess.RegisterMe(theLightIonBC);
		
		
		
	}
	else { //The implementation of the abrasion model is taken following  "Truscott P., Nuclear-Nuclear Interaction Models in Geant4: Software User Manual" 
		G4WilsonAbrasionModel* theAM = new G4WilsonAbrasionModel(true);
		theAM->SetMinEnergy(100.*MeV);
		theAM->SetMaxEnergy(10.1*GeV); // coirrespond to limit in G4WilsonAbrasionModel code
		theGenIonBC->SetMinEnergy(0*MeV);
        	theGenIonBC->SetMaxEnergy(100*MeV);
		
		
		G4ElementTable::iterator iter;
		G4ElementTable *elementTable =
			const_cast<G4ElementTable*> (G4Element::GetElementTable());
		for (iter = elementTable->begin(); iter !=elementTable->end(); ++iter){
			 theAM->ActivateFor(*iter);
		}	
		
		//Deuteron
		theDeuteronInelasticProcess.RegisterMe(theGenIonBC);
	        theDeuteronInelasticProcess.RegisterMe(theAM);
  		
		//Triton
		theTritonInelasticProcess.RegisterMe(theGenIonBC);
  		theTritonInelasticProcess.RegisterMe(theAM);
  		
		//Alpha
		theAlphaInelasticProcess.RegisterMe(theGenIonBC);
		theAlphaInelasticProcess.RegisterMe(theAM);

		//Generic ion
		//------------
        	theIonInelasticProcess.RegisterMe(theGenIonBC);
       		theIonInelasticProcess.RegisterMe(theAM);
		
		//He3
		theHe3InelasticProcess.RegisterMe(theGenIonBC);
		theHe3InelasticProcess.RegisterMe(theAM);
	}	
       
       
      
   }
 	
   //Register the different processes
   //-------------------------------
   
   // Generic Ion
   pManager = G4GenericIon::GenericIon()->GetProcessManager();
   pManager->AddDiscreteProcess(&theIonElasticProcess);
   if (mode =="BIC"|| mode=="ABRASION" ) pManager->AddDiscreteProcess(&theIonInelasticProcess);


  
   // Deuteron
   pManager = G4Deuteron::Deuteron()->GetProcessManager();
  // pManager->AddDiscreteProcess(&theDElasticProcess);
   pManager->AddDiscreteProcess(&theDeuteronInelasticProcess);
   
   // Triton
   pManager = G4Triton::Triton()->GetProcessManager();
  // pManager->AddDiscreteProcess(&theTElasticProcess);
   pManager->AddDiscreteProcess(&theTritonInelasticProcess);
   
   //Alpha
   pManager = G4Alpha::Alpha()->GetProcessManager();
  // pManager->AddDiscreteProcess(&theAElasticProcess);
   pManager->AddDiscreteProcess(&theAlphaInelasticProcess);
   
   
   
   //He3
   pManager = G4He3::He3()->GetProcessManager();
   pManager->AddDiscreteProcess(&theHe3ElasticProcess);
   if (mode =="BIC" || mode=="ABRASION")pManager->AddDiscreteProcess(&theHe3InelasticProcess);
   

   
}
////////////////////////////////////////////////////////////////////////////////
