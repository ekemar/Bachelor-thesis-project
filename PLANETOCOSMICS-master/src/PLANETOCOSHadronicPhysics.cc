////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSHadronicPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

//Binary cascade
#include "MYBinaryNeutronBuilder.hh"
#include "MYBinaryProtonBuilder.hh"


////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSHadronicPhysics::PLANETOCOSHadronicPhysics (const G4String& name)
  : G4VPhysicsConstructor(name), HadronicPhysicsType(name)
{

 theBinaryProton = 0;
 theBinaryNeutron = 0;
 SetMaxAandZForFermiBreakUp(13,7);
 //SetMinEForMultiFrag(5.*MeV); 

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSHadronicPhysics::~PLANETOCOSHadronicPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSHadronicPhysics::ConstructParticle ()
{
  
  G4cout<<"Hello"<<std::endl;
  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

}
////////////////////////////////////////////////////////////////////////////////
//
//builders

//LHEP
#include "G4LHEPPiKBuilder.hh"
#include "G4LHEPProtonBuilder.hh"
#include "G4LHEPNeutronBuilder.hh"


//LEP
#include "G4LEPPiKBuilder.hh"
#include "G4LEPProtonBuilder.hh"
#include "G4LEPNeutronBuilder.hh"

//HP
#include "G4NeutronHPBuilder.hh"

//Precompound model
#include "G4PrecoNeutronBuilder.hh"
#include "G4PrecoProtonBuilder.hh"

//QGSP
#include "G4QGSPNeutronBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4QGSPPiKBuilder.hh"

//QGSC
#include "G4QGSCNeutronBuilder.hh"
#include "G4QGSCProtonBuilder.hh"
#include "G4QGSCPiKBuilder.hh"

//FTFP
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4FTFPPiKBuilder.hh"

//FTFC
#include "G4FTFCNeutronBuilder.hh"
#include "G4FTFCProtonBuilder.hh"
#include "G4FTFCPiKBuilder.hh"


//Bertini
#include "G4BertiniNeutronBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

//Elastic physics

#include "G4HadronElasticPhysics.hh"
void PLANETOCOSHadronicPhysics::ConstructProcess()
{ 
  NeutronHPElasticFlag =true;
  
  
  
  
  G4cout<<"You have selected the "+HadronicPhysicsType+" hadronic physics model"<<std::endl;
  if (HadronicPhysicsType == "LHEP")			ConstructLHEPBasedPhysics(false,false,false); 		
  else if (HadronicPhysicsType == "LHEP_HP") 		ConstructLHEPBasedPhysics(false,false,true);
  else if (HadronicPhysicsType == "LHEP_BIC")		ConstructLHEPBasedPhysics(true,false,false);
  else if (HadronicPhysicsType == "LHEP_BIC_HP")	ConstructLHEPBasedPhysics(true,false,true);
  else if (HadronicPhysicsType == "LHEP_BERT")		ConstructLHEPBasedPhysics(false,true,false);
  else if (HadronicPhysicsType == "LHEP_BERT_HP")	ConstructLHEPBasedPhysics(false,true,true);
 
  if (HadronicPhysicsType == "QGSP")			ConstructPartonStringBasedPhysics(false,false,false,false,false); 		
  else if (HadronicPhysicsType == "QGSP_HP") 		ConstructPartonStringBasedPhysics(false,false,false,false,true);
  else if (HadronicPhysicsType == "QGSP_BIC")		ConstructPartonStringBasedPhysics(false,false,true,false,false);
  else if (HadronicPhysicsType == "QGSP_BIC_HP")	ConstructPartonStringBasedPhysics(false,false,true,false,true);
  else if (HadronicPhysicsType == "QGSP_BERT")		ConstructPartonStringBasedPhysics(false,false,false,true,false);
  else if (HadronicPhysicsType == "QGSP_BERT_HP")	ConstructPartonStringBasedPhysics(false,false,false,true,true);
  
  if (HadronicPhysicsType == "QGSC")			ConstructPartonStringBasedPhysics(false,true,false,false,false); 		
  else if (HadronicPhysicsType == "QGSC_HP") 		ConstructPartonStringBasedPhysics(false,true,false,false,true);
  else if (HadronicPhysicsType == "QGSC_BIC")		ConstructPartonStringBasedPhysics(false,true,true,false,false);
  else if (HadronicPhysicsType == "QGSC_BIC_HP")	ConstructPartonStringBasedPhysics(false,true,true,false,true);
  else if (HadronicPhysicsType == "QGSC_BERT")		ConstructPartonStringBasedPhysics(false,true,false,true,false);
  else if (HadronicPhysicsType == "QGSC_BERT_HP")	ConstructPartonStringBasedPhysics(false,true,false,true,true);
  
  if (HadronicPhysicsType == "FTFP")			ConstructPartonStringBasedPhysics(true,false,false,false,false); 		
  else if (HadronicPhysicsType == "FTFP_HP") 		ConstructPartonStringBasedPhysics(true,false,false,false,true);
  else if (HadronicPhysicsType == "FTFP_BIC")		ConstructPartonStringBasedPhysics(true,false,true,false,false);
  else if (HadronicPhysicsType == "FTFP_BIC_HP")	ConstructPartonStringBasedPhysics(true,false,true,false,true);
  else if (HadronicPhysicsType == "FTFP_BERT")		ConstructPartonStringBasedPhysics(true,false,false,true,false);
  else if (HadronicPhysicsType == "FTFP_BERT_HP")	ConstructPartonStringBasedPhysics(true,false,false,true,true);
  
  if (HadronicPhysicsType == "FTFC")			ConstructPartonStringBasedPhysics(true,true,false,false,false); 		
  else if (HadronicPhysicsType == "FTFC_HP") 		ConstructPartonStringBasedPhysics(true,true,false,false,true);
  else if (HadronicPhysicsType == "FTFC_BIC")		ConstructPartonStringBasedPhysics(true,true,true,false,false);
  else if (HadronicPhysicsType == "FTFC_BIC_HP")	ConstructPartonStringBasedPhysics(true,true,true,false,true);
  else if (HadronicPhysicsType == "FTFC_BERT")		ConstructPartonStringBasedPhysics(true,true,false,true,false);
  else if (HadronicPhysicsType == "FTFC_BERT_HP")	ConstructPartonStringBasedPhysics(true,true,false,true,true);
  
  //Build
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  
  //ElasticPhysics
  //---------------
  ConstructElasticPhysics();
}


////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSHadronicPhysics::SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ)
{//G4cout<<anA<<" "<<aZ<<" Max A and Max Z"<<std::endl;
 if (!theBinaryProton) {
 	theBinaryProton = new MYBinaryProtonBuilder;
  	theBinaryNeutron = new MYBinaryNeutronBuilder;
 }	
 theBinaryProton->SetMaxAandZForFermiBreakUp(anA,aZ);
 theBinaryNeutron->SetMaxAandZForFermiBreakUp(anA,aZ);
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSHadronicPhysics::SetMinEForMultiFrag(G4double anE)
{if (!theBinaryProton) {
 	theBinaryProton = new MYBinaryProtonBuilder;
  	theBinaryNeutron = new MYBinaryNeutronBuilder;
 }
 theBinaryProton->SetMinEForMultiFrag(anE);
 theBinaryNeutron->SetMinEForMultiFrag(anE);
}

///////////////////////////////////////////////////////////////////////////////
void PLANETOCOSHadronicPhysics::ConstructLHEPBasedPhysics(G4bool bic_flag,G4bool bert_flag,G4bool hp_flag)
{  G4LHEPNeutronBuilder* theLHEPNeutron = new G4LHEPNeutronBuilder();
   theNeutrons.RegisterMe(theLHEPNeutron);
   
   if (hp_flag) theLHEPNeutron->SetMinEnergy(19.9*MeV);
   if (bic_flag){
   	if (!theBinaryProton) SetMaxAandZForFermiBreakUp(13,7);
   	theLHEPNeutron->SetMinInelasticEnergy(9.5*GeV);
 	theBinaryNeutron->SetMaxEnergy(9.9*GeV);
	theNeutrons.RegisterMe(theBinaryNeutron);
	if (hp_flag) theBinaryNeutron->SetMinEnergy(19.9*MeV);
   }
   else if(bert_flag){
   	
	G4BertiniNeutronBuilder* theBertiniNeutron = new G4BertiniNeutronBuilder();
 	theLHEPNeutron->SetMinInelasticEnergy(9.5*GeV);
 	theBertiniNeutron->SetMaxEnergy(9.9*GeV);
 	theNeutrons.RegisterMe(theBertiniNeutron);
	if (hp_flag)	theBertiniNeutron->SetMinEnergy(19.9*MeV);
   }
   if (hp_flag) theNeutrons.RegisterMe(new G4NeutronHPBuilder());
   
   
   G4LHEPProtonBuilder* theLHEPProton = new G4LHEPProtonBuilder();
   thePro.RegisterMe(theLHEPProton);
   if (bic_flag){
   	theLHEPProton->SetMinEnergy(9.5*GeV);
 	theBinaryProton->SetMaxEnergy(9.9*GeV);
	thePro.RegisterMe(theBinaryProton);
   }
   else if(bert_flag){
   	G4BertiniProtonBuilder* theBertiniProton = new G4BertiniProtonBuilder();
 	theLHEPProton->SetMinEnergy(9.5*GeV);
 	theBertiniProton->SetMaxEnergy(9.9*GeV);
 	thePro.RegisterMe(theBertiniProton);
   }
   
   G4LHEPPiKBuilder* theLHEPPiK =new G4LHEPPiKBuilder();
   thePiK.RegisterMe(theLHEPPiK);
   
   if(bert_flag){
   	G4BertiniPiKBuilder* theBertiniPiK=new G4BertiniPiKBuilder();
   	thePiK.RegisterMe(theBertiniPiK);
   	theLHEPPiK->SetMinEnergy(9.5*GeV);
   	theBertiniPiK->SetMaxEnergy(9.9*GeV);
   }	
   
   NeutronHPElasticFlag =hp_flag;
   
 
 
 
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSHadronicPhysics::ConstructPartonStringBasedPhysics(G4bool ftf_flag, G4bool chips_flag,G4bool bic_flag,G4bool bert_flag,G4bool hp_flag)
{ 
  
  //flag for the different models
  //-----------------------------
  G4bool QGSC_flag = ((!ftf_flag) && chips_flag);
  G4bool QGSP_flag = ((!ftf_flag) && (!chips_flag));
  G4bool FTFC_flag = ( ftf_flag && chips_flag);
  G4bool FTFP_flag = ( ftf_flag && (!chips_flag));
  
  G4bool cascade_flag = (bic_flag ||  bert_flag); //one cacade model is used eithe Binary or Bertini
  G4bool LEP_flag = (!QGSC_flag ||  !cascade_flag);  //the LEP model will also be used
  
  G4bool Only_LEP_PiK_flag = (!LEP_flag && bic_flag); 
  
 
 
 
 
 //high energy models
 //---------------
  
  
  if (FTFC_flag) {
  	theNeutrons.RegisterMe(new G4FTFCNeutronBuilder());
	thePro.RegisterMe(new G4FTFCProtonBuilder());
	thePiK.RegisterMe(new G4FTFCPiKBuilder());
 }	
  else if (FTFP_flag) {
  	theNeutrons.RegisterMe(new G4FTFPNeutronBuilder());
	thePro.RegisterMe(new G4FTFPProtonBuilder());
	thePiK.RegisterMe(new G4FTFPPiKBuilder());
 }
  else if (QGSC_flag) {
  	theNeutrons.RegisterMe(new G4QGSCNeutronBuilder());
	thePro.RegisterMe(new G4QGSCProtonBuilder());
	thePiK.RegisterMe(new G4QGSCPiKBuilder());
  }	
  else if (QGSP_flag) {
  	theNeutrons.RegisterMe(new G4QGSPNeutronBuilder());
	thePro.RegisterMe(new G4QGSPProtonBuilder());
	thePiK.RegisterMe(new G4QGSPPiKBuilder());
  }
  
  
  //LEP model if needed
  //--------------------
  G4LEPNeutronBuilder*  theLEPNeutron =0;
  G4LEPProtonBuilder*  theLEPProton =0;
  G4LEPPiKBuilder*  theLEPPiK =0;
  if (LEP_flag) {
  	theLEPNeutron=new G4LEPNeutronBuilder;
	theLEPNeutron->SetMaxInelasticEnergy(25*GeV);
  	theNeutrons.RegisterMe(theLEPNeutron);
	if (hp_flag) theLEPNeutron->SetMinEnergy(19.9*MeV);
	
	theLEPProton=new G4LEPProtonBuilder;
	theLEPProton->SetMaxEnergy(25*GeV);
  	thePro.RegisterMe(theLEPProton);
	
	theLEPPiK=new G4LEPPiKBuilder;
	theLEPPiK->SetMaxEnergy(25*GeV);
  	thePiK.RegisterMe(theLEPPiK);
	
  }
  else if (Only_LEP_PiK_flag) {
  	theLEPPiK=new G4LEPPiKBuilder;
	theLEPPiK->SetMaxEnergy(25*GeV);
  	thePiK.RegisterMe(theLEPPiK);
  
  
  }
  
  //cascade models
  //-------------
  if (cascade_flag){
   	if (LEP_flag) {
	 	theLEPNeutron->SetMinInelasticEnergy(9.5*GeV);
		theLEPProton->SetMinEnergy(9.5*GeV);
		if (bert_flag) theLEPPiK->SetMinEnergy(9.5*GeV);
	}	
	if (bic_flag) { //binary
		
		if (!theBinaryProton) SetMaxAandZForFermiBreakUp(13,7);
		theBinaryNeutron->SetMaxEnergy(9.9*GeV);
		theNeutrons.RegisterMe(theBinaryNeutron);
		if (hp_flag) theBinaryNeutron->SetMinEnergy(19.9*MeV);
		
		theBinaryProton->SetMaxEnergy(9.9*GeV);
		thePro.RegisterMe(theBinaryProton);
	}
	else { //bertini
		G4BertiniNeutronBuilder* theBertiniNeutron = new G4BertiniNeutronBuilder();
		theBertiniNeutron->SetMaxEnergy(9.9*GeV);
		theNeutrons.RegisterMe(theBertiniNeutron);
		if (hp_flag) theBertiniNeutron->SetMinEnergy(19.9*MeV);
		
		G4BertiniProtonBuilder* theBertiniProton = new G4BertiniProtonBuilder();
		theBertiniProton->SetMaxEnergy(9.9*GeV);
		thePro.RegisterMe(theBertiniProton);
		
		G4BertiniPiKBuilder* theBertiniPiK = new G4BertiniPiKBuilder();
		theBertiniPiK->SetMaxEnergy(9.9*GeV);
		thePiK.RegisterMe(theBertiniPiK);
		
		
	}	
  }
  
  if (hp_flag) theNeutrons.RegisterMe(new G4NeutronHPBuilder());
  
}

///////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSHadronicPhysics::ConstructElasticPhysics()
{ G4HadronElasticPhysics* theElasticPhysics =
			new G4HadronElasticPhysics("elastic",1,NeutronHPElasticFlag);
  //theElasticPhysics->ConstructParticle();
  theElasticPhysics->ConstructProcess();
}
