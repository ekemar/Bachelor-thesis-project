//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
 #include "MYBinaryNeutronBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"
 #include "G4PreCompoundModel.hh"

 MYBinaryNeutronBuilder::
 MYBinaryNeutronBuilder() 
 {
   theMin = 0;
   theMax = 9.9*GeV;
   theModel =0;
   theExcitationHandler = new G4ExcitationHandler;
//  theExcitationHandler->SetMaxAandZForFermiBreakUp(17,9);
//  theExcitationHandler->SetMinEForMultiFrag(5.*MeV);
 }
#ifdef  BEFORE_V8 
 void MYBinaryNeutronBuilder::
 Build(G4NeutronInelasticProcess & aP)
 {
   if (!theModel) {
  	theModel = new G4BinaryCascade;
  	theModel->SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
  }
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP.RegisterMe(theModel);
   aP.AddDataSet(&theXSec);  
 }

 MYBinaryNeutronBuilder::
 ~MYBinaryNeutronBuilder() {}

 void MYBinaryNeutronBuilder::
 Build(G4HadronElasticProcess & )
 {
 }

 void MYBinaryNeutronBuilder::
 Build(G4HadronFissionProcess & )
 {
 }

 void MYBinaryNeutronBuilder::
 Build(G4HadronCaptureProcess & )
 {
 }
#else
 void MYBinaryNeutronBuilder::
 Build(G4NeutronInelasticProcess * aP)
 { if (!theModel) {
  	theModel = new G4BinaryCascade;
  	theModel->SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
  }
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP->RegisterMe(theModel);
   aP->AddDataSet(&theXSec);  
 }

 MYBinaryNeutronBuilder::
 ~MYBinaryNeutronBuilder() {}

 void MYBinaryNeutronBuilder::
 Build(G4HadronElasticProcess * )
 {
 }

 void MYBinaryNeutronBuilder::
 Build(G4HadronFissionProcess * )
 {
 }

 void MYBinaryNeutronBuilder::
 Build(G4HadronCaptureProcess * )
 {
 } 
#endif 

 // 2002 by J.P. Wellisch
