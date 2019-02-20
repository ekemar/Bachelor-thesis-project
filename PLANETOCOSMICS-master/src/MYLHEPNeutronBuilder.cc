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
 #include "MYLHEPNeutronBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"
 #include "G4Element.hh"

 MYLHEPNeutronBuilder::
 MYLHEPNeutronBuilder() 
 {
   theMin = 0;
   theIMin = 0;
 }

 MYLHEPNeutronBuilder::
 ~MYLHEPNeutronBuilder() {}
#ifdef   BEFORE_V8 
 void MYLHEPNeutronBuilder::
 Build(G4NeutronInelasticProcess & aP)
 {
   theLENeutronModel = new G4LENeutronInelastic();
   theHENeutronModel = new G4HENeutronInelastic();
   theLENeutronModel->SetMinEnergy(theIMin);
   theLENeutronModel->SetMaxEnergy(55*GeV);
   theHENeutronModel->SetMinEnergy(25*GeV);
   aP.RegisterMe(theLENeutronModel);
   aP.RegisterMe(theHENeutronModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronFissionProcess & aP)
 {
   theNeutronFissionModel = new G4LFission();
   theNeutronFissionModel->SetMinEnergy(theMin);
   aP.RegisterMe(theNeutronFissionModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronElasticProcess & aP)
 {
//replaced by L. Desorgher to consider a realistic np elastic physics 
   theLElasticModel = new G4LElastic();
   theLElasticModel->SetMinEnergy(theMin);
   theLElasticModel->SetMinEnergy(0,G4Element::GetElement ("Hydrogen"));
   theLElasticModel->SetMaxEnergy(0,G4Element::GetElement ("Hydrogen"));
   aP.RegisterMe(theLElasticModel);
  theLEnpModel = new G4LEnp();
  theLEnpModel->SetMinEnergy(0.);
  theLEnpModel->SetMaxEnergy(0.);
  theLEnpModel->SetMinEnergy(theMin,G4Element::GetElement ("Hydrogen"));
  theLEnpModel->SetMaxEnergy(20.*TeV,G4Element::GetElement ("Hydrogen"));
  aP.RegisterMe(theLEnpModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronCaptureProcess & aP)
 {
   theNeutronCaptureModel = new G4LCapture();
   theNeutronCaptureModel->SetMinEnergy(theMin);
   aP.RegisterMe(theNeutronCaptureModel);
 }
#else
void MYLHEPNeutronBuilder::
 Build(G4NeutronInelasticProcess * aP)
 {
   theLENeutronModel = new G4LENeutronInelastic();
   theHENeutronModel = new G4HENeutronInelastic();
   theLENeutronModel->SetMinEnergy(theIMin);
   theLENeutronModel->SetMaxEnergy(55*GeV);
   theHENeutronModel->SetMinEnergy(25*GeV);
   aP->RegisterMe(theLENeutronModel);
   aP->RegisterMe(theHENeutronModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronFissionProcess * aP)
 {
   theNeutronFissionModel = new G4LFission();
   theNeutronFissionModel->SetMinEnergy(theMin);
   aP->RegisterMe(theNeutronFissionModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronElasticProcess * aP)
 {
//replaced by L. Desorgher to consider a realistic np elastic physics 
   theLElasticModel = new G4LElastic();
   theLElasticModel->SetMinEnergy(theMin);
   theLElasticModel->SetMinEnergy(0,G4Element::GetElement ("Hydrogen"));
   theLElasticModel->SetMaxEnergy(0,G4Element::GetElement ("Hydrogen"));
   aP->RegisterMe(theLElasticModel);
  theLEnpModel = new G4LEnp();
  theLEnpModel->SetMinEnergy(0.);
  theLEnpModel->SetMaxEnergy(0.);
  theLEnpModel->SetMinEnergy(theMin,G4Element::GetElement ("Hydrogen"));
  theLEnpModel->SetMaxEnergy(20.*TeV,G4Element::GetElement ("Hydrogen"));
  aP->RegisterMe(theLEnpModel);
 }

 void MYLHEPNeutronBuilder::
 Build(G4HadronCaptureProcess * aP)
 {
   theNeutronCaptureModel = new G4LCapture();
   theNeutronCaptureModel->SetMinEnergy(theMin);
   aP->RegisterMe(theNeutronCaptureModel);
 }
#endif
 // 2002 by J.P. Wellisch
