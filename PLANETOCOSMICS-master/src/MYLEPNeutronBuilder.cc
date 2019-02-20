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
#include "MYLEPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh" 
#include "G4Element.hh"

MYLEPNeutronBuilder::
MYLEPNeutronBuilder() 
{
  theMin = 0;
  theIMin = theMin;
  theMax = 20*TeV;
  theIMax = 55*GeV;
}

MYLEPNeutronBuilder::
~MYLEPNeutronBuilder() {}

#ifdef   BEFORE_V8 
void MYLEPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
  theLElasticModel = new G4LElastic();
  theLElasticModel->SetMinEnergy(theMin);
  theLElasticModel->SetMaxEnergy(theMax);
  theLElasticModel->SetMinEnergy(0.,G4Element::GetElement ("Hydrogen"));
  theLElasticModel->SetMaxEnergy(0.,G4Element::GetElement ("Hydrogen"));
  
  theLEnpModel = new G4LEnp();
  theLEnpModel->SetMinEnergy(0.);
  theLEnpModel->SetMaxEnergy(0.);
  theLEnpModel->SetMinEnergy(theMin,G4Element::GetElement ("Hydrogen"));
  theLEnpModel->SetMaxEnergy(theMax,G4Element::GetElement ("Hydrogen"));
  

  
  aP.RegisterMe(theLElasticModel);
  aP.RegisterMe(theLEnpModel);
}

void MYLEPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  theNeutronFissionModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theNeutronFissionModel);
}

void MYLEPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  theNeutronCaptureModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theNeutronCaptureModel);
}

void MYLEPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theLENeutronModel = new G4LENeutronInelastic();
  theLENeutronModel->SetMinEnergy(theIMin);
  theLENeutronModel->SetMaxEnergy(theIMax);
  aP.RegisterMe(theLENeutronModel);
}
#else
void MYLEPNeutronBuilder::
Build(G4HadronElasticProcess * aP)
{
  theLElasticModel = new G4LElastic();
  theLElasticModel->SetMinEnergy(theMin);
  theLElasticModel->SetMaxEnergy(theMax);
  theLElasticModel->SetMinEnergy(0.,G4Element::GetElement ("Hydrogen"));
  theLElasticModel->SetMaxEnergy(0.,G4Element::GetElement ("Hydrogen"));
  
  theLEnpModel = new G4LEnp();
  theLEnpModel->SetMinEnergy(0.);
  theLEnpModel->SetMaxEnergy(0.);
  theLEnpModel->SetMinEnergy(theMin,G4Element::GetElement ("Hydrogen"));
  theLEnpModel->SetMaxEnergy(theMax,G4Element::GetElement ("Hydrogen"));
  

  
  aP->RegisterMe(theLElasticModel);
  aP->RegisterMe(theLEnpModel);
}

void MYLEPNeutronBuilder::
Build(G4HadronFissionProcess * aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  theNeutronFissionModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theNeutronFissionModel);
}

void MYLEPNeutronBuilder::
Build(G4HadronCaptureProcess * aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  theNeutronCaptureModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theNeutronCaptureModel);
}

void MYLEPNeutronBuilder::
Build(G4NeutronInelasticProcess * aP)
{
  theLENeutronModel = new G4LENeutronInelastic();
  theLENeutronModel->SetMinEnergy(theIMin);
  theLENeutronModel->SetMaxEnergy(theIMax);
  aP->RegisterMe(theLENeutronModel);
}
#endif

// 2002 by J.P. Wellisch
