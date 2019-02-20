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
#include "MYBinaryProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4PreCompoundModel.hh"

MYBinaryProtonBuilder::
MYBinaryProtonBuilder() 
{
  
  theMin = 0;
  theMax = 9.9*GeV;
  theModel = 0;
  theExcitationHandler = new G4ExcitationHandler;
  //theExcitationHandler->SetMaxAandZForFermiBreakUp(17,9);
  //theExcitationHandler->SetMinEForMultiFrag(5.*MeV);
}
#ifdef  BEFORE_V8
void MYBinaryProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{ 
  if (!theModel) {
  	theModel = new G4BinaryCascade;
  	theModel->SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
  }
  aP.AddDataSet(&theXSec);  
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
}

MYBinaryProtonBuilder::
~MYBinaryProtonBuilder() {}

void MYBinaryProtonBuilder::
Build(G4HadronElasticProcess & )
{
}
#else
void MYBinaryProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
{
  if (!theModel) {
  	theModel = new G4BinaryCascade;
  	theModel->SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
 }	
	
  aP->AddDataSet(&theXSec);  
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theModel);
}

MYBinaryProtonBuilder::
~MYBinaryProtonBuilder() {}

void MYBinaryProtonBuilder::
Build(G4HadronElasticProcess * )
{
}
#endif

// 2002 by J.P. Wellisch
