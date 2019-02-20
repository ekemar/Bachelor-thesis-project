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
#ifndef MYBinaryNeutronBuilder_h
#define MYBinaryNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4BinaryCascade.hh"   
#include "G4NeutronInelasticCrossSection.hh"
#include "G4ExcitationHandler.hh"
class MYBinaryNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    MYBinaryNeutronBuilder();
    virtual ~MYBinaryNeutronBuilder();

  public: 
#ifdef  BEFORE_V8 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
#else
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4HadronFissionProcess * aP);
    virtual void Build(G4HadronCaptureProcess * aP);
    virtual void Build(G4NeutronInelasticProcess * aP);
#endif    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}
    inline void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ)
    	 {theExcitationHandler->SetMaxAandZForFermiBreakUp(anA,aZ);}
    inline void SetMinEForMultiFrag(G4double anE)
         {theExcitationHandler->SetMinEForMultiFrag(anE);} 
  private:
    G4BinaryCascade * theModel;    
    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;
    G4double theMax;
    G4ExcitationHandler* theExcitationHandler ;
};

// 2002 by J.P. Wellisch

#endif

