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
#ifndef MYBinaryProtonBuilder_h
#define MYBinaryProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4BinaryCascade.hh"   
#include "G4ExcitationHandler.hh"   
#include "G4ProtonInelasticCrossSection.hh"

class MYBinaryProtonBuilder : public G4VProtonBuilder
{
  public: 
    MYBinaryProtonBuilder();
    virtual ~MYBinaryProtonBuilder();

  public: 
#ifdef  BEFORE_V8
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
#else 
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4ProtonInelasticProcess * aP); 
#endif    
    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}
    inline void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ)
    	 {theExcitationHandler->SetMaxAandZForFermiBreakUp(anA,aZ);}
    inline void SetMinEForMultiFrag(G4double anE)
         {theExcitationHandler->SetMinEForMultiFrag(anE);}  
    

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4BinaryCascade * theModel;   
    G4ExcitationHandler* theExcitationHandler ;
    G4double theMin;
    G4double theMax;
};

// 2002 by J.P. Wellisch

#endif

