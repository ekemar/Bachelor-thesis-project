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
#ifndef MYLHEPNeutronBuilder_h
#define MYLHEPNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4LElastic.hh" 
#include "G4LEnp.hh"  
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LENeutronInelastic.hh"
#include "G4HENeutronInelastic.hh"

class MYLHEPNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    MYLHEPNeutronBuilder();
    virtual ~MYLHEPNeutronBuilder();

  public: 
#ifdef   BEFORE_V8
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
    void SetMinEnergy(G4double aM) 
    {
      theMin = aM;
      theIMin = aM;
    }
    void SetMinInelasticEnergy(G4double aM) 
    {
      theIMin = aM;
    }

  private:
   //added by L. Desorgher to consider a realistic np elastic physics 
  
    G4LEnp* theLEnpModel;
    G4LElastic * theLElasticModel;
    G4LENeutronInelastic * theLENeutronModel;
    G4HENeutronInelastic * theHENeutronModel;
    G4LFission * theNeutronFissionModel;
    G4LCapture * theNeutronCaptureModel;
    
    G4double theMin;
    G4double theIMin;

};

// 2002 by J.P. Wellisch

#endif

