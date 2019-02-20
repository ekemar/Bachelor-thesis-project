#ifndef PLANETOCOSPhysicsList_h
#define PLANETOCOSPhysicsList_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSVModularPhysicsList.hh"
#include "globals.hh"
#include <vector>
#include "PLANETOCOSHadronicPhysics.hh"


class PLANETOCOSGeneralPhysics;
class G4ProductionCuts;
class PLANETOCOSPhysicsListMessenger;
class PLANETOCOSGeometryConstruction;

/*
class G4Regin;

struct Cuts {
  G4String regName;
  G4double gamma;
  G4double electron;
  G4double positron;
};
typedef  std::vector<Cuts> RegCuts;

*/
////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSPhysicsList: public PLANETOCOSVModularPhysicsList
{
public:
  PLANETOCOSPhysicsList();
  virtual ~PLANETOCOSPhysicsList();
  
public:

  // Type of physics
  inline void SelectHadronicPhysicsType (G4String aString)
                                {HadronicPhysicsType=aString;}
  inline void SelectEMPhysicsType (G4String aString)
                                {EMPhysicsType=aString;}
  inline void SelectLightIonHadronicPhysicsType(G4String aString)
                                {LightIonHadronicPhysicsType=aString;}
  inline void SetConsiderElectromagneticNuclearPhysics(G4bool aBool) 
  				{ConsiderElectroNucPhysics=aBool;}
  
  inline void SetConsiderMuonNuclearPhysics(G4bool aBool) 
  				{ConsiderMuonNucPhysics=aBool;}
 
  inline void SetConsiderSynchrotronRadiation(G4bool aBool) 
  				{ConsiderSyncPhysics=aBool;}				
  				
  inline void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ)
  			{//G4cout<<" test"<<std::endl;
			 theHadronicPhysics->SetMaxAandZForFermiBreakUp(anA,aZ);}
  inline void SetMinEForMultiFrag(G4double anE)
  			{theHadronicPhysics->SetMinEForMultiFrag(anE);}				
  
  void ShowTypeOfPhysics(); 
  
  
  // Cut by regions
  void SetCutInDepthForAllLayers(G4double );
  void SetCutInDepthForASpecificLayer(G4double , G4String);
  void SetCutInLengthForASpecificLayer(G4double , G4String);
  void SetCutInDepthForAllLayersAndForParticle(G4double, G4String );
  void DeleteAllRegions();
 
 // cuts
  virtual void SetCuts ();

  inline void SetGlobalDefault(G4double d) { cutForPositron 
					= cutForElectron = cutForGamma 
					=  defaultCutValue = d;};
  inline void SetGammaDefault (G4double d) {cutForGamma = d;};
  inline void SetElectronDefault (G4double d) {cutForElectron = d;};
  inline void SetPositronDefault (G4double d) {cutForPositron = d;};


public:

  void BuildList ();

  void SetGlobalCuts ();
  void SetCutsByRegion ();
  
  inline void SetSteppingAlgorithmParameters(bool aBool, double fac){
    	SteppingAlgorithmMsc=aBool;
	facrange=fac;
  }

private:

  PLANETOCOSPhysicsListMessenger* physListMessenger;
  
  G4String EMPhysicsType;
  G4String HadronicPhysicsType;
  G4String LightIonHadronicPhysicsType;
  G4bool ConsiderElectroNucPhysics;
  G4bool ConsiderMuonNucPhysics;
  G4bool ConsiderSyncPhysics;
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  
  std::vector<G4ProductionCuts*> production_cuts;
  //  std::vector<G4Regin*> mlRegion;
  //  RegCuts  mlRegCuts;

  PLANETOCOSGeometryConstruction* geometry;
  PLANETOCOSHadronicPhysics* theHadronicPhysics;
  PLANETOCOSGeneralPhysics* theGeneralPhysics;
  G4bool SteppingAlgorithmMsc;
  G4double facrange;

};
////////////////////////////////////////////////////////////////////////////////
#endif
