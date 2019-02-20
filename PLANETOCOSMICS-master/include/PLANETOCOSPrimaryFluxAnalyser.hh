#ifndef PLANETOCOSPrimaryFluxAnalyser_HH
#define PLANETOCOSPrimaryFluxAnalyser_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"PLANETOCOSPrimaryHit.hh"
#include"IHistoDefinition.hh"


class G4Step;
class PLANETOCOSPrimaryGeneratorAction;
class PLANETOCOSAnalysisManager;
class PLANETOCOSPrimaryFluxAnalyser
{
public:
   
   PLANETOCOSPrimaryFluxAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager);
   ~PLANETOCOSPrimaryFluxAnalyser();
   
   //Analyse
   //--------
   void Analyse(PLANETOCOSPrimaryHitsCollection* FluxHC);
   
   //Get methods
   //-----------
   inline bool GetDetectPrimary(){return DetectPrimary;}
   
   //Set methods
   //-----------
  
   
   
  //Primary histograms
  //----------------			 
  
   void CreateEnergyVsZenithPrimaryHisto(G4String aParticleName,G4String label,
			         G4int nE,G4double E1,G4double E2, G4String ScaleType, 
			         G4int nZenith, G4double cosZen1, G4double cosZen2);
  
  
   void CreateZenithVsAzimuthPrimaryHisto(G4String aParticleName,G4String label,
			          G4int nZenith, G4double cosZen1, G4double cosZen2,
				  G4int nAzim, G4double Azim1, G4double Azim2);
			 	 
   void OmnifluxFluxPrimaryHisto(G4String aParticleName,G4String label,
			         G4int nE,G4double E1,G4double E2, 
				 G4String ScaleType); 			 
   
   void CreateZenithPrimaryHisto(G4String aParticleName,G4String label,
			         G4int nZenith, G4double cosZen1, 
				 G4double cosZen2);
			 
   void CreateAzimuthPrimaryHisto(G4String aParticleName,G4String label,
			          G4int nAzim, G4double Azim1, 
				  G4double Azim2);  
   
   void Initialise();
   
   
   
private:

  PLANETOCOSAnalysisManager* theAnalysisManager; 
  
  
  
  
  //Primary histograms
  //------------------
  G4bool DetectPrimary;
  std::vector< HISTO2D*>  EnergyVsZenithPrimaryHisto;
  std::vector< G4int>  EnergyVsZenithPrimaryPDGCode;
  G4bool DetectEnergyVsZenithPrimary;
  
  
  std::vector< HISTO2D*>  ZenithVsAzimuthPrimaryHisto;
  std::vector< G4int>  ZenithVsAzimuthPrimaryPDGCode;
  G4bool DetectZenithVsAzimuthPrimary;
  
  std::vector< HISTO1D*>  OmniFluxPrimaryHisto;
  std::vector< G4int>   OmniFluxPrimaryPDGCode;
  G4bool DetectOmniFluxPrimary;
  
  std::vector< HISTO1D*>  ZenithPrimaryHisto;
  std::vector< G4int>   ZenithPrimaryPDGCode;
  G4bool DetectZenithPrimary;
  
  std::vector< HISTO1D*>  AzimuthPrimaryHisto;
  std::vector< G4int>   AzimuthPrimaryPDGCode;
  G4bool DetectAzimuthPrimary;
  
  
  //Bin content for histogram 
  // 1 is for unscaled case
  // 2 is for scaling to primary flux
  // 3 is for scaling per pimary part
  //-------------------
  
  
  G4String EdepBin1;
  G4String EdepBin2;
  G4String EdepBin3;
  
  G4String DirDifFluxBin1;
  G4String DirDifFluxBin2;
  G4String DirDifFluxBin3;
  
  G4String OmniDifFluxBin1;
  G4String OmniDifFluxBin2;
  G4String OmniDifFluxBin3;
  
  G4String IntegFluxBin1;
  G4String IntegFluxBin2;
  G4String IntegFluxBin3;
  
  //Common actions
  //------------------   
   
  G4int GetPDGCodeAndInitialise(G4String aParticleName); 
   
   
};
#endif
