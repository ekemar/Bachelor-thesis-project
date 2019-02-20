#ifndef PLANETOCOSEdepAnalyser_HH
#define PLANETOCOSEdepAnalyser_HH 1
#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"PLANETOCOSEdepHit.hh"
#include"IHistoDefinition.hh"


class G4Step;
class PLANETOCOSPrimaryGeneratorAction;
class PLANETOCOSAnalysisManager;

class PLANETOCOSEdepAnalyser
{
public: 

   PLANETOCOSEdepAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager);
   ~PLANETOCOSEdepAnalyser();
   
   
   //Analyse
   //--------
   void Analyse(PLANETOCOSEdepHitsCollection* EdepHC); 
    
   
   //Get methods
   //-----------
   inline bool GetDetectEdep(){return DetectEdep;}
     
   //Set methods
   //-----------
   inline void SetAtmosphericLayerAltitudes(std::vector<double> *pvec)
                                                {AtmosphericLayerAltitudes = pvec;} 
   inline void SetAtmosphericLayerDepths(std::vector<double> *pvec)
  		                                {AtmosphericLayerDepths = pvec;}
   
   //Edep histograms
   //----------------
   void CreateEdepVsAltitude(G4String label,G4int n_Alt);
   void CreateEdepVsDepth(G4String label,G4int n_Depth);

private:

  PLANETOCOSAnalysisManager* theAnalysisManager; 
  
 
  
   //deposited energy histograms
   //---------------------------
   G4bool DetectEdep;
   HISTO1D* EdepVsAltitudeHisto;
   HISTO1D* EdepVsDepthHisto;
#ifdef  TEST_EDEP    
   
   std::vector< HISTO1D* > TestEdepVsDepthHisto;
#endif   
   
   G4bool DetectEdepvsAltitude;
   G4bool DetectEdepvsDepth;
 
  
   std::vector<G4double> DepthForAltitudeHisto;
   G4double delta_depth;
   std::vector<G4double>* AtmosphericLayerAltitudes;
   std::vector<G4double>* AtmosphericLayerDepths; 
   
   
  //Bin content for histogram 
  // 1 is for unscaled case
  // 2 is for scaling to primary flux
  // 3 is for scaling per pimary part
  //-------------------
  
  
  G4String EdepBin1;
  G4String EdepBin2;
  G4String EdepBin3;
  
};
#endif
