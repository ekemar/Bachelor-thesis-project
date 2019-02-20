#ifndef PLANETOCOSSoilEdepAnalyser_HH
#define PLANETOCOSSoilEdepAnalyser_HH 1
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

class PLANETOCOSSoilEdepAnalyser
{
public: 

   PLANETOCOSSoilEdepAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager);
   ~PLANETOCOSSoilEdepAnalyser();
   
   
   //Analyse
   //--------
   void Analyse(PLANETOCOSEdepHitsCollection* EdepHC); 
    
   
   //Get methods
   //-----------
   inline bool GetDetectEdep(){return DetectEdep;}
     
   //Set methods
   //-----------
   inline void SetSoilLayerDepthsInLengthUnit(std::vector<double> *pvec)
                                                { SoilLayerDepthsInLengthUnit= pvec;} 
   inline void SetSoilLayerDepthsInDepthUnit(std::vector<double> *pvec)
  		                                { SoilLayerDepthsInDepthUnit= pvec;}
   
   //Edep histograms
   //----------------
   void CreateEdepVsThickness(G4String label,G4int n_Bins);
   void CreateEdepVsDepth(G4String label,G4int n_Bins);

private:

  PLANETOCOSAnalysisManager* theAnalysisManager; 
  
 
  
   //deposited energy histograms
   //---------------------------
   G4bool DetectEdep;
   HISTO1D* EdepVsThicknessHisto;
   HISTO1D* EdepVsDepthHisto;
  
   
   G4bool DetectEdepVsThickness;
   G4bool DetectEdepVsDepth;
 
  
   std::vector<G4double> DepthForThicknessHisto;
   G4double delta_depth;
   std::vector<G4double>* SoilLayerDepthsInLengthUnit;
   std::vector<G4double>* SoilLayerDepthsInDepthUnit; 
   
   
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
