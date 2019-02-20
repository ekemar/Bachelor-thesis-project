#ifndef PLANETOCOSPseudoTrappingAnalyser_HH
#define PLANETOCOSPseudoTrappingAnalyser_HH 1
#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"PLANETOCOSPostTrackHit.hh"
#include"IHistoDefinition.hh"
/*namespace AIDA
{class ITree;
 class HISTOFactory;
 class HISTO1D;
 class HISTO2D;
 class HISTO;
 class IAnalysisFactory;
 class ITupleFactory;
 class ITuple;}

using namespace AIDA;
*/

class G4Step;
class PLANETOCOSPrimaryGeneratorAction;
class PLANETOCOSAnalysisManager;

class PLANETOCOSPseudoTrappingAnalyser
{
public: 

   PLANETOCOSPseudoTrappingAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager);
   ~PLANETOCOSPseudoTrappingAnalyser();
   
   
   //Analyse
   //--------
   void Analyse(PLANETOCOSPostTrackHitsCollection* PostTrackHC); 
    
   
   //Get methods
   //-----------
   inline bool GetDetectPseudoTrapping(){return DetectPseudoTrapping;}
     
   //Set methods
   //-----------
 
   
   //PostTrackHit Histograms
   //-----------------------
   void CreateNbOfPlanetTurnVsStartEkinHisto(G4String aParticleName,G4String label,
			                     G4int nE,G4double E1,G4double E2, 
					     G4String ScaleType,
			                     G4int nBinTurn,G4double nTurnMax );
   void CreateNbOfPlanetTurnVsLifeTimeHisto(G4String aParticleName,G4String label,
			                    G4int nTime,G4double time1,G4double time2, 
					    G4String ScaleType,
			                    G4int nBinTurn,G4double nTurnMax );
   void CreateNbOfEquatorCrossingVsStartEkinHisto(G4String aParticleName,G4String label,
			                     G4int nE,G4double E1,G4double E2, 
					     G4String ScaleType,
			                     G4int nCrossingMax );
   void CreateNbOfEquatorCrossingVsLifeTimeHisto(G4String aParticleName,G4String label,
			                    G4int nTime,G4double time1,G4double time2, 
					    G4String ScaleType,
			                    G4int nCrossingMax );
   
   
   void CreateNbOfEquatorCrossingVsNbOfPlanetTurnHisto(G4String aParticleName,G4String label,
			                    G4int nBinTurn,G4double nTurnMax, 
			                    G4int nCrossingMax );
   
   void CreateLifeTimeVsStartEkinHisto(G4String aParticleName,G4String label,
			                G4int nE,G4double E1,G4double E2, 
					G4String ScaleType1,
			                G4int nTime, G4double time1,G4double time2, 
					G4String ScaleType2);

private:

  PLANETOCOSAnalysisManager* theAnalysisManager; 
  
  
  
  //Pseudo trapping  histograms
  //--------------- 
  G4bool DetectPseudoTrapping;
  
  std::vector< HISTO2D* > NbOfPlanetTurnVsStartEkinHisto;
  std::vector< G4int > NbOfPlanetTurnVsStartEkinPDGCode;
  G4bool DetectNbOfPlanetTurnVsStartEkin;
  
  std::vector< HISTO2D* > NbOfPlanetTurnVsLifeTimeHisto;
  std::vector< G4int >NbOfPlanetTurnVsLifeTimePDGCode;
  G4bool DetectNbOfPlanetTurnVsLifeTime;
  
  std::vector< HISTO2D* > NbOfEquatorCrossingVsStartEkinHisto;
  std::vector< G4int > NbOfEquatorCrossingVsStartEkinPDGCode;
  G4bool DetectNbOfEquatorCrossingVsStartEkin;
  
  std::vector< HISTO2D* > NbOfEquatorCrossingVsLifeTimeHisto;
  std::vector< G4int >NbOfEquatorCrossingVsLifeTimePDGCode;
  G4bool DetectNbOfEquatorCrossingVsLifeTime;
  
  std::vector< HISTO2D* > NbOfEquatorCrossingVsNbOfPlanetTurnHisto;
  std::vector< G4int >NbOfEquatorCrossingVsNbOfPlanetTurnPDGCode;
  G4bool DetectNbOfEquatorCrossingVsNbOfPlanetTurn;
  
  std::vector< HISTO2D* > LifeTimeVsStartEkinHisto;
  std::vector< G4int > LifeTimeVsStartEkinPDGCode;
  G4bool DetectLifeTimeVsStartEkin;
  
   
  
};
#endif
