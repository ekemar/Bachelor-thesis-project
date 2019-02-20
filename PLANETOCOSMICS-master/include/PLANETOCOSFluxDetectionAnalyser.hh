#ifndef PLANETOCOSFluxDetectionAnalyser_HH
#define PLANETOCOSFluxDetectionAnalyser_HH 1
#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"PLANETOCOSFluxHit.hh"
#include"IHistoDefinition.hh"



class G4Step;
class PLANETOCOSPrimaryGeneratorAction;
class PLANETOCOSAnalysisManager;
//class PLANETOCOSFluxHitsCollection;

class PLANETOCOSFluxDetectionAnalyser
{
public: 

   PLANETOCOSFluxDetectionAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager);
   ~PLANETOCOSFluxDetectionAnalyser();
   
   
   //Analyse
   //--------
   void Analyse(PLANETOCOSFluxHitsCollection* FluxHC);
   void RegisterLowestAltitudeNeeded(G4int PDGcode, 
   				   G4double Energy, 
				   G4double alt); 
   
   
   //Initialisation of histograms
   //----------------------------
   
   void InitialiseHistograms(G4int nb_detectors);
   
   
   //Get methods
   //-----------
   inline bool GetDetectFlux(){return DetectFlux;}
   inline void SetDetectFlux()
		{
		//DetectFlux = aVal;
		SetDetectionFlux();
		}
   inline bool GetDetectLowestAltitudeNeeded(){return DetectLowestAltitudeNeeded;}
     
   //Set methods
   //-----------
   inline void SetDetectorAltitudes(std::vector<double>* pVec) {DetectorAltitudes =pVec;}
   inline void SetDetectorDepths(std::vector<double>* pVec) {DetectorDepths =pVec;}
   inline void SetEkinMin(G4double aVal){EkinMin=aVal;}
   inline void SetEkinMax(G4double aVal){EkinMax=aVal;}
   inline void SetMinLat(G4double aVal){MinLat=aVal;}
   inline void SetMaxLat(G4double aVal){MaxLat=aVal;}
   inline void SetMinLong(G4double aVal){MinLong=aVal;}
   inline void SetMaxLong(G4double aVal){MaxLong=aVal;}
   inline void SetDivideByCosTh(bool aBool){DivideByCosTh=aBool;}
   
   
   
   // Flux detector histograms
   //----------------
   inline void SelectDetector(G4int n){
   	if (n>0 && n <= nb_flux_detector) selected_detectors[n-1] = true;
       	else G4cout<<"wrong detector number"<<std::endl; 
   }
   inline void DeselectDetector(G4int n){
     	if (n>0 && n <= nb_flux_detector) selected_detectors[n-1] = false;
       	else G4cout<<"wrong detector number"<<std::endl; 
   }
   inline void SelectAllDetectors(){
   	selected_detectors.clear();
       	selected_detectors.insert(selected_detectors.end(),(unsigned int) nb_flux_detector,true);
   }
   inline void DeselectAllDetectors(){
   	selected_detectors.clear();
       	selected_detectors.insert(selected_detectors.end(),(unsigned int) nb_flux_detector,false);
   }
   
   
   
			 	 
   void CreateDownwardFluxVsEkinHisto(G4String aParticleName,G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType); 			 
   
   void CreateEnergyVsCosZenithHisto(G4String aParticleName,
                                     G4String label,
			             G4int nE,G4double E1,G4double E2, G4String ScaleType, 
			             G4int nCosZenith, G4double cosZen1, G4double cosZen2);
   void CreateUpwardFluxVsEkinHisto(G4String aParticleName,G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType);
   void CreateCosZenithHisto(G4String aParticleName,G4String label,
			 G4int nZenith, G4double cosZen1, G4double cosZen2);
			 
   void CreateAzimuthHisto(G4String aParticleName,G4String label,
			   G4int nAzim, G4double Azim1, G4double Azim2);
   
   void CreateDownwardFluxVsPosHisto(G4String aParticleName,G4String label,
                              G4int nlon,G4double lon_min,G4double lon_max,
			      G4int nlat,G4double lat_min,G4double lat_max); 
   void CreateUpwardFluxVsPosHisto(G4String aParticleName,G4String label,
                              G4int nlon,G4double lon_min,G4double lon_max,
			      G4int nlat,G4double lat_min,G4double lat_max); 
   void CreateDownwardFluxVsPosFlatHisto(G4String aParticleName,G4String label,
                              G4int nx, G4int ny); 
   void CreateUpwardFluxVsPosFlatHisto(G4String aParticleName,G4String label,
                              G4int nx, G4int ny); 
   
   G4double FindGeoFactor(HISTOBASE* anHisto);
   void CreateLowestAltitudeNeededHisto(G4String aParticleName,G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType,
			      G4int nAlt,G4double AltitudeMin, G4double AltitudeMax); 			      
	 
  
private:
  PLANETOCOSAnalysisManager* theAnalysisManager; 
  
  
  //Flux histograms
  //--------------- 
  
  G4int nb_flux_detector;
  std::vector<bool> selected_detectors;
  
  
  std::vector<double>* DetectorAltitudes;
  std::vector<double>* DetectorDepths;
  
  
  G4String DirDifFluxBin1;
  G4String DirDifFluxBin2;
  G4String DirDifFluxBin3;
  
  G4String OmniDifFluxBin1;
  G4String OmniDifFluxBin2;
  G4String OmniDifFluxBin3;
  
  G4String IntegFluxBin1;
  G4String IntegFluxBin2;
  G4String IntegFluxBin3;
  
  
  
  G4double EkinMin, EkinMax, MinLat, MaxLat, MinLong, MaxLong;
  G4bool DivideByCosTh;
 
  G4bool DetectFlux;
  std::vector< std::vector< HISTO1D*> > DownwardFluxVsEkinHisto;
  std::vector< std::vector< G4int> > DownwardFluxVsEkinPDGCode;
  std::vector< std::vector< G4double > >  DownwardFluxVsEkinMinLat;
  std::vector< std::vector< G4double > >  DownwardFluxVsEkinMaxLat;
  std::vector< std::vector< G4double > >  DownwardFluxVsEkinMinLon;
  std::vector< std::vector< G4double > >  DownwardFluxVsEkinMaxLon;
  std::vector< std::vector< G4bool > >    DownFluxVsEkinDivByCosTh;
  G4bool DetectDownwardFluxVsEkin;
  
  std::vector< std::vector< HISTO1D*> > UpwardFluxVsEkinHisto;
  std::vector< std::vector< G4int> > UpwardFluxVsEkinPDGCode;
  std::vector< std::vector< G4double > >  UpwardFluxVsEkinMinLat;
  std::vector< std::vector< G4double > >  UpwardFluxVsEkinMaxLat;
  std::vector< std::vector< G4double > >  UpwardFluxVsEkinMinLon;
  std::vector< std::vector< G4double > >  UpwardFluxVsEkinMaxLon;
  std::vector< std::vector< G4bool > >    UpFluxVsEkinDivByCosTh;
  G4bool DetectUpwardFluxVsEkin;
  
  std::vector< std::vector< HISTO2D*> > EnergyVsCosZenithHisto;
  std::vector< std::vector< G4int> > EnergyVsCosZenithPDGCode;
  std::vector< std::vector< G4double > >  EnergyVsCosZenithMinLat;
  std::vector< std::vector< G4double > >  EnergyVsCosZenithMaxLat;
  std::vector< std::vector< G4double > >  EnergyVsCosZenithMinLon;
  std::vector< std::vector< G4double > >  EnergyVsCosZenithMaxLon;
  std::vector< std::vector< G4bool > >    EnergyVsCosZenithDivByCosTh;
  G4bool DetectEnergyVsCosZenith;
  
  
  
  
  std::vector< std::vector< HISTO1D*> > CosZenithHisto;
  std::vector< std::vector< G4int> > CosZenithPDGCode;
  std::vector< std::vector< G4double > >  CosZenithHistoMinLat;
  std::vector< std::vector< G4double > >  CosZenithHistoMaxLat;
  std::vector< std::vector< G4double > >  CosZenithHistoMinLon;
  std::vector< std::vector< G4double > >  CosZenithHistoMaxLon;
  std::vector< std::vector< G4bool > >    CosZenDivByCosTh;
  G4bool DetectCosZenith;
  
  std::vector< std::vector< HISTO1D*> > AzimuthHisto;
  std::vector< std::vector< G4int> > AzimuthPDGCode;
  std::vector< std::vector< G4double > >  AzimuthHistoMinLat;
  std::vector< std::vector< G4double > >  AzimuthHistoMaxLat;
  std::vector< std::vector< G4double > >  AzimuthHistoMinLon;
  std::vector< std::vector< G4double > >  AzimuthHistoMaxLon;
  std::vector< std::vector< G4bool > >    AzDivByCosTh;
  G4bool DetectAzimuth;
  
  
  std::vector< std::vector< HISTO2D*> > DownwardFluxVsPosHisto;
  std::vector< std::vector< G4int> > DownwardFluxVsPosPDGCode;
  std::vector< std::vector< G4double > > DownwardFluxVsPosMinE;
  std::vector< std::vector< G4double > > DownwardFluxVsPosMaxE;
  std::vector< std::vector< G4bool > >    DownFluxVsPosDivByCosTh;
  G4bool DetectDownwardFluxVsPos;
  
  std::vector< std::vector< HISTO2D*> > UpwardFluxVsPosHisto;
  std::vector< std::vector< G4int> > UpwardFluxVsPosPDGCode;
  std::vector< std::vector< G4double > > UpwardFluxVsPosMinE;
  std::vector< std::vector< G4double > > UpwardFluxVsPosMaxE;
  std::vector< std::vector< G4bool > >    UpFluxVsPosDivByCosTh;
  G4bool DetectUpwardFluxVsPos;
  
  std::vector< std::vector< HISTO2D*> > DownwardFluxVsPosFlatHisto;
  std::vector< std::vector< G4int> > DownwardFluxVsPosFlatPDGCode;
  std::vector< std::vector< G4double > > DownwardFluxVsPosFlatMinE;
  std::vector< std::vector< G4double > > DownwardFluxVsPosFlatMaxE;
  std::vector< std::vector< G4bool > >    DownFluxVsPosFlatDivByCosTh;
  G4bool DetectDownwardFluxVsPosFlat;
  
  std::vector< std::vector< HISTO2D*> > UpwardFluxVsPosFlatHisto;
  std::vector< std::vector< G4int> > UpwardFluxVsPosFlatPDGCode;
  std::vector< std::vector< G4double > > UpwardFluxVsPosFlatMinE;
  std::vector< std::vector< G4double > > UpwardFluxVsPosFlatMaxE;
  std::vector< std::vector< G4bool > >    UpFluxVsPosFlatDivByCosTh;
  G4bool DetectUpwardFluxVsPosFlat;
  
  std::vector< HISTOBASE*> HistoList;
  std::vector< G4double>   GeoFactors;
  
  
  
  
  
  std::vector< HISTO2D*> LowestAltitudeNeededHisto;
  std::vector< G4int> LowestAltitudeNeededPDGCode;
  G4bool DetectLowestAltitudeNeeded;

private:
   
  void SetDetectionFlux();
  

};
#endif
