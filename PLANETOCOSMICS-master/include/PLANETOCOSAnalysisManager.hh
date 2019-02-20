#ifndef PLANETOCOSAnalysisManager_HH
#define PLANETOCOSAnalysisManager_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"IHistoDefinition.hh"

#include"PLANETOCOSFluxDetectionAnalyser.hh"
#include"PLANETOCOSPrimaryFluxAnalyser.hh"
#include"PLANETOCOSEdepAnalyser.hh"
#include"PLANETOCOSSoilEdepAnalyser.hh"
#include"PLANETOCOSPseudoTrappingAnalyser.hh"

 class G4Step;
 class PLANETOCOSPrimaryGeneratorAction;
 class PLANETOCOSAnalysisMessenger;
 class PLANETOCOSMagneticShieldingAnalyser;
 class PLANETOCOSEventAction;
 class PLANETOCOSSteppingAction;
 class PLANETOCOSAnalysisManager
{
public:
  
  ~PLANETOCOSAnalysisManager();
   static PLANETOCOSAnalysisManager* GetInstance();
   
   
   //Event actions
  //-------------
   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);
   void finish();
    
   //Methods to reset , save, copy and add histograms and trees 
   //--------------------------------------------------------------
  
   void ResetHistograms(); 
   void SaveTree(G4String aNameFile,G4String StoreType, G4String aNameDir, G4bool normalise);
   void ResetTree();
   void PrintHistogramTree(G4String file_name, 
                           G4String path,
			   G4String tree_file="",
			   G4String store_type="");
   void InitialiseHistogramTree(G4int nb_detectors);
   
   void SetBaseFilenameRootFile(G4String filename);
   void CreateFluxRootFile();
   void FillFluxRootFile(PLANETOCOSPrimaryHitsCollection* PrimaryHC, PLANETOCOSFluxHitsCollection* FluxHC);
   void WriteFluxRoot();
   
#ifndef USE_ANALYSIS_ROOT
   void AddTreeFiles(G4String aNameFile1, G4String store_type1,G4String aNameDir1,
                    G4String aNameFile2, G4String store_type2,G4String aNameDir2);
   void AddFileToTree(G4String aNameFile, G4String store_type,G4String aNameDir);		    
   void AddTreeToFile(G4String aNameFile, G4String store_type,G4String aNameDir);
   void CopyList(ITree* aTree,std::vector<std::string> aList, std::vector<std::string> aList1);
#endif    
   
   
   
   
   //Get Methods
   //-------------
   inline PLANETOCOSMagneticShieldingAnalyser* GetMagneticShieldingAnalyser()
   							{return theMagneticShieldingAnalyser;}
   
   inline PLANETOCOSFluxDetectionAnalyser* GetFluxDetectionAnalyser()
   							{return theFluxDetectionAnalyser;}
   
   inline PLANETOCOSPrimaryFluxAnalyser* GetPrimaryFluxAnalyser()
   							{return thePrimaryFluxAnalyser;}
   
   inline PLANETOCOSEdepAnalyser* GetEdepAnalyser()
   							{return theEdepAnalyser;}
   
   
   inline PLANETOCOSSoilEdepAnalyser* GetSoilEdepAnalyser()
   							{return theSoilEdepAnalyser;}
   
   inline PLANETOCOSPseudoTrappingAnalyser* GetPseudoTrappingAnalyser()
   							{return thePseudoTrappingAnalyser;}
  
  //Set methods
  //---------------
  
   inline void SetPrimaryAction(PLANETOCOSPrimaryGeneratorAction* aPrimaryAction)
                                           		{thePrimaryAction = aPrimaryAction;}
   inline void SetSecuritySave(G4bool aBool) {security_save= aBool;}
   inline void SetNeventsSecuritySave(G4int n) {nevents_security_save = n;}
   inline void SetSecuritySaveFilename(G4String aString)
                               {security_save_filename = aString;}
   inline void SetSecuritySaveType(G4String aString)
                               {security_save_type = aString;}
  // inline void SetPrimaryIntegralFluxAboveThreshold(G4double aVal){inc_int_flux_above_thresh=aVal;}
   inline void SetPrimaryIntegralFlux(G4double aVal){inc_int_flux=aVal;}
 
   inline void SetTypeOfNormalisation(G4String aType)
                                                {type_of_normalisation = aType;}					
   inline void SwitchOffRegisterResults()
                           {RegisterResults =false;}
   inline void SetBugVerbose(G4int n){bug_verbose=n;}
   
   
   
   inline void SetCallEventActions(G4bool aBool ){CallEventActions = aBool;}
   inline void SetDetectorAltitudes(std::vector<double>* pVec) {DetectorAltitudes =pVec;}
   inline void SetDetectorDepths(std::vector<double>* pVec) {DetectorDepths =pVec;}
   inline void SetEventAction(PLANETOCOSEventAction* anEventAction){theEventAction = anEventAction;}
   inline void SetSteppingAction(PLANETOCOSSteppingAction* aSteppingAction){theSteppingAction = aSteppingAction;}
   
   //make directory
  void CreateAndMoveOnHistoDirectory(G4String dir);
   
   //create 1D histograms with linear or logaritmic scale
  HISTO1D* Create1DHisto(G4String label,G4String dir,G4String title,G4String title_x,
                               G4int nbin, G4double low, G4double up,
			       G4String type);
  //create 2D histograms with 1st dimension with linear or logaritmic scale			      
  HISTO2D* Create2DHisto(G4String label,G4String dir,G4String title,G4String title_x,G4String title_y,
                               G4int nbin1, G4double low1, G4double up1,
			       G4String type1,
			       G4int nbin2, G4double low2, G4double up2); 
  HISTO2D* Create2DHisto(G4String label,G4String dir,G4String title,G4String title_x,G4String title_y, 
                               G4int nbin1, G4double low1, G4double up1,
			       G4String type1,
			       G4int nbin2, G4double low2, G4double up2,
			       G4String type2);
  			       			       
  //fill 1D histogram with weigth divided by bin size
  void FillHistogram1D(G4double x, G4double weight,HISTO1D* anHisto);
  
  //fill 2D histogram with weigth divided by bin size
  void FillHistogram2D(G4double x,G4double y, G4double weight,HISTO2D* anHisto);


#ifdef TEST_ELEC_AT_BOUNDARY
  void RegisterElecProducedAtBoundary(G4double cos_th,G4String ProcessName, G4int n_step);
  
#endif  

#ifdef DETECT_AFTER_SOIL
  HISTO1D* gamma_after_soil_histo;
  HISTO1D* electron_after_soil_histo;
  HISTO1D* positron_after_soil_histo;
  HISTO1D* proton_after_soil_histo;
  HISTO1D* neutron_after_soil_histo;
  HISTO1D* muminus_after_soil_histo;
  HISTO1D* muplus_after_soil_histo;
  void DetectParticleAfterSoil(const G4Step* aStep);

#endif 	

#ifdef TEST_ELEC_AT_BOUNDARY
  HISTO2D* IONELEC_AT_BOUNDARY;
  HISTO2D* PHOTOELEC_AT_BOUNDARY;
  HISTO2D* COMPTONELEC_AT_BOUNDARY;
  HISTO2D* PAIRELEC_AT_BOUNDARY;
  HISTO2D* PAIRPOS_AT_BOUNDARY;
#endif


  G4int GetParticlePDGCode(G4String aParticleName);
  G4bool CheckIfSphericalGeometryWithErrorMessage();
  G4bool CheckIfSphericalGeometry();
  G4bool CheckIfFlatGeometry();
  
  void DetectNuclideProducedInAtmosphere(G4int N, G4int Z, G4double weight);


private:
  static PLANETOCOSAnalysisManager* instance;
 

private:
  PLANETOCOSAnalysisManager(); 
  
  
 
#ifndef USE_ANALYSIS_ROOT   
  //Add histograms  
  void PrintHistogram(HISTOBASE* aHistogram,G4double scaling_factor,std::fstream & FileOutput);
  bool AddHistograms(HISTOBASE* histo1, HISTOBASE* histo2);
  bool Add1DHistograms(HISTO1D* histo1, HISTO1D* histo2);
  bool Add2DHistograms(HISTO2D* histo1, HISTO2D* histo2);
#else
  void PrintHistogram(TH1* aHistogram,
   		       G4double scaling_factor,
		       std::fstream & FileOutput,
		       G4String directory,
		       std::vector<G4String> * Info =0);  
#endif   
       
  void SetBinUnitInHistoTitle(HISTOBASE* aHistogram, G4String bin1, G4String bin2, G4String bin3,
                               G4String file_type="");
  G4String DefineBinUnit(G4String bin1, G4String bin2, G4String bin3); 
  

   
 

private:
  PLANETOCOSAnalysisMessenger* theMessenger;
  PLANETOCOSEventAction* theEventAction;
  PLANETOCOSSteppingAction* theSteppingAction;
  
  
  //The different analyser
  //---------------------
  PLANETOCOSMagneticShieldingAnalyser* theMagneticShieldingAnalyser;
  PLANETOCOSFluxDetectionAnalyser* theFluxDetectionAnalyser;
  PLANETOCOSPrimaryFluxAnalyser* thePrimaryFluxAnalyser;
  PLANETOCOSEdepAnalyser* theEdepAnalyser;
  PLANETOCOSSoilEdepAnalyser* theSoilEdepAnalyser;
  PLANETOCOSPseudoTrappingAnalyser* thePseudoTrappingAnalyser;
  
  PLANETOCOSPrimaryGeneratorAction* thePrimaryAction;
  
  G4int edepCollID;
  G4int edepSoilCollID;
  G4int primaryCollID;
  G4int fluxdetectorCollID;
  G4int post_trackCollID;
  
 
  //security save 
  G4int nevents_security_save;
  G4bool security_save;
  G4String  security_save_filename;
  G4bool security_copy_in_action;
  G4String  security_save_type;

  //Random seed
  G4bool RandomSeedNeeded; 

  //Tree, Histograms
  //----------------
#ifndef USE_ANALYSIS_ROOT 
  IAnalysisFactory  *aFact;
  ITree             *theTree;
  IHistogramFactory *histFact;
#endif  
 
 
 
  //atmosphere cosmonuc histo 
  //---------------------------
 HISTO2D* theCosmonucHisto; 
 
  
  //Bin content for histogram 
  // 1 is for unscaled case
  // 2 is for scaling to primary flux
  // 3 is for scaling per pimary part
  //-------------------
  
  
  G4String EdepBin1;
  G4String EdepBin2;
  G4String EdepBin3;
  
  G4String FluxBin1;
  G4String FluxBin2;
  G4String FluxBin3;
  
  G4String CosmoBin1;
  G4String CosmoBin2;
  G4String CosmoBin3;
  
  G4String QuasiTrapped1;
  G4String QuasiTrapped2;
  G4String QuasiTrapped3;
  
  
 
    
  //variable for normalisation
  //-------------------------------
  
  G4double inc_int_flux_above_thresh;//nb of incident at the top of the atmosphere  
                                  // per atmosphere surface per sec  
				  // with energy greater than energy_threshold
  
  G4double inc_int_flux_without_cutoff;//nb of incident at the top of the atmosphere  
                                  // per atmosphere surface per sec  
				  // without considering cutoff
  G4double inc_int_flux; // integral flux of primary at the top of the atmosphere
                         // nb of particle at the top of atmosphere per cm2 pers
			 
  //G4double energy_threshold; 
  //G4double nb_of_primaries_above_tresh; //nb of primary filling ther conting condition for scaling purpose
  G4double nb_of_primaries; // nb of generated primaries
  G4double nb_of_primaries_for_scaling;
 
                                        
  G4int    nb_of_events; // 
  G4String type_of_normalisation; //none, to incident_flux, per incident particle 
  
  
 
    
  //variable for avoid registering when a particle was blocked at a boundary
  //----------------------------------------------------------------------
  G4bool RegisterResults;
  G4int bug_verbose;
  
  //nb_of_histograms
  //-----------------
  G4int nb_of_histograms;
  
  //flux detectors
  //----------------------
  G4int nb_flux_detector;
  std::vector<bool> selected_detectors;
  std::vector<double>* DetectorAltitudes;
  std::vector<double>* DetectorDepths;
  
  
  
  G4bool CallEventActions;
    
  //variables for saving of complete event
  //--------------------------------------  
  
  TFile * complete_Event_File;
  TTree * complete_Event_Tree;
  G4String base_filename;
  bool checkRootFile;
  bool check_spherical_or_cartesic;
  
  std::vector<int> * primary_PDGcode;
  std::vector<double> * primary_Energy;
  std::vector<double> * primary_zenith;
  std::vector<double> * primary_azimuth;
  
  std::vector<double> * primary_latitude;
  std::vector<double> * primary_longitude;
  std::vector<double> * primary_posY;
  std::vector<double> * primary_posX;   
  
  std::vector<int> * secondary_PDGcode;
  std::vector<double> * secondary_Energy;
  std::vector<double> * secondary_zenith;
  std::vector<double> * secondary_azimuth;
  
  std::vector<double> * secondary_latitude;
  std::vector<double> * secondary_longitude;
  std::vector<double> * secondary_posY;
  std::vector<double> * secondary_posX; 
  std::vector<double> * secondary_BoundaryDetector; 
};

#endif




