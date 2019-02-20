#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSAnalysisMessenger.hh"

//Analyser
#include "PLANETOCOSMagneticShieldingAnalyser.hh"

#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinatePlanet.hh"
#include "myfunctions.hh"
#include "G4ParticleTable.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "PLANETOCOSEventAction.hh"

#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "PlanetUnits.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSPrimaryHit.hh"
#include "PLANETOCOSEdepHit.hh"

#include "PLANETOCOSFluxHit.hh"
#include "PLANETOCOSPostTrackHit.hh"
#include "PLANETOCOSStackingAction.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "PLANETOCOSSD.hh"
#include "PlanetManager.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

PLANETOCOSAnalysisManager* PLANETOCOSAnalysisManager::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSAnalysisManager::PLANETOCOSAnalysisManager()  
{

  theMessenger = new PLANETOCOSAnalysisMessenger(this);  
  
  //initialize for root file to save complete events 	//PVD
  checkRootFile = false;				//PVD
  base_filename = "complete_event";			//PVD

#ifndef USE_ANALYSIS_ROOT  
  
  aFact=0; 
  theTree=0;
  histFact=0;  
  
 
 //Basic structure of the histogram tree
 //---------------------------------------
  aFact = AIDA_createAnalysisFactory(); 
  ITreeFactory* treeFact = aFact->createTreeFactory();
  theTree = treeFact->create();
  delete treeFact;
  histFact = aFact->createHistogramFactory( *theTree ); 

#endif  
 //define the # analysers
 //------------------------- 
  
  theMagneticShieldingAnalyser = new PLANETOCOSMagneticShieldingAnalyser();
  
 
  thePrimaryFluxAnalyser = new PLANETOCOSPrimaryFluxAnalyser(this);
  
  theEdepAnalyser = new PLANETOCOSEdepAnalyser(this);
  
  theSoilEdepAnalyser = new PLANETOCOSSoilEdepAnalyser(this);
  
  thePseudoTrappingAnalyser = new PLANETOCOSPseudoTrappingAnalyser(this);
  
  theFluxDetectionAnalyser = new PLANETOCOSFluxDetectionAnalyser(this); 
  
  //Security save
  //----------------
  nevents_security_save =10000;
  security_save = false;
#ifdef USE_ANALYSIS_ROOT
  security_save_filename ="security_save.root";
  security_save_type   ="root";
#else
  security_save_filename ="security_save.xml";
  security_save_type   ="xml";
#endif 
  security_copy_in_action =false; 
  
  //normalisation variable
  //----------------------
  inc_int_flux =0.;
  nb_of_primaries=0;
  nb_of_events=0;
  nb_of_primaries_for_scaling =  0;
  	
	
  type_of_normalisation="TO_PRIMARY_FLUX"; 
  
  //Bin content for histogram 
  // 1 is for unnormalised case
  // 2 is for normalisation to primary flux
  // 3 is for normalisation to one  pimary particle
  //-------------------
  EdepBin1 ="Edep[rad*cm2]";
  EdepBin2 ="Edep[rad/s]";
  EdepBin3 ="Edep[rad*cm2/nb primaries]";
  
  FluxBin1 ="Flux[nb particles]";
  FluxBin2 ="Flux[nb particles/cm2/s]";
  FluxBin3 ="Flux[nb particles/nb primaries]";
  
  CosmoBin1 ="Nb Produced nuclei [nb nuclei]";
  CosmoBin2 ="Nb Produced nuclei [nb nuclei/s/cm2]";
  CosmoBin3 ="Nb Produced nuclei [nb nuclei/nb primaries]";
  
  QuasiTrapped1="Nb quasi trapped particles [nb particles]"; 
  QuasiTrapped2="Nb quasi trapped particles [nb particles]"; 
  QuasiTrapped3="Nb quasi trapped particles [nb particles]";
  
  bug_verbose =1;
  
  CallEventActions = true;
  
  edepCollID =-1;
  edepSoilCollID =-1;
  primaryCollID =-1;
  fluxdetectorCollID =-1;
  post_trackCollID=-1;	

/*#ifdef DETECT_AFTER_SOIL
  #ifndef USE_ANALYSIS_ROOT 

  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/gamma/");
  theTree->cd("/AFTERSOIL/gamma/");
  gamma_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("gammas that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  
  theTree->mkdir("/AFTERSOIL/proton/");
  theTree->cd("/AFTERSOIL/proton/");
  proton_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("protons that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  
  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/neutron/");
  theTree->cd("/AFTERSOIL/neutron/");
  neutron_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("neutrons that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  
  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/e-/");
  theTree->cd("/AFTERSOIL/e-/");
  electron_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("e-s that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  
  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/e+/");
  theTree->cd("/AFTERSOIL/e+/");
  positron_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("e+s that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  
  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/mu-/");
  theTree->cd("/AFTERSOIL/mu-/");
  muminus_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("mu-s that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
  theTree->mkdir("/AFTERSOIL/");
  theTree->mkdir("/AFTERSOIL/mu+/");
  theTree->cd("/AFTERSOIL/mu+/");
  muplus_after_soil_histo = Create1DHisto(G4String("1"),
                                         G4String("mu+s that reach the planet core"),
                                         18, 1.e-3, 1.e6,G4String("log"));
#endif
#endif
*/
  theCosmonucHisto = Create2DHisto("1","/COSMONUC",
                                   "Production of nuclides in atmosphere","N","Z", 
                                   81, -0.5 , 80.5, "Lin",  81, -0.5 , 80.5);


#ifdef TEST_ELEC_AT_BOUNDARY
  IONELEC_AT_BOUNDARY =Create2DHisto("1","/TESTELEC",
                                   "Cos th of ion e- produced at boundary","cos_th", "n_step" ,
                                   40, -1. , 1., "Lin",50,-0.5,50.5);
  PHOTOELEC_AT_BOUNDARY =Create2DHisto("2","/TESTELEC",
                                   "Cos th of photoelectric e- produced at boundary","cos_th","n_step",  
                                   40, -1. , 1., "Lin",50,-0.5,50.5);
  COMPTONELEC_AT_BOUNDARY =Create2DHisto("3","/TESTELEC",
                                   "Cos th of compton e- produced at boundary","cos_th","n_step",  
                                   40, -1. , 1., "Lin", 50,-0.5,50.5);
  PAIRELEC_AT_BOUNDARY =Create2DHisto("4","/TESTELEC",
                                   "Cos th of pair e- produced at boundary","cos_th","n_step", 
                                   40, -1. , 1., "Lin",50,-0.5,50.5);
  PAIRPOS_AT_BOUNDARY =Create2DHisto("5","/TESTELEC",
                                   "Cos th of pair positron produced at boundary","cos_th","n_step" , 
                                   40, -1. , 1., "Lin", 50,-0.5,50.5);
#endif
			           

}
////////////////////////////////////////////////////////////////////////////////
//  
PLANETOCOSAnalysisManager::~PLANETOCOSAnalysisManager() 
{
#ifndef USE_ANALYSIS_ROOT 
  delete histFact;
  histFact=0;

  delete theTree;
  theTree=0;
  
  delete aFact;
  aFact = 0;
#endif  

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSAnalysisManager* PLANETOCOSAnalysisManager::GetInstance()
{
  if (instance == 0) instance = new PLANETOCOSAnalysisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::BeginOfEventAction(const G4Event*)
{if (!CallEventActions) return;
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(edepCollID<0||primaryCollID<0||fluxdetectorCollID<0 ||post_trackCollID <0 || edepSoilCollID<0){
  	edepCollID = SDman->GetCollectionID("edepCol");
	edepSoilCollID = SDman->GetCollectionID("edepSoilCol");
	primaryCollID = SDman->GetCollectionID("primaryCol");
   	fluxdetectorCollID = SDman->GetCollectionID("detCol");
	post_trackCollID = SDman->GetCollectionID("post_trackCol");
  } 
  nb_of_events++;
  RegisterResults =true; 
 //Set stop_all_particles = false in Stacking action
  const G4UserStackingAction* theStackingAction = 
     G4RunManager::GetRunManager()->GetUserStackingAction();
  G4UserStackingAction* theStackingAction1 =    
           const_cast<G4UserStackingAction*> (theStackingAction);
 (dynamic_cast<PLANETOCOSStackingAction*> (theStackingAction1))->SetStopAllParticles(false);	 
}
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// 
void PLANETOCOSAnalysisManager::EndOfEventAction(const G4Event* evt)
{ if (!CallEventActions) return;
  //G4cout<<"ANALYSE1"<<std::endl;
  theEventAction->SetVisSecondaryMode(true); 
 
  if(edepCollID<0||primaryCollID<0||fluxdetectorCollID<0 || post_trackCollID<0 || edepSoilCollID<0) return;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  PLANETOCOSPrimaryHitsCollection* PrimaryHC = 0;
  PLANETOCOSFluxHitsCollection* FluxHC = 0;
  PLANETOCOSEdepHitsCollection* EdepHC = 0;
  PLANETOCOSEdepHitsCollection* SoilEdepHC = 0;
  
  PLANETOCOSPostTrackHitsCollection* PostTrackHC = 0;
 
  //G4cout<<"ANALYSE2"<<std::endl;
  if (RegisterResults){
        
	nb_of_primaries_for_scaling =  thePrimaryAction->GetNbOfPrimariesForScaling();
  	nb_of_primaries++;
	
  	if(HCE){
    		PrimaryHC = (PLANETOCOSPrimaryHitsCollection*)(HCE->GetHC(primaryCollID));
    		FluxHC = (PLANETOCOSFluxHitsCollection*)(HCE->GetHC(fluxdetectorCollID));
    		//FluxHC = (PLANETOCOSFluxHitsCollection*)(HCE->GetHC(1));
    		//G4cout<<"ANALYSE1"<<std::endl;
		EdepHC = (PLANETOCOSEdepHitsCollection*)(HCE->GetHC(edepCollID));
		SoilEdepHC = (PLANETOCOSEdepHitsCollection*)(HCE->GetHC(edepSoilCollID));
		
    		//EdepHC = (PLANETOCOSEdepHitsCollection*)(HCE->GetHC(2));
		PostTrackHC = (PLANETOCOSPostTrackHitsCollection*)(HCE->GetHC(post_trackCollID));
		//PostTrackHC = (PLANETOCOSPostTrackHitsCollection*)(HCE->GetHC(3));
		
  	}

	//Register primaries
 	if(PrimaryHC && thePrimaryFluxAnalyser->GetDetectPrimary()){
		thePrimaryFluxAnalyser->Analyse(PrimaryHC);
	}
  	//G4cout<<"ANALYSE3"<<std::endl; 	
 
 	//registering at flux detector
 	//----------------------------
	//G4cout<<"test"<<std::endl;
  	if(FluxHC && theFluxDetectionAnalyser->GetDetectFlux()){
	       	theFluxDetectionAnalyser->Analyse(FluxHC);
  	}

	if(PrimaryHC && PrimaryHC->entries() == 1 && FluxHC&& theFluxDetectionAnalyser->GetDetectFlux()) FillFluxRootFile(PrimaryHC, FluxHC);
	
	//G4cout<<"ANALYSE3"<<std::endl; 
	//G4cout<<"test1"<<std::endl;
	//registering Edep in atmosphere
 	//----------------------------
  	if(EdepHC && theEdepAnalyser->GetDetectEdep()){
	       	theEdepAnalyser->Analyse(EdepHC);
  	}
	//registering Edep in soil
 	//----------------------------
  	if(SoilEdepHC && theSoilEdepAnalyser->GetDetectEdep()){
	       	theSoilEdepAnalyser->Analyse(SoilEdepHC);
  	}
	//G4cout<<"ANALYSE4"<<std::endl;
	//register pseudo trapping information
	//------------------------------------
	if (PostTrackHC && thePseudoTrappingAnalyser->GetDetectPseudoTrapping()){
	    	thePseudoTrappingAnalyser->Analyse(PostTrackHC);
	}
 
 	
  
	//security save
	//-----------------
  	if (security_save){
  		G4int l=nb_of_events-nevents_security_save*(nb_of_events/nevents_security_save);
   		if ( l == 0 &&  nb_of_events!=0 ){ // the tree will be saved
    			G4cout<<"Security save"<<std::endl;
     			SaveTree(security_save_filename,security_save_type,"/security",false);
    		} 
  	}

	// print nb_of events
	//-------------------  
//	G4double lowest_altitude_needed = theSteppingAction->GetLowestAltitudeNeeded();
  	if (nb_of_events<=10){
   		G4cout<<">>Events nb "<<nb_of_events<<std::endl;
	//	G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
	}	
  	else if (nb_of_events<=50 && nb_of_events - nb_of_events/10 * 10 == 0){
   		G4cout<<">>Events nb "<<nb_of_events<<std::endl;
	//	G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
	}
  	else if (nb_of_events<=500 && nb_of_events -nb_of_events/50 *50 == 0){
   		G4cout<<">>Events nb "<<nb_of_events<<std::endl;
	//	G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
	}
  	else if (nb_of_events<=1000 && nb_of_events -nb_of_events/100 *100 == 0){
   		G4cout<<">>Events nb "<<nb_of_events<<std::endl;
	//	G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
	}
  	else if (nb_of_events -nb_of_events/500 *500 == 0){
   	G4cout<<">>Events nb "<<nb_of_events<<std::endl;
	//	G4cout<<"Lowest altitude needed "<<lowest_altitude_needed/km<<std::endl;
	}
	//G4cout<<"ANALYSE5"<<std::endl;
  }
  else {
        //G4cout<<"ANALYSE6"<<std::endl;
   	if (bug_verbose >0){
		G4cout<<" A bug occured during this event"<<std::endl;
    		G4cout<<" The results of this  event will not be registered"<<std::endl;
    	}
	thePrimaryAction->SetNbOfPrimariesForScaling(int(nb_of_primaries_for_scaling));
	nb_of_events+=-1;
	nb_of_primaries+=-1;
  } 
}

///////////////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSAnalysisManager::CreateAndMoveOnHistoDirectory(G4String dir)
{ 

  G4String aNameDir =dir;
  G4String sub_dir;
  G4String file_path ="";
#ifndef USE_ANALYSIS_ROOT
  theTree->cd("/");
#else
  gDirectory->cd("/");
#endif 
  if (aNameDir !="/" || aNameDir !=""){
  	aNameDir.remove(0,1);
  	while (aNameDir.contains('/')){
		G4int n= aNameDir.first('/');
		sub_dir =aNameDir.substr(0,n);
		file_path=file_path+"/"+sub_dir;
		aNameDir.remove(0,n+1);
#ifndef USE_ANALYSIS_ROOT
		//G4cout<<"test_mk"<<std::endl;
		//G4cout<<sub_dir<<std::endl;
		//theTree->mkdir("/"+sub_dir);
		theTree->mkdir(file_path);
		//G4cout<<"File_path "<<file_path<<std::endl;
		//theTree->cd(file_path);
#else
		if (!gDirectory->Get(file_path))  gDirectory->mkdir(sub_dir);
		gDirectory->cd(file_path);	
#endif		
	}
	if (aNameDir != ""){
		
#ifndef USE_ANALYSIS_ROOT
 		//G4cout<<aNameDir<<std::endl;	
		//theTree->mkdir(aNameDir);
		file_path=file_path+"/"+aNameDir;
		//G4cout<<file_path<<std::endl;
		theTree->mkdir(file_path);
		theTree->cd(dir);
#else
		if (!gDirectory->Get(dir))  gDirectory->mkdir(aNameDir);
		gDirectory->cd(dir);
#endif
  	}	
  }
}
#ifndef USE_ANALYSIS_ROOT
#define  NEWHISTO1D_FIXBIN histFact->createHistogram1D
#define  NEWHISTO1D_USERDEFBIN(arg1,arg2,arg3,arg4) histFact->createHistogram1D(arg1,arg2,arg4)
#define  NEWHISTO2D_FIXBIN histFact->createHistogram2D
#define  NEWHISTO2D_USERDEFBIN(arg1,arg2,arg3,arg4,arg5,arg6) histFact->createHistogram2D(arg1,arg2,arg4,arg6)
#else
#define  NEWHISTO1D_FIXBIN new TH1F
#define  NEWHISTO2D_FIXBIN new TH2F
#define  NEWHISTO1D_USERDEFBIN new TH1F
#define  NEWHISTO2D_USERDEFBIN new TH2F
#endif
///////////////////////////////////////////////////////////////////////////////////////////
//
HISTO1D* PLANETOCOSAnalysisManager::Create1DHisto(G4String label,G4String dir,
                                                  G4String title,G4String title_x,
                                                  G4int nbin, G4double low, G4double up,
			                          G4String type)
{ CreateAndMoveOnHistoDirectory(dir);
  G4String TYPE = type;
  TYPE.toUpper();
  HISTO1D* theHisto;
  
  if (TYPE =="LOG"){
#ifndef USE_ANALYSIS_ROOT  
   	std::vector< G4double> bins;
#else
	double bins[1000];
#endif	
       
    	for ( int i=0; i <= nbin; i++) {
	        G4double val_bin=low * std::pow(10., i * std::log10(up/low)/nbin);
		G4double exp_10=4.-int(std::log10(val_bin));
		G4double factor =std::pow(10., exp_10);
		val_bin=int(factor*val_bin)/factor;
		;
		
#ifndef USE_ANALYSIS_ROOT         	
		bins.push_back(val_bin);
#else
		bins[i] = val_bin;
#endif		
		
     	}
	//G4cout<<"Analysis21"<<std::endl;
	theHisto =NEWHISTO1D_USERDEFBIN(label,title,nbin,bins);
	//G4cout<<"Analysis22"<<std::endl;	
  }
  else  theHisto = NEWHISTO1D_FIXBIN(label,title,nbin,low,up);
#ifndef USE_ANALYSIS_ROOT         	
	//G4cout<<title_x<<std::endl;
	//G4cout<<theHisto<<std::endl;
	theHisto->annotation().addItem("xaxis_title",title_x);
	//G4cout<<"Analysis23"<<std::endl;
#else
	theHisto->GetXaxis()->SetTitle(title_x); 
#endif	    
  //G4cout<<theHisto<<std::endl;
  //G4cout<<"Analysis3"<<std::endl;
  return theHisto;
 
}
///////////////////////////////////////////////////////////////////////////////////////////
//	
HISTO2D* PLANETOCOSAnalysisManager::Create2DHisto(G4String label,G4String dir,
                                                    G4String title,G4String title_x,G4String title_y, 
                                                    G4int nbin1, G4double low1, G4double up1,
			                            G4String type1,
			                            G4int nbin2, G4double low2, G4double up2)
{
 return Create2DHisto(label, dir,
                      title,title_x,title_y,
                      nbin1, low1, up1,type1,
		      nbin2, low2, up2,"LIN");
 
}
///////////////////////////////////////////////////////////////////////////////////////////
//	
HISTO2D* PLANETOCOSAnalysisManager::Create2DHisto(G4String label,G4String dir,
                                                    G4String title,G4String title_x,G4String title_y,
                                                    G4int nbin1, G4double low1, G4double up1,
			                            G4String type1,
			                            G4int nbin2, G4double low2, G4double up2,
						    G4String type2)
{ CreateAndMoveOnHistoDirectory(dir);
  G4String TYPE1 = type1;
  G4String TYPE2 = type2;
  TYPE1.toUpper();
  TYPE2.toUpper();
  HISTO2D* theHisto;
  if ( TYPE1 =="LOG" || TYPE2 =="LOG" ||  TYPE2=="SIN"){
#ifndef USE_ANALYSIS_ROOT  
  	std::vector< G4double> bins1;
   	std::vector< G4double> bins2;
#else
	double bins1[5000];
	double bins2[5000];
#endif
	if  (TYPE1 =="LOG") {
    		for ( int i=0; i <= nbin1; i++){
			G4double val_bin=low1 * std::pow(10., i * std::log10(up1/low1)/nbin1);
			G4double exp_10=4.-int(std::log10(val_bin));
			G4double factor =std::pow(10., exp_10);
			val_bin=int(factor*val_bin)/factor; 
#ifndef USE_ANALYSIS_ROOT 		
         		bins1.push_back(val_bin);
#else
			bins1[i]=float(val_bin);	
#endif			
		}	
	}
		
	else {
		for ( int i=0; i <= nbin1; i++){ 
         		G4double val_bin = low1 +i*(up1-low1)/nbin1;
#ifndef USE_ANALYSIS_ROOT 		
         		bins1.push_back(val_bin);
#else
			bins1[i]=float(val_bin);	
#endif				
		}	
	}		
    	if  (TYPE2 =="LOG") {
    		for ( int i=0; i <= nbin2; i++){
			G4double val_bin=low2 * std::pow(10., i * std::log10(up2/low2)/nbin2);
			G4double exp_10=4.-int(std::log10(val_bin));
			G4double factor =std::pow(10., exp_10);
			val_bin=int(factor*val_bin)/factor;
#ifndef USE_ANALYSIS_ROOT 		
         		bins2.push_back(val_bin);
#else
			bins2[i]=float(val_bin);	
#endif				
		}
	}
	else if  (TYPE2 =="SIN") { //works only if -90 degree<low2, up2<90degree
		double sin1= std::sin(low2*degree);
		double sin2= std::sin(up2*degree);
		for ( int i=0; i <= nbin2; i++){ 

			G4double val_sin = sin1 +i*(sin2-sin1)/nbin2;
			G4double val_bin= std::asin(val_sin)/degree;
			val_bin=int(1.e4*val_bin)/1.e4;
#ifndef USE_ANALYSIS_ROOT 		
         		bins2.push_back(val_bin);
#else
			bins2[i]=float(val_bin);	
#endif         		
		}
			
	}
	else {
		for ( int i=0; i <= nbin2; i++){ 
			G4double val_bin = low2 +i*(up2-low2)/nbin2;
#ifndef USE_ANALYSIS_ROOT 		
         		bins2.push_back(val_bin);
#else
			bins2[i]=float(val_bin);	
#endif         		
		}
	}
	 theHisto = NEWHISTO2D_USERDEFBIN(label,title,nbin1,bins1,nbin2,bins2);

  }     
  else  theHisto = NEWHISTO2D_FIXBIN(label,title,nbin1,low1,up1,
                                               nbin2,low2,up2);	
#ifndef USE_ANALYSIS_ROOT         	
  theHisto->annotation().addItem("xaxis_title",title_x);
  theHisto->annotation().addItem("yaxis_title",title_y);
#else
  theHisto->GetXaxis()->SetTitle(title_x);
  theHisto->GetYaxis()->SetTitle(title_y); 
#endif	 		
  
  return theHisto;
}					      
///////////////////////////////////////////////////////////////////////////////////////////
//	
void PLANETOCOSAnalysisManager::FillHistogram1D(G4double x, G4double weight,HISTO1D* anHisto)
{
//if (weight <=0 ) G4cout<< weight<<std::endl; 
#ifndef USE_ANALYSIS_ROOT
   anHisto->fill(x,weight); 
#else 
  anHisto->Fill(x,weight); 
#endif  
  return;
 
}
/////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::FillHistogram2D(G4double x, G4double y, 
                                               G4double weight,
					       HISTO2D* anHisto)
{
#ifndef USE_ANALYSIS_ROOT
 anHisto->fill(x,y,weight);
#else 
 anHisto->Fill(x,y,weight);
#endif  
 return;
}
#ifndef USE_ANALYSIS_ROOT
//////////////////////////////////////////////////////////////////////////////////////////
//
bool PLANETOCOSAnalysisManager::AddHistograms(IHistogram* histo1, IHistogram* histo2)
{ G4int nev1,nev2;
  G4double scale1, scale2;
  
  if (theTree->findPath(*(dynamic_cast<IManagedObject*> (histo1))) == ""){
   	std::stringstream astream;
   	astream<< histo1->annotation().value("nb_of_events")<<'\t'
	       << histo1->annotation().value("normalisation_factor");
    	astream>>nev1>>scale1;
    
  }
  else {
   	nev1 = nb_of_events;
    	scale1=1;
  } 
 
  if (theTree->findPath(*(dynamic_cast<IManagedObject*> (histo2))) == ""){
   	std::stringstream astream;
    	astream << histo2->annotation().value("nb_of_events" )<<'\t'
	   	<< histo2->annotation().value("normalisation_factor");
    	astream>>nev2>>scale2;
   
  }
  else {
   	nev2 = nb_of_events;
    	scale2=1;
  }  	   

  
  G4int n1 = histo1->dimension();
  G4int n2 = histo2->dimension();
 
 

  G4bool aBool=true;
  if (n1 != n2) return false;
  else if (n1 == 1){
  	HISTO1D* histo1_1D =dynamic_cast<HISTO1D*>(histo1);
   	HISTO1D* histo2_1D =dynamic_cast<HISTO1D*>(histo2);
   	histo1_1D->scale(1./scale1); 
   	histo2_1D->scale(1./scale2);
   	aBool = histo1_1D->add(*histo2_1D);
	//G4cout<<aBool<<std::endl;
   	if (!aBool) aBool =Add1DHistograms(histo1_1D,histo2_1D); 
  }
  else if (n1 == 2){
  	HISTO2D* histo1_2D =dynamic_cast<HISTO2D*>(histo1);
   	HISTO2D* histo2_2D =dynamic_cast<HISTO2D*>(histo2);
   	histo1_2D->scale(1./scale1); 
   	histo2_2D->scale(1./scale2);
   	aBool = histo1_2D->add(*histo2_2D);
	//G4cout<<aBool<<std::endl;
	if (!aBool) aBool =Add2DHistograms(histo1_2D,histo2_2D); 
  }
  histo2->scale(scale2); 
  if (!aBool){
     	histo1->scale(scale1);
	return aBool;
  }
  if (theTree->findPath(*(dynamic_cast<IManagedObject*> (histo1))) == ""){
   	std::stringstream astream;
    	astream<<nev1+nev2;
    	G4String nev;
    	astream>>nev;
  
    	histo1->annotation().setValue("nb_of_events",nev);
    	
    	if (scale1 > 0.){
     		G4double scale = scale1*nev1/(nev1+nev2);
      		histo1->scale(scale);
      		std::stringstream astream;
      		astream<<scale;
      		G4String sc;
      		astream>>sc;
      		histo1->annotation().setValue("normalisation_factor",sc);
     	}
  } 
 

 
   
  return true;
 
}
///////////////////////////////////////////////////////////////////////////////////////////
//	  					      
bool PLANETOCOSAnalysisManager::Add1DHistograms(HISTO1D* histo1, HISTO1D* histo2)
{ G4int n1=histo1->axis().bins();
  G4int n2=histo1->axis().bins();
  if (n1!= n2) return false;
  for (int i=0;i<n1;i++){
   	G4double binHeight = histo2->binHeight(i);
    	G4double x=histo2->binMean(i);
    	if (binHeight >0){
       		if (binHeight >1){
          		histo1->scale(1/binHeight);
	   		histo1->fill(x,1.);
	   		histo1->scale(binHeight);
	  	}
	else  histo1->fill(x,binHeight);  
       }
   }
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////
//
bool PLANETOCOSAnalysisManager::
                        Add2DHistograms(HISTO2D* histo1, HISTO2D* histo2)
{ G4int nx1=histo1->xAxis().bins();
  G4int ny1=histo1->yAxis().bins();
  G4int nx2=histo2->xAxis().bins();
  G4int ny2=histo2->yAxis().bins();
  
  if (nx1!= nx2 ||  ny1!= ny2) return false;
  for (int i=0;i<nx1;i++){
  	for (int j=0;j<ny1;j++){
 		G4double binHeight = histo2->binHeight(i,j);
    	 	G4double x=histo2->binMeanX(i,j);
		G4double y=histo2->binMeanY(i,j);
    		if (binHeight >0){
       			if (binHeight >1){
          			histo1->scale(1/binHeight);
	   			histo1->fill(x,y,1.);
	   			histo1->scale(binHeight);
	  		}
		else  histo1->fill(x,binHeight);  
		}
   	}
  }	
  
 return true;
}
#endif
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::ResetHistograms()
{

#ifndef USE_ANALYSIS_ROOT
 std::vector<std::string> aListOfName = theTree->listObjectNames("/",true);
 std::vector<std::string> aListOfType = theTree->listObjectTypes("/",true);
 for (unsigned i=0;i<aListOfName.size();i++){
  	if (aListOfType[i]!= "dir") {
    		IHistogram* anHisto =dynamic_cast<IHistogram*>
                         		(theTree->find("/"+aListOfName[i]));
     		anHisto->reset();			 
    	}
 }
#else
 unsigned i=0;
 std::vector<G4String> dir_vector;
 dir_vector.push_back("/");
 
 while (i <  dir_vector.size()){
 	gDirectory->cd(dir_vector[i]);
	TObjLink* lnk = gDirectory->GetList()->FirstLink();
	while (lnk){
		TObject*  object =lnk->GetObject();
		G4String class_name =G4String(object->ClassName());
		if (class_name == "TDirectory"){
			G4String path = G4String(
				(dynamic_cast<TDirectory*>(object))->GetPath());
			unsigned j=0;
			bool dir_exist =false;
			while (j<dir_vector.size()){
				if (path == dir_vector[j]) {
					j = dir_vector.size();
					dir_exist =true;
				}	
				j++;
			}
			if (!dir_exist) dir_vector.push_back(path);
		}
		else {
			(dynamic_cast<TH1*>(object))->Reset();
		}
		lnk=lnk->Next();
	}
	i++;	
	 
 }
   
#endif 
 //nb_of_primaries_above_tresh=0;
 nb_of_primaries=0;
 nb_of_events=0; 
 nb_of_primaries_for_scaling =0;  
 thePrimaryAction->ResetNbOfPrimariesForScaling();
 
 return; 
}
////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#ifndef USE_ANALYSIS_ROOT
void PLANETOCOSAnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}
#endif
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::SaveTree(G4String aNameFile, 
                                        G4String StoreType, 
                                        G4String aNameDir,
					G4bool normalize) 
{
 

  nb_of_primaries_for_scaling =  thePrimaryAction->GetNbOfPrimariesForScaling();
  G4String up_string = StoreType;
  up_string.toUpper();
  G4String low_string = StoreType;
  low_string.toLower();
  if ( up_string== "ASCII")  return PrintHistogramTree(aNameFile,G4String("/")); 
  
  //normalisation factor
  
  G4double normalisation_factor =1.;
  if (normalize) {
  	if (type_of_normalisation == "TO_PRIMARY_FLUX" && nb_of_primaries_for_scaling >0)
	  			normalisation_factor = inc_int_flux*cm2*s/nb_of_primaries_for_scaling;

  	else if (type_of_normalisation == "PER_PRIMARY" && nb_of_events>0)
         			normalisation_factor = 1./nb_of_primaries;
	  			
  }
  if (normalisation_factor <= 0) normalisation_factor = 1.;
  
  // geometry type
  G4String geometry_type = ((PLANETOCOSGeometryConstruction*)
						(G4RunManager::GetRunManager()->GetUserDetectorConstruction()))
  									->GetGeometryType();
  //Source altitude 
 
  PLANETOCOSPrimaryGeneratorAction* thePrimaryAction 
	 = (PLANETOCOSPrimaryGeneratorAction*)
	       (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4double source_altitude =thePrimaryAction->GetSourceAltitude();
  
  //Rplanet 
  G4double Rplanet = PlanetManager::GetInstance()->GetRplanet();
 
 
  G4String Bin1,Bin2,Bin3;
  
#ifndef USE_ANALYSIS_ROOT
   
  //AIDA case 
  //---------
   
  //create the new tree
  
  ITreeFactory* treeFact1 = aFact->createTreeFactory();
  ITree* aTree1 = treeFact1->create(aNameFile,StoreType,false, true);
  delete treeFact1;
  
  G4String aNameDir1 = aNameDir;
 
  G4bool aBool = aTree1->mkdir("/temporary");
 
  if (aNameDir1 != "" && aNameDir1 != "/" ){ 
         aBool = aTree1->mkdir(aNameDir1);
	 size_t l=aNameDir1.length();
 	 if (aNameDir1[l-1]=='/'){
	 	aNameDir1.remove(aNameDir1.last('/'),1);
	 } 
  }				
  if (aNameDir1 == "/") aNameDir1="";
 

  aBool = aTree1->mount("/temporary",*theTree,"/");

  std::vector<std::string> aListOfName = aTree1->listObjectNames("/temporary",true);
  std::vector<std::string> aListOfType = aTree1->listObjectTypes("/temporary",true);
 
  //G4bool test;
  G4cout<<"SaveTree1"<<std::endl;
  for (unsigned i=0;i<aListOfName.size();i++){
  
   	if (aListOfType[i]== "dir") {
	     //G4cout<<aListOfName[i]<<std::endl;
             aTree1->mkdir(aNameDir1+"/"+aListOfName[i]);
	     //G4cout<<aListOfName[i]<<std::endl;
	}     
   	else {  
		//G4cout<<aNameDir1+"/"+aListOfName[i]<<std::endl;
		aTree1->cp("/temporary/"+aListOfName[i],aNameDir1+"/"+aListOfName[i]);
      		
		IHistogram* theHisto = dynamic_cast<IHistogram*>
                                  (aTree1->find(aNameDir1+"/"+aListOfName[i]));
      		
		//G4cout<<theHisto<<std::endl;
		Bin1=FluxBin1;
		Bin2=FluxBin2;
		Bin3=FluxBin3;
		G4double geometry_factor= 1.;
		G4String dir = aListOfName[i];
		if (dir.contains("FLUX")){
			G4int n= dir.first('D');
			G4String str_i = dir.substr(n+3,1);
			std::stringstream astream;
			G4int ndet;
			astream<<str_i;
			astream>>ndet;
			ndet=ndet-1;
			if (geometry_type == "SPHERICAL" && type_of_normalisation == "TO_PRIMARY_FLUX"){
					geometry_factor= 
						(source_altitude+Rplanet)/((*DetectorAltitudes)[ndet]+Rplanet);
					geometry_factor = geometry_factor*geometry_factor;	
					
			}
			if (type_of_normalisation == "TO_PRIMARY_FLUX")
					geometry_factor *= 
						theFluxDetectionAnalyser->FindGeoFactor(theHisto);
			
			
			
			
			
		}
		else if  (dir.contains("EDEP")){
			Bin1=EdepBin1;
			Bin2=EdepBin2;
			Bin3=EdepBin3;
		}
		else if  (dir.contains("QUASI")){
			Bin1=QuasiTrapped1;
			Bin2=QuasiTrapped2;
			Bin3=QuasiTrapped3;
			geometry_factor = 1./normalisation_factor;
			
		}
     		
		std::stringstream astream;
      		G4String nb_prim, nb_ev, prim_flux, norm_fac;
      		astream<<nb_of_primaries<<'\t'
             		<<nb_of_events<<'\t'
	     		<<inc_int_flux*cm2*s*nb_of_primaries/double(nb_of_primaries_for_scaling)<<'\t'
			<<normalisation_factor*geometry_factor;
      		astream>>nb_prim>>nb_ev>>prim_flux>>norm_fac;
      		theHisto->annotation().addItem("nb_of_events",nb_ev);
		theHisto->annotation().addItem("primary_downward_flux",prim_flux);
		theHisto->annotation().addItem("nb_primaries",nb_prim);
      		theHisto->annotation().addItem("normalisation",type_of_normalisation);
		theHisto->annotation().addItem("normalisation_factor",norm_fac);
		
		theHisto->scale(normalisation_factor*geometry_factor);				
		SetBinUnitInHistoTitle(theHisto,Bin1, Bin2, Bin3,low_string);		
     	}
  }    

  
  aTree1->commit();   
  aTree1->close();
  
  if (aTree1) delete aTree1;
  aListOfName.clear();
  aListOfType.clear();
#else 

  
  //root case
  //---------

  TFile* file;
  if( security_copy_in_action) 
    		file = new TFile(aNameFile,"RECREATE","Planetocosmics results");
//else file = new TFile(aNameFile,"UPDATE","Planetocosmics results");
  else file = new TFile(aNameFile,"RECREATE","Planetocosmics results"); //PVD 060220
  
  
  //check the directory name
  
  if (aNameDir !="/" || aNameDir !=""){
  	while (aNameDir.contains('/')){
		G4int n= aNameDir.first('/');
		aNameDir.remove(n,1);
	}
	file->mkdir(aNameDir);
	aNameDir="/"+aNameDir;
  	
  }	
  
  	

  
  //General Information is given in a ntuple
  
  file->cd(aNameDir);
  G4String path = aNameDir+G4String("/INFO");
  if (!gDirectory->Get(path)) gDirectory->mkdir("INFO");
  file->cd(path);
  //G4cout<<"Try1"<<std::endl;
  TNtuple* InfoNtuple = new TNtuple("General Info","ntuple used for some general information ","nb_events:nb_primaries:norm_type");
  if (normalize) {
  	if (type_of_normalisation == "TO_PRIMARY_FLUX") InfoNtuple->Fill(nb_of_events,nb_of_primaries,2);  
  	else if (type_of_normalisation == "PER_PRIMARY") InfoNtuple->Fill(nb_of_events,nb_of_primaries,3);
  }
  else InfoNtuple->Fill(nb_of_events,nb_of_primaries,1);
  //G4cout<<"Try2"<<std::endl;
  InfoNtuple->Write();
  //G4cout<<"Try3"<<std::endl;
  if (DetectorAltitudes){
  	TNtuple* Info2Ntuple = new TNtuple("Detector altitudes and depths","ntuple used for registering the altitude and depth of detectors","nb_det:altitude:depth"); 
  	G4cout<<DetectorAltitudes<<std::endl;
  	G4cout<<DetectorAltitudes->size()<<std::endl;
  	for (unsigned i=0;i<DetectorAltitudes->size();i++){
  			/*G4cout<<"Try4"<<std::endl;
			G4cout<<((*DetectorAltitudes)[i])/km<<std::endl;
			G4cout<<((*DetectorDepths)[i])*cm2/g<<std::endl;*/
  			Info2Ntuple->Fill(i+1,((*DetectorAltitudes)[i])/km,((*DetectorDepths)[i])*cm2/g);
			//G4cout<<"Try5"<<std::endl;
  
  	}
	Info2Ntuple->Write();
  }
  
  
  

  //Copy each histogram from the permanent tree in the root file
  
  unsigned i=0;
  std::vector<G4String> dir_vector;
  dir_vector.push_back("/");
 
 
  while (i <  dir_vector.size()){
	G4String file_path =aNameDir;
	file->cd(file_path);
	G4String dir = dir_vector[i];
	G4int n= dir.first('/');
	dir.remove(0,n+1);
	G4String sub_dir;
	while (dir.contains('/')){
			n= dir.first('/');
			sub_dir =dir.substr(0,n);
			dir.remove(0,n+1);
			file_path=file_path+"/"+sub_dir;
			if (!gDirectory->Get(file_path)) gDirectory->mkdir(sub_dir);
			file->cd(file_path);	
	}
	if (dir !="") file_path=file_path+"/"+dir;
	if (!gDirectory->Get(file_path)) gDirectory->mkdir(dir);
	gROOT->cd(dir_vector[i]);
	TObjLink* lnk = gDirectory->GetList()->FirstLink();
	file->cd(file_path);
	while (lnk){
		TObject*  object =lnk->GetObject();
		G4String class_name =G4String(object->ClassName());
		if (class_name == "TDirectory"){
			G4String path = G4String(
				(dynamic_cast<TDirectory*>(object))->GetPath());
			unsigned j=0;
			bool dir_exist =false;
			while (j<dir_vector.size()){
				if (path == dir_vector[j]) {
					j = dir_vector.size();
					dir_exist =true;
				}	
				j++;
			}
			if (!dir_exist) dir_vector.push_back(path);
		
			
		}
		else {  TH1* theHisto = dynamic_cast<TH1*>(object);
			G4double geometry_factor= 1.;
			Bin1=FluxBin1;
			Bin2=FluxBin2;
			Bin3=FluxBin3;
			if (dir_vector[i].contains("FLUX")){
				G4int n= dir_vector[i].first('D');
				G4String str_i = dir_vector[i].substr(n+3,1);
				std::stringstream astream,astream1;
				G4int ndet;
				astream<<str_i;
				astream>>ndet;
				ndet=ndet-1;
				if (geometry_type == "SPHERICAL" &&
					type_of_normalisation == "TO_PRIMARY_FLUX" ){
					geometry_factor= 
						(source_altitude+Rplanet)/((*DetectorAltitudes)[ndet]+Rplanet);
					geometry_factor = geometry_factor*geometry_factor;	
				}
				if (type_of_normalisation == "TO_PRIMARY_FLUX")
					geometry_factor *= 
						theFluxDetectionAnalyser->FindGeoFactor(theHisto);
			
			}
			else if  (dir_vector[i].contains("EDEP")){
				Bin1=EdepBin1;
			 	Bin2=EdepBin2;
			 	Bin3=EdepBin3;
			}
			else if  (dir_vector[i].contains("QUASI")){
				Bin1=QuasiTrapped1;
			 	Bin2=QuasiTrapped2;
			 	Bin3=QuasiTrapped3;
				geometry_factor = 1./normalisation_factor;
			}
			if (normalisation_factor != 1.) theHisto->Scale(normalisation_factor*geometry_factor);	
			SetBinUnitInHistoTitle(theHisto,Bin1, Bin2, Bin3);
			file->cd(file_path);
			theHisto->Write();
			 
			if (normalisation_factor != 1.) theHisto->Scale(1./geometry_factor/normalisation_factor);	
			
			
		}
		lnk=lnk->Next();
	}
	i++;	
	 
 }  
 file->Close(); 
#endif  

}
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::ResetTree( ) 
{

  
  //nb_of_primaries_above_tresh=0;
  nb_of_primaries=0;
  nb_of_events=0;
  nb_of_primaries_for_scaling =0;  
  thePrimaryAction->ResetNbOfPrimariesForScaling();
  
}




/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void PLANETOCOSAnalysisManager::InitialiseHistogramTree(G4int nb_detectors) 
{ResetTree();
 theFluxDetectionAnalyser->InitialiseHistograms(nb_detectors);
}

/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PLANETOCOSAnalysisManager::SetBaseFilenameRootFile(G4String filename)
{
base_filename = filename;
}

/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PLANETOCOSAnalysisManager::CreateFluxRootFile()
{
//check for geometry type

if (CheckIfSphericalGeometry()) check_spherical_or_cartesic = true;
else check_spherical_or_cartesic = false;

char total_filename[500];
sprintf(total_filename,"%s_%d.root",base_filename.c_str(),int(time(0)));

complete_Event_File = new TFile(total_filename, "RECREATE");
complete_Event_Tree = new TTree("properties","");

complete_Event_Tree->Branch("primary_PDGcode", &primary_PDGcode);
complete_Event_Tree->Branch("primary_Energy", &primary_Energy);
complete_Event_Tree->Branch("primary_zenith", &primary_zenith);
complete_Event_Tree->Branch("primary_azimuth", &primary_azimuth);

if(check_spherical_or_cartesic)
	{
	complete_Event_Tree->Branch("primary_latitude", &primary_latitude);
	complete_Event_Tree->Branch("primary_longitude", &primary_longitude);
	}
else
	{
	complete_Event_Tree->Branch("primary_posY", &primary_posY);
	complete_Event_Tree->Branch("primary_posX", &primary_posX);
	}

complete_Event_Tree->Branch("secondary_PDGcode", &secondary_PDGcode);
complete_Event_Tree->Branch("secondary_Energy", &secondary_Energy);
complete_Event_Tree->Branch("secondary_zenith", &secondary_zenith);
complete_Event_Tree->Branch("secondary_azimuth", &secondary_azimuth);

if(check_spherical_or_cartesic)
	{
	complete_Event_Tree->Branch("secondary_latitude", &secondary_latitude);
	complete_Event_Tree->Branch("secondary_longitude", &secondary_longitude);
	}
else
	{
	complete_Event_Tree->Branch("secondary_posY", &secondary_posY);
	complete_Event_Tree->Branch("secondary_posX", &secondary_posX);
	}

complete_Event_Tree->Branch("secondary_BoundaryDetector", &secondary_BoundaryDetector);

checkRootFile = true;
}

void PLANETOCOSAnalysisManager::FillFluxRootFile(PLANETOCOSPrimaryHitsCollection* PrimaryHC, PLANETOCOSFluxHitsCollection* FluxHC)
{
if(checkRootFile)
	{
	primary_PDGcode = new std::vector<int>;
	primary_Energy = new std::vector<double>;
	primary_zenith = new std::vector<double>;
	primary_azimuth = new std::vector<double>;

	if(check_spherical_or_cartesic)
		{
		primary_latitude = new std::vector<double>;
		primary_longitude = new std::vector<double>;
		}
	else 
		{
		primary_posX = new std::vector<double>;
		primary_posY = new std::vector<double>;
		}

	secondary_PDGcode = new std::vector<int>;
	secondary_Energy = new std::vector<double>;
	secondary_zenith = new std::vector<double>;
	secondary_azimuth = new std::vector<double>;

	if(check_spherical_or_cartesic)
		{
		secondary_latitude = new std::vector<double>;
		secondary_longitude = new std::vector<double>;
		}
	else	
		{
		secondary_posX = new std::vector<double>;
		secondary_posY = new std::vector<double>;
		}
		
	secondary_BoundaryDetector = new std::vector<double>;	
		
	for(int i = 0; i < PrimaryHC->entries(); i++)
		{
		primary_Energy->push_back((*PrimaryHC)[i]->GetEnergy());
		primary_zenith->push_back((*PrimaryHC)[i]->GetTheta()); 
		primary_azimuth->push_back((*PrimaryHC)[i]->GetPhi());
		primary_PDGcode->push_back((*PrimaryHC)[i]->GetPDGCode());

		if(check_spherical_or_cartesic)
			{
			primary_latitude->push_back((*PrimaryHC)[i]->GetLatitudeOrY()/degree);
			primary_longitude->push_back((*PrimaryHC)[i]->GetLongitudeOrX()/degree);
			}
		else
			{
			primary_posY->push_back((*PrimaryHC)[i]->GetLatitudeOrY()/km);
			primary_posX->push_back((*PrimaryHC)[i]->GetLongitudeOrX()/km); 
			}
		}
	for(int i = 0; i < FluxHC->entries(); i++)
		{
		secondary_Energy->push_back((*FluxHC)[i]->GetEnergy());
		secondary_zenith->push_back((*FluxHC)[i]->GetTheta()); 
		secondary_azimuth->push_back((*FluxHC)[i]->GetPhi());
		secondary_PDGcode->push_back((*FluxHC)[i]->GetPDGCode());

		if(check_spherical_or_cartesic)
			{
			secondary_latitude->push_back((*FluxHC)[i]->GetLatitudeOrY()/degree);
			secondary_longitude->push_back((*FluxHC)[i]->GetLongitudeOrX()/degree);
			}
		else
			{
			secondary_posY->push_back((*FluxHC)[i]->GetLatitudeOrY()/km);
			secondary_posX->push_back((*FluxHC)[i]->GetLongitudeOrX()/km); 
			}
			
		secondary_BoundaryDetector->push_back((*DetectorAltitudes)[(*FluxHC)[i]->GetnBoundaryDetector()]);
		}
	complete_Event_Tree->Fill();
	
	primary_Energy->clear();
	primary_PDGcode->clear();
	primary_zenith->clear();
	primary_azimuth->clear();

	delete primary_Energy; 		primary_Energy = 0;
	delete primary_PDGcode; 	primary_PDGcode = 0;
	delete primary_zenith; 		primary_zenith = 0;
	delete primary_azimuth;		primary_azimuth = 0;

	if(check_spherical_or_cartesic)
		{
		primary_latitude->clear();
		primary_longitude->clear();
	
		delete primary_latitude;	primary_latitude = 0;
		delete primary_longitude;	primary_longitude = 0;
		}
	else
		{
		primary_posY->clear();
		primary_posX->clear();
		delete primary_posY;		primary_posY = 0;
		delete primary_posX;		primary_posX = 0;
		}	

	secondary_Energy->clear();
	secondary_PDGcode->clear();
	secondary_zenith->clear();
	secondary_azimuth->clear();
	
	delete secondary_Energy; 	secondary_Energy = 0;
	delete secondary_PDGcode; 	secondary_PDGcode = 0;
	delete secondary_zenith; 	secondary_zenith = 0;
	delete secondary_azimuth;	secondary_azimuth = 0;

	if(check_spherical_or_cartesic)
		{
		secondary_latitude->clear();
		secondary_longitude->clear();
	
		delete secondary_latitude;	secondary_latitude = 0;
		delete secondary_longitude;	secondary_longitude = 0;
		}
	else
		{
		secondary_posY->clear();
		secondary_posX->clear();
		delete secondary_posY;		secondary_posY = 0;
		delete secondary_posX;		secondary_posX = 0;
		}	

	secondary_BoundaryDetector->clear();
	delete secondary_BoundaryDetector;	secondary_BoundaryDetector = 0;
	}
}



void PLANETOCOSAnalysisManager::WriteFluxRoot()
{
if(checkRootFile)
	{
	complete_Event_File->Write();
	complete_Event_File->Close();
	}
}


/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////




#ifndef USE_ANALYSIS_ROOT
//////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::PrintHistogram(IHistogram* aHistogram,
                                              G4double scaling_factor,
					      std::fstream& File_Output)
{
  File_Output<<std::setiosflags(std::ios::scientific);
  File_Output<<std::setprecision(6);
  for (int i=0;i<aHistogram->annotation().size();i++){
  	File_Output<<aHistogram->annotation().key(i)<<" :"<<'\t'
	      <<aHistogram->annotation().value(i)<<std::endl;;
  } 
  
 //G4cout<<"Scaling factor "<<scaling_factor<<std::endl;
 // print histogram   
  G4int ndim = aHistogram->dimension();
  if (ndim ==1){
  	IHistogram1D* aHist = dynamic_cast<IHistogram1D*> (aHistogram);
   	G4int nbins = aHist->axis().bins();
   	for (int i =0;i<nbins;i++) {
    		File_Output<<std::setw(6)<<aHist->axis().binLowerEdge(i)<<'\t'<<aHist->axis().binUpperEdge(i)
                           <<'\t'<<aHist->binMean(i)
			   <<'\t'<<aHist->binHeight(i)*scaling_factor
			   <<'\t'<<aHist->binError(i)*scaling_factor<<std::endl;
    	}
  } 
 
  else if (ndim ==2){
  	IHistogram2D* aHist = dynamic_cast<IHistogram2D*> (aHistogram);
   	G4int nxbins = aHist->xAxis().bins();
   	G4int nybins = aHist->yAxis().bins();
   	
	for (int i =0;i<nxbins;i++) {
    		
		for (int j =0;j<nybins;j++) {
			File_Output<<aHist->xAxis().binLowerEdge(i)
			   <<'\t'<<aHist->xAxis().binUpperEdge(i)<<'\t';
			File_Output<<aHist->yAxis().binLowerEdge(j)<<'\t'<<aHist->yAxis().binUpperEdge(j)
       			           <<'\t'<<aHist->binHeight(i,j)*scaling_factor
				   <<'\t'<<aHist->binError(i,j)*scaling_factor<<std::endl;
   		}
    	}
  }
   
}
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::AddTreeFiles(G4String aNameFile1, G4String store_type1, G4String aNameDir1,
                                            G4String aNameFile2, G4String store_type2, G4String aNameDir2)
            
{ 
  ITreeFactory* treeFact1 = aFact->createTreeFactory();
  ITreeFactory* treeFact2 = aFact->createTreeFactory();
  
  ITree* aTree1 = treeFact1->create(aNameFile1,store_type1,false, false);
  ITree* aTree2 = treeFact2->create(aNameFile2,store_type2,false, false);
  delete treeFact1;
  delete treeFact2;
  
  if (store_type1 =="hbook") aNameDir1.toUpper();
  if (store_type2 =="hbook") aNameDir2.toUpper();
  
  if (!aTree1 || !aTree2) return;
  
  G4bool aBool =true;
  
  std::vector<std::string> aListOfName1 = aTree1->listObjectNames(aNameDir1,true);
  std::vector<std::string> aListOfType1 = aTree1->listObjectTypes(aNameDir1,true);
  std::vector<std::string> aListOfName2 = aTree2->listObjectNames(aNameDir2,true);
  std::vector<std::string> aListOfType2 = aTree2->listObjectTypes(aNameDir2,true);

  if  (aListOfName1.size() !=0 &&  aListOfName2.size() !=0){ 
  	for (unsigned int  i=0;i<aListOfName1.size();i++){
   		IHistogram* Histo1;
    		IHistogram* Histo2;
    		if (aListOfType1[i]!= "dir"){ 
    			G4String name1 = aNameDir1+"/"+aListOfName1[i];
     			G4String name2 = aNameDir2+"/"+aListOfName1[i];
     			if (store_type1 =="hbook") name1.toUpper();
     			if (store_type2 =="hbook") name2.toUpper(); 
     			Histo1=dynamic_cast< IHistogram*>(aTree1->find(name1));
     			Histo2=dynamic_cast< IHistogram*>(aTree2->find(name2));
     			if (Histo1 && Histo2) aBool = AddHistograms(Histo1, Histo2);
     			if (!aBool){ 
      				G4cout<<"A problem occurs when adding your selected trees"<<std::endl;
       				i = aListOfName1.size();
      			}
    		}
   	}
   
  }
   aTree1->commit();   
   aTree1->close();
   aTree2->commit();   
   aTree2->close();
}
//////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::AddFileToTree
                 (G4String aNameFile, G4String store_type, G4String aNameDir)            
{ 
  ITreeFactory* treeFact1 = aFact->createTreeFactory();
  ITree* aTree1 = treeFact1->create(aNameFile,store_type,true, false);
  
  delete treeFact1;

  if (!aTree1) return;
  
  G4bool aBool =true;
  if (store_type == "hbook" ) aNameDir.toUpper();
  
  std::vector<std::string> aListOfName1 = aTree1->listObjectNames(aNameDir,true);
  std::vector<std::string> aListOfType1 = aTree1->listObjectTypes(aNameDir,true);
  std::vector<std::string> aListOfName = theTree->listObjectNames("/",true);
  std::vector<std::string> aListOfType = theTree->listObjectTypes("/",true);
  IHistogram* Histo1 = 0;
 // G4cout<< aListOfName1.size()<<" "<< aListOfName.size()<<G4endl;
  if  (aListOfName1.size() !=0 &&  aListOfName.size() !=0) { 
  	for (unsigned int  i=0;i<aListOfName.size();i++){
   		if (aListOfType[i]!= "dir"){
		       
    			Histo1=dynamic_cast< IHistogram*>(aTree1->find(aNameDir+"/"+aListOfName[i]));
     			IHistogram* Histo=dynamic_cast< IHistogram*>(theTree->find("/"+aListOfName[i]));
     			if (Histo1 && Histo) aBool = AddHistograms(Histo, Histo1);
     			if (!aBool) { 
      				G4cout<<"A problem occurs when adding your selected trees"<<std::endl;
       				i = aListOfName1.size();
      			}
    		}
   	}
   
  }
  
  if (Histo1 && aBool){
  	G4int nev1;
   	std::stringstream astream;
   	astream<< Histo1->annotation().value("nb_of_events" );
   	astream>>nev1;
   	nb_of_events += nev1;
  }
 
  
  aTree1->commit();   
  aTree1->close();
}
//////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::AddTreeToFile
                 (G4String aNameFile, G4String store_type, G4String aNameDir)            
{ 
  ITreeFactory* treeFact1 = aFact->createTreeFactory();
  ITree* aTree1 = treeFact1->create(aNameFile,store_type,false, false);
  
  delete treeFact1;
  
  if (!aTree1) return;
  
  G4bool aBool = true;
  if (store_type == "hbook" ) aNameDir.toUpper();
  std::vector<std::string> aListOfName1 = aTree1->listObjectNames(aNameDir,true);
  std::vector<std::string> aListOfType1 = aTree1->listObjectTypes(aNameDir,true);
  std::vector<std::string> aListOfName = theTree->listObjectNames("/",true);
  std::vector<std::string> aListOfType = theTree->listObjectTypes("/",true);
  
  // G4cout<< aListOfName1.size() << aListOfName.size()<<G4endl;
  if  (aListOfName1.size() !=0 &&  aListOfName.size() !=0){ 
  	for (unsigned int  i=0;i<aListOfName1.size();i++){
   		if (aListOfType1[i]!= "dir"){ 
    			IHistogram* Histo1=dynamic_cast< IHistogram*>(aTree1->find(aNameDir+"/"+aListOfName1[i]));
     			IHistogram* Histo=dynamic_cast<IHistogram*>(theTree->find("/"+aListOfName1[i]));
     			if (Histo1 && Histo) aBool = AddHistograms(Histo1, Histo);
     			if (!aBool) { 
      				G4cout<<"A problem occurs when adding your selected trees"<<std::endl;
       				i = aListOfName1.size();
      			}
    		}
   	}
   }
   aTree1->commit();   
   aTree1->close();
}
//////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::CopyList(ITree* aTree,std::vector<std::string> aList, 
                                                         std::vector<std::string> aList1)
{G4bool test;
 for (unsigned i=0;i<aList.size();i++)
 {G4cout<<aList[i]<<G4endl;
  G4cout<<aList1[i]<<G4endl;
  if (aList1[i]== "dir") 
  {test=aTree->mkdir(aList[i]);}
  else 
   {test=aTree->cp("/temporary/"+aList[i],"/"+aList[i]);}
 };  
}
#else
void PLANETOCOSAnalysisManager::PrintHistogram(TH1* aHistogram,
					    G4double scaling_factor,
					    std::fstream & FileOutput,
					    G4String dir,
					    std::vector<G4String> * Info)
{// print titles
  int ndim =aHistogram->GetDimension();

  G4String ClassType;
  if (ndim == 1) ClassType = "Histogram1D";
  else ClassType = "Histogram2D";
  FileOutput<<"//////////////////////////////////////////////////"<<std::endl;
  FileOutput<<ClassType<<'\t'<<dir<<aHistogram->GetName()<<std::endl;
  FileOutput<<"//////////////////////////////////////////////////"<<std::endl;
  FileOutput<<"Title : "<<aHistogram->GetTitle()<<std::endl;
  FileOutput<<"Xaxis : "<<aHistogram->GetXaxis()->GetTitle()<<std::endl;
  if (ndim == 2)   FileOutput<<"Yaxis : "<<aHistogram->GetYaxis()->GetTitle()<<std::endl;
        	
  if (Info){
  	for (unsigned int i=0;i<Info->size();i++) {
  		FileOutput<<(*Info)[i]<<std::endl;
		
	}
  }
  
  FileOutput<<std::setiosflags(std::ios::scientific);
  FileOutput<<std::setprecision(6);
 
 // print histogram   
  
  if (ndim ==1){
  	G4int nbins = aHistogram->GetXaxis()->GetNbins();
   	for (int i =1;i<=nbins;i++) {
    		FileOutput<<std::setw(6)<<aHistogram->GetXaxis()->GetBinLowEdge(i)
			   <<'\t'<<aHistogram->GetXaxis()->GetBinUpEdge(i)
                           <<'\t'<<aHistogram->GetBinCenter(i)
			   <<'\t'<<aHistogram->GetBinContent(i)*scaling_factor
			   <<'\t'<<aHistogram->GetBinError(i)*scaling_factor<<std::endl;
    	}
  } 
 
  else if (ndim ==2){
   	G4int nxbins = aHistogram->GetXaxis()->GetNbins();
   	G4int nybins = aHistogram->GetYaxis()->GetNbins();
   	
	for (int i =1;i<=nxbins;i++) {
    		
		for (int j =1;j<=nybins;j++) {
			FileOutput<<aHistogram->GetXaxis()->GetBinLowEdge(i)
			           <<'\t'<<aHistogram->GetXaxis()->GetBinUpEdge(i)
				   <<'\t'<<aHistogram->GetYaxis()->GetBinLowEdge(j)
				   <<'\t'<<aHistogram->GetYaxis()->GetBinUpEdge(j)
				   <<'\t'<<aHistogram->GetBinContent(i,j)*scaling_factor
				   <<'\t'<<aHistogram->GetBinError(i,j)*scaling_factor<<std::endl;
   		}
    	}
  }
  
  
  
}

#endif
//////////////////////////////////////////////////////////////////////////////////////////
//   
void  PLANETOCOSAnalysisManager::SetBinUnitInHistoTitle(HISTOBASE* aHistogram, G4String bin1, G4String bin2, G4String bin3,
							G4String file_type)
{ G4String bin = DefineBinUnit(bin1, bin2, bin3);
  G4int n1 = bin.first('[');
  bin.remove(0,n1);
  if  (file_type == "hbook") {
  	G4int n1 = bin.first('[');
	bin.replace(n1,1,"(");
	n1 = bin.first(']');
	bin.replace(n1,1,")");
	
	
  } 
  
  
  G4String title;
#ifdef USE_ANALYSIS_ROOT  
  title = G4String(aHistogram->GetTitle());
#else
  title = G4String(aHistogram->title());
#endif  
  if (title.contains("[")) {
  	G4int n2 = title.first('[');
	title.remove(n2);	
  }   
  title+= " "+bin;
#ifdef USE_ANALYSIS_ROOT  
  aHistogram->SetTitle(title);
#else
  aHistogram->setTitle(title);
  aHistogram->annotation().addItem("Title",title);
#endif 
  
  
  
}
//////////////////////////////////////////////////////////////////////////////////////////
//
#ifdef TEST_ELEC_AT_BOUNDARY
void PLANETOCOSAnalysisManager::RegisterElecProducedAtBoundary(G4double cos_th,G4String ProcessName, G4int n_step)   
{  G4double w=1.;
    
   if (std::abs(cos_th)>0.1) {
   	w=1/std::abs(cos_th);
   }
   else {
   	w=10.;
   }
   //G4cout<<w<<std::endl;

  if (ProcessName.contains("Ioni")){
  	FillHistogram2D(cos_th,n_step,w,IONELEC_AT_BOUNDARY);
  }
  else if (ProcessName == "phot"){
  	FillHistogram2D(cos_th,n_step,w,PHOTOELEC_AT_BOUNDARY);
  }
  else if (ProcessName == "compt"){
  	FillHistogram2D(cos_th,n_step,w,COMPTONELEC_AT_BOUNDARY);
	
  }
  else if (ProcessName == "conve-"){
  	FillHistogram2D(cos_th,n_step,w,PAIRELEC_AT_BOUNDARY);
  }
  else if (ProcessName == "conve+"){
  	FillHistogram2D(cos_th,n_step,w,PAIRPOS_AT_BOUNDARY);
  }
  
}
#endif  
//////////////////////////////////////////////////////////////////////////////////////////
//
G4String PLANETOCOSAnalysisManager::DefineBinUnit(G4String bin1, G4String bin2, G4String bin3) 
{ if (type_of_normalisation == "TO_PRIMARY_FLUX"){
 	return bin2;
  }
  else if (type_of_normalisation == "PER_PRIMARY"){
 	return bin3;
  }
  return bin1;
}
//////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::DetectNuclideProducedInAtmosphere(G4int N, G4int Z, G4double weight)
{FillHistogram2D(N,Z, weight,theCosmonucHisto);
}




////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSAnalysisManager::PrintHistogramTree(G4String file_name, 
                                                G4String path,
						G4String tree_file,
						G4String store_type
						)
{
 G4double norm_factor=1.;
 nb_of_primaries_for_scaling = 
                     thePrimaryAction->GetNbOfPrimariesForScaling();
 
 // geometry type
 G4String geometry_type = ((PLANETOCOSGeometryConstruction*)
						(G4RunManager::GetRunManager()->GetUserDetectorConstruction()))
									->GetGeometryType();
 
 
 
 if  (type_of_normalisation == "TO_PRIMARY_FLUX"){
      		norm_factor = inc_int_flux * cm2 * s /double( nb_of_primaries_for_scaling);   
 }
 else if (type_of_normalisation == "PER_PRIMARY"){
 		norm_factor = 1. / double( nb_of_primaries);	
 }
 
 
 //Source altitude 
 
 PLANETOCOSPrimaryGeneratorAction* thePrimaryAction 
	 = (PLANETOCOSPrimaryGeneratorAction*)
	       (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 G4double source_altitude =thePrimaryAction->GetSourceAltitude();
 
 

 


#ifndef USE_ANALYSIS_ROOT
 ITree* aTree =0;
 
 G4String Bin1,Bin2,Bin3;
 
 if (tree_file != ""){
 	ITreeFactory* treeFact1 = aFact->createTreeFactory();
  	aTree = treeFact1->create(tree_file,store_type,true, false);
  	delete treeFact1;
 }
 else aTree = theTree;
 
 if (!aTree) {
  	G4cout<<"The tree file that you have specified does not exist"<<std::endl;
   	G4cout<<"Procedure interrupted"<<std::endl;
   	return;
 }  
 
 G4cout<<path<<std::endl;
 std::vector<std::string> aListOfName = aTree->listObjectNames(path,true);
 std::vector<std::string> aListOfType = aTree->listObjectTypes(path,true);
 if (aListOfName.size() <1){ 
  	G4cout<<"No histograms have been found  under the path that you have specified";
   	G4cout<<std::endl<<"Procedure interrupted"<<std::endl;
   	return;
 } 
 
 //open file
 std::fstream FileOutput(file_name, std::ios::out);
 if (aTree==theTree){
 	FileOutput<<"nb_of_events: "<<nb_of_events<<std::endl;
 	FileOutput<<"nb_of_primaries: "<<nb_of_primaries_for_scaling<<std::endl;
 
 	if (geometry_type == "EUCLIDEAN"){
 		FileOutput<<"type_of_geometry: "<<"Flat"<<std::endl;
 	}	
 	else {	
 		FileOutput<<"type_of_geometry: "<<"Spherical"<<std::endl;	
 	}
 
 
 
 	if  (type_of_normalisation == "TO_PRIMARY_FLUX"){
    		FileOutput<<"downward_primary_flux: "<<inc_int_flux* cm2 * s*nb_of_primaries/double(nb_of_primaries_for_scaling)<<std::endl;
     		FileOutput<<"normalisation_type : TO_PRIMARY_FLUX"<<std::endl;
	}
 	else if (type_of_normalisation == "PER_PRIMARY"){
    		FileOutput<<"normalisation_type : PER_PRIMARY"<<std::endl;
	}
 	else  FileOutput<<"normalisation_type : NO_NORMALISATION"<<std::endl;
 
	 
 } 
 
 
 std::string dir = path;
 char c='/';
 if (aListOfType.size() == 1 && dir == aListOfName[0]) dir=""; 
 else if (dir[dir.size()-1] != c) dir =dir+"/";
 for (unsigned int i=0; i<aListOfType.size(); i++){
   	if (aListOfType[i] != "dir"){
        	FileOutput<<"//////////////////////////////////////////////////"<<std::endl;
        	FileOutput<<aListOfType[i]<<'\t'<<dir<<aListOfName[i]<<std::endl;
        	FileOutput<<"//////////////////////////////////////////////////"<<std::endl;
        	IHistogram* aHisto = dynamic_cast<IHistogram*> 
                                  (aTree->find(dir+aListOfName[i]));
		if (aTree==theTree){
			G4double geometry_factor= 1.;
			G4String dir = aListOfName[i];
			Bin1=FluxBin1;
			Bin2=FluxBin2;
			Bin3=FluxBin3;
			if (dir.contains("FLUX")){
				G4int n= dir.first('D');
				G4String str_i = dir.substr(n+3,1);
				std::stringstream astream;
				G4int ndet;
				astream<<str_i;
				astream>>ndet;
				ndet=ndet-1;
				if (geometry_type == "SPHERICAL" && type_of_normalisation == "TO_PRIMARY_FLUX"){
					geometry_factor= 
						(source_altitude+Re)/((*DetectorAltitudes)[ndet]+Re);
					geometry_factor = geometry_factor*geometry_factor;	
				}
			
			}
			else if (dir.contains("EDEP")){
				Bin1=EdepBin1;
			 	Bin2=EdepBin2;
			 	Bin3=EdepBin3;
			}
			else if (dir.contains("QUASI")){
				Bin1=QuasiTrapped1;
			 	Bin2=QuasiTrapped2;
			 	Bin3=QuasiTrapped3;
				geometry_factor= 1/norm_factor; 
			
			}
			else if (dir.contains("COSMONUC")){
				Bin1=CosmoBin1;
			 	Bin2=CosmoBin2;
			 	Bin3=CosmoBin3;
			}	
			FileOutput<<"normalisation_factor :"<<norm_factor*geometry_factor<<std::endl;
			SetBinUnitInHistoTitle(aHisto,Bin1, Bin2, Bin3);
			PrintHistogram(aHisto,norm_factor*geometry_factor,FileOutput);
			
		}		  
        	else PrintHistogram(aHisto,1,FileOutput);
	}	
 }
 if (aTree!=theTree){ 
   	aTree->commit();
    	aTree->close();
 }
#else
 //dummy commands to avoid warning
 G4String tree_file1, store_type1, path1;
 tree_file1 = tree_file;
 store_type1 = store_type;
 path1 = path;
 
 //open file
 std::fstream FileOutput(file_name, std::ios::out);
 FileOutput<<"nb_of_events: "<<nb_of_events<<std::endl;
 FileOutput<<"nb_of_primaries: "<<nb_of_primaries<<std::endl;
 
 if (geometry_type == "EUCLIDEAN"){
 	FileOutput<<"type_of_geometry: "<<"Flat"<<std::endl;
 }	
 else {	
 	FileOutput<<"type_of_geometry: "<<"Spherical"<<std::endl;	
 }
 
 
 
 if  (type_of_normalisation == "TO_PRIMARY_FLUX"){
    		FileOutput<<"downward_primary_flux: "<<inc_int_flux* cm2 * s*nb_of_primaries/double(nb_of_primaries_for_scaling)<<std::endl;
     		FileOutput<<"normalisation_type : TO_PRIMARY_FLUX"<<std::endl;
}
 else if (type_of_normalisation == "PER_PRIMARY"){
    		FileOutput<<"normalisation_type : PER_PRIMARY"<<std::endl;
}
 else  FileOutput<<"normalisation_type : NO_NORMALISATION"<<std::endl;
 
 
 
  //Copy each histogram from the permanent tree in the ASCII file
  
  unsigned i=0;
  std::vector<G4String> dir_vector;
  dir_vector.push_back("/");
  G4String Bin1,Bin2,Bin3;
 
  while (i <  dir_vector.size()){
	G4String dir = dir_vector[i];
	G4int n= dir.first('/');
	dir.remove(0,n+1);
	G4String file_path="/"+dir+"/";
	gROOT->cd(dir_vector[i]);
	TObjLink* lnk = gDirectory->GetList()->FirstLink();
	while (lnk){
		TObject*  object =lnk->GetObject();
		G4String class_name =G4String(object->ClassName());
		if (class_name == "TDirectory"){
			G4String path = G4String(
				(dynamic_cast<TDirectory*>(object))->GetPath());
			unsigned j=0;
			bool dir_exist =false;
			while (j<dir_vector.size()){
				if (path == dir_vector[j]) {
					j = dir_vector.size();
					dir_exist =true;
				}	
				j++;
			}
			if (!dir_exist) dir_vector.push_back(path);
		
			
		}
		else {  std::vector<G4String> info;
			info.clear();
			G4double geometry_factor= 1.;
			
			if (file_path.contains("FLUX") ){
				Bin1=FluxBin1;
				Bin2=FluxBin2;
				Bin3=FluxBin3;
				G4int n= file_path.first('D');
				G4String str_i = file_path.substr(n+3,1);
				std::stringstream astream,astream1;
				G4int ndet;
				astream<<str_i;
				astream>>ndet;
				ndet=ndet-1;
				G4String str_alt,str_depth;
				astream1<<((*DetectorAltitudes)[ndet])/km<<'\t'<<((*DetectorDepths)[ndet])*cm2/g;
         			astream1>>str_alt>>str_depth;
				info.push_back("Altitude : "+str_alt+" km");
				info.push_back("Depth : "+str_depth+" g/cm2");
				if (geometry_type == "SPHERICAL" && type_of_normalisation == "TO_PRIMARY_FLUX"){
					geometry_factor= 
						(source_altitude+Re)/((*DetectorAltitudes)[ndet]+Re);
					geometry_factor = geometry_factor*geometry_factor;	
				
				}
			}
			else if (file_path.contains("PRIMARY")){
				Bin1=FluxBin1;
				Bin2=FluxBin2;
				Bin3=FluxBin3;
			}
			else if  (file_path.contains("EDEP")){
				Bin1=EdepBin1;
			 	Bin2=EdepBin2;
			 	Bin3=EdepBin3;
			}
			else if (dir.contains("QUASI")){
				Bin1=QuasiTrapped1;
			 	Bin2=QuasiTrapped2;
			 	Bin3=QuasiTrapped3;
				geometry_factor= 1/norm_factor; 
			
			}
			else if (dir.contains("COSMONUC")){
				Bin1=CosmoBin1;
			 	Bin2=CosmoBin2;
			 	Bin3=CosmoBin3;
			}
			std::stringstream astream;
			G4String str_norm_factor;
			astream<<geometry_factor*norm_factor;
			astream>>str_norm_factor;
			info.push_back("normalisation_factor : "+str_norm_factor);
				
			TH1* theHisto = dynamic_cast<TH1*>(object);
			SetBinUnitInHistoTitle(theHisto,Bin1, Bin2, Bin3);
			PrintHistogram(theHisto,geometry_factor*norm_factor,
		       			FileOutput,file_path,&info);
			
		}
		lnk=lnk->Next();
	}
	i++;	
	 
 }  
 FileOutput.close();

 
 
#endif
} 
///////////////////////////////////////////////////////////////////////////////
//
#ifdef DETECT_AFTER_SOIL
#ifndef USE_ANALYSIS_ROOT
void PLANETOCOSAnalysisManager::DetectParticleAfterSoil(const G4Step* aStep)
{ G4Track* aTrack=aStep->GetTrack();
  G4double energy =aTrack->GetKineticEnergy();
  HISTO1D* theHisto= 0;
  if (aTrack->GetDefinition()->GetParticleName() == "gamma")  theHisto= gamma_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "proton")  theHisto= proton_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "neutron")  theHisto= neutron_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "e-")  theHisto= electron_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "e+")  theHisto= positron_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "mu-")  theHisto= muminus_after_soil_histo;
  else if (aTrack->GetDefinition()->GetParticleName() == "mu+")  theHisto= muplus_after_soil_histo;
  if (theHisto){
  	G4double weight  = aTrack->GetWeight();
	FillHistogram1D(energy/MeV, weight,theHisto);
  	
  }



}
#endif
#endif
////////////////////////////////////////////////////////////////////////////////
//
G4int PLANETOCOSAnalysisManager::GetParticlePDGCode(G4String aParticleName)
{
  //Check if particle exists
  //------------------------ 
  G4int PDGCode;
  PDGCode =-999999; 
  if (aParticleName != "all"){
 	G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    	G4ParticleDefinition* aPartDefinition = 
                              theParticleTable->FindParticle(aParticleName);
    	if 	(aPartDefinition) 
		{
		PDGCode = aPartDefinition->GetPDGEncoding();	
		//G4cout<<aPartDefinition->GetPDGEncoding()<<" "<<aPartDefinition->GetParticleName()<<" "<<aPartDefinition->GetPDGMass()<<" "<<aPartDefinition->GetPDGWidth ()<<" "<<aPartDefinition->GetPDGCharge ()<<" "<<aPartDefinition->GetParticleType ()<<" "<<aPartDefinition->GetLeptonNumber ()<<" "<<aPartDefinition->GetBaryonNumber ()<<std::endl;
		}
    	else { 
     		G4cout<<"The particle "<<aParticleName<<" is not defined"<<std::endl;
      		G4cout<<"No  histogram will be created!"<<std::endl;
     	}	      
  }
  return PDGCode;
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSAnalysisManager::CheckIfSphericalGeometryWithErrorMessage()
{ G4String geometry_type=
	dynamic_cast<const PLANETOCOSGeometryConstruction*>(
			G4RunManager::GetRunManager()
				->GetUserDetectorConstruction())
						->GetGeometryType();
  if  (geometry_type !="SPHERICAL"){
 	G4cout<<"You are in the case of a flat geometry!"<<std::endl;
	G4cout<<"The type of histogram that you have selected can be created only for a spherical geometry"<<std::endl; 
	G4cout<<"No histograms will be created!"<<std::endl; 
	return false; 
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSAnalysisManager::CheckIfSphericalGeometry()
{ G4String geometry_type=
	dynamic_cast<const PLANETOCOSGeometryConstruction*>(
			G4RunManager::GetRunManager()
				->GetUserDetectorConstruction())
						->GetGeometryType();
  if  (geometry_type !="SPHERICAL"){
	return false; 
  }
  return true;
}  
////////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSAnalysisManager::CheckIfFlatGeometry()
{ G4String geometry_type=
	dynamic_cast<const PLANETOCOSGeometryConstruction*>(
			G4RunManager::GetRunManager()
				->GetUserDetectorConstruction())
						->GetGeometryType();
  if  (geometry_type =="SPHERICAL"){
 	G4cout<<"You are in the case of a spherical geometry!"<<std::endl;
	G4cout<<"The type of histogram that you have selected can be created only for a flat geometry"<<std::endl; 
	G4cout<<"No histograms will be created!"<<std::endl; 
	return false; 
  }
  return true;
}
