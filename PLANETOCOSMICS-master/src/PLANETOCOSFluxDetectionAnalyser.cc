#include "PLANETOCOSFluxDetectionAnalyser.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4EventManager.hh"
#include "PLANETOCOSEventAction.hh"
#include "PLANETOCOSFluxHit.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "G4SolidStore.hh"
#include "G4Box.hh"


////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSFluxDetectionAnalyser::PLANETOCOSFluxDetectionAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager):
theAnalysisManager(anAnalysisManager)

{     
  //Flux detection
  //--------------
  DetectFlux = false;
  DetectDownwardFluxVsEkin = false;
  DetectUpwardFluxVsEkin = false;
  DetectEnergyVsCosZenith= false;	
  DetectDownwardFluxVsPos = false;
  DetectUpwardFluxVsPos = false;	
  DetectCosZenith = false;
  DetectAzimuth = false;
  DetectLowestAltitudeNeeded = false;
  nb_flux_detector =0;
  EkinMin = 0.0001*eV;
  EkinMax = 100*GeV;
  MinLat = -90; 
  MaxLat =  90; 
  MinLong = -180; 
  MaxLong =  180;
  DivideByCosTh = true;
  
  
  
  //Bin content for histogram 
  // 1 is for unnormalised case
  // 2 is for normalisation to primary flux
  // 3 is for normalisation to one  pimary particle
  //-------------------
 
  DirDifFluxBin1 ="Flux[nb particles]";
  DirDifFluxBin2 ="Flux[nb particles/cm2/s]";
  DirDifFluxBin3 ="Flux[nb particles/nb_primary_part]";
  
  
  OmniDifFluxBin1 ="Flux[nb particles]";
  OmniDifFluxBin2 ="Flux[nb particles/cm2/s]";
  OmniDifFluxBin3 ="Flux[nb particles/nb_primary_part]";
  
  IntegFluxBin1 ="Flux[nb particles]";
  IntegFluxBin2 ="Flux[nb particles/cm2/s]";
  IntegFluxBin3 ="Flux[nb particles/nb_primary_part]";
  
  
  //create lowest altitude needed histo
  
  //G4int nAlt= 100;
  /*G4double AltitudeMin=14.995*km;
  G4double AltitudeMax=15*km;
  CreateLowestAltitudeNeededHisto(G4String("neutron"),G4String("1"),
			      15,1e-9*MeV,1e6*MeV,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax); 			      
  G4int nE = 8;
  G4double Emin =1e-2 *MeV;
  G4double Emax =1e6 *MeV;
  CreateLowestAltitudeNeededHisto(G4String("e-"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax); 
  CreateLowestAltitudeNeededHisto(G4String("e+"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax);
  CreateLowestAltitudeNeededHisto(G4String("mu-"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax);	 
  CreateLowestAltitudeNeededHisto(G4String("mu+"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax);
  CreateLowestAltitudeNeededHisto(G4String("gamma"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax);
  CreateLowestAltitudeNeededHisto(G4String("proton"),G4String("1"),
			      nE,Emin,Emax,G4String("log") ,
 			      nAlt,AltitudeMin,AltitudeMax);*/
}
////////////////////////////////////////////////////////////////////////////////
//  
PLANETOCOSFluxDetectionAnalyser::~PLANETOCOSFluxDetectionAnalyser() 
{
}
////////////////////////////////////////////////////////////////////////////////
//  
void PLANETOCOSFluxDetectionAnalyser::Analyse(PLANETOCOSFluxHitsCollection* FluxHC) 
{ int n_hit = FluxHC->entries();
  //G4cout<<"Test detection1"<<n_hit<<std::endl;
  for (int i=0;i<n_hit;i++) {
  
      	G4double weight = (*FluxHC)[i]->GetWeight();
       	G4double energy = (*FluxHC)[i]->GetEnergy();
       	G4double zenith = (*FluxHC)[i]->GetTheta();
  	G4double azimuth = (*FluxHC)[i]->GetPhi();
	G4double latitude = (*FluxHC)[i]->GetLatitudeOrY()/degree;
	G4double longitude = (*FluxHC)[i]->GetLongitudeOrX()/degree;
	G4double posY = (*FluxHC)[i]->GetLatitudeOrY()/km;
	G4double posX = (*FluxHC)[i]->GetLongitudeOrX()/km;		
			
       	G4int PDGCode = (*FluxHC)[i]->GetPDGCode();
       	G4int n =  (*FluxHC)[i]->GetnBoundaryDetector();
	//G4cout<<n<<std::endl;
       	G4double cos_zen = std::cos(zenith);
	
       	if (cos_zen < 1.e-2 && cos_zen>= 0) cos_zen=1.e-2;
       	if (cos_zen > -1.e-2 && cos_zen< 0) cos_zen= -1.e-2;
	G4double weight_DivideByCosTh  = weight/ std::abs(cos_zen);
     
       
       	//register in histograms
       			
	//G4cout<<"Test detection2"<<std::endl;		
       	if (DetectDownwardFluxVsEkin && zenith<90.*degree){
		//G4cout<<"Test detection3"<<std::endl;	
		for (unsigned  int j=0; j < DownwardFluxVsEkinPDGCode[n].size();j++){
	    		G4int code = DownwardFluxVsEkinPDGCode[n][j];
			if (code == PDGCode ) {
				G4bool detect= true;
				if (theAnalysisManager->CheckIfSphericalGeometry()){
					detect= false;
					G4double MinLat = DownwardFluxVsEkinMinLat[n][j];
					G4double MaxLat = DownwardFluxVsEkinMaxLat[n][j];
					if ( latitude >= MinLat && latitude <= MaxLat){ 
						G4double MinLon = DownwardFluxVsEkinMinLon[n][j];
						G4double MaxLon = DownwardFluxVsEkinMaxLon[n][j];
	         				if (longitude >= MinLon && longitude <= MaxLon) detect= true;
					}
				}	
				if (detect){
					G4double w=weight;
					if  (DownFluxVsEkinDivByCosTh[n][j]) w = weight_DivideByCosTh;
					theAnalysisManager->FillHistogram1D(energy/MeV,w, DownwardFluxVsEkinHisto[n][j]);
				}
			}
		}
	 }
       	if (DetectUpwardFluxVsEkin && zenith>90.*degree){
		//G4cout<<"Test detection4"<<std::endl;	
         	for (unsigned  int j=0; j < UpwardFluxVsEkinPDGCode[n].size();j++){
	    		G4int code = UpwardFluxVsEkinPDGCode[n][j];
			if (code == PDGCode ) {
				G4bool detect= true;
				if (theAnalysisManager->CheckIfSphericalGeometry()){
					detect= false;
					G4double MinLat = UpwardFluxVsEkinMinLat[n][j];
					G4double MaxLat = UpwardFluxVsEkinMaxLat[n][j];
					if ( latitude >= MinLat && latitude <= MaxLat){ 
						G4double MinLon = UpwardFluxVsEkinMinLon[n][j];
						G4double MaxLon = UpwardFluxVsEkinMaxLon[n][j];
	         				if (longitude >= MinLon && longitude <= MaxLon) detect= true;
					}
				}	
				if (detect){
					G4double w=weight;
					if  (UpFluxVsEkinDivByCosTh[n][j]) w = weight_DivideByCosTh;
					theAnalysisManager->FillHistogram1D(energy/MeV,w, UpwardFluxVsEkinHisto[n][j]);
				}
			}
		}
				
	 }
	if (DetectEnergyVsCosZenith ){
		//G4cout<<"Test detection5"<<std::endl;
         	for (unsigned  int j=0; j < EnergyVsCosZenithPDGCode[n].size();j++){
	    		G4int code = EnergyVsCosZenithPDGCode[n][j];
			if (code == PDGCode ) {
				G4bool detect= true;
				if (theAnalysisManager->CheckIfSphericalGeometry()){
					detect= false;
					G4double MinLat = EnergyVsCosZenithMinLat[n][j];
					G4double MaxLat = EnergyVsCosZenithMaxLat[n][j];
					if ( latitude >= MinLat && latitude <= MaxLat){
						G4double MinLon = EnergyVsCosZenithMinLon[n][j];
						G4double MaxLon = EnergyVsCosZenithMaxLon[n][j];
	         				if (longitude >= MinLon && longitude <= MaxLon)detect= true;
					}
				}
				if (detect){	
					G4double w=weight;
					if  (EnergyVsCosZenithDivByCosTh[n][j]) 
									w = weight_DivideByCosTh;
					theAnalysisManager->FillHistogram2D(energy/MeV,cos_zen,w, EnergyVsCosZenithHisto[n][j]);
	    			}
				
			}
		}
				
	 }		
       	if (DetectCosZenith){
		//G4cout<<"Test detection6"<<std::endl;
         	for (unsigned  int j=0; j < CosZenithPDGCode[n].size();j++){
	    		G4int code =CosZenithPDGCode[n][j];
	     		if (code == PDGCode || code == -1){
				G4bool detect= true;
				if (theAnalysisManager->CheckIfSphericalGeometry()){
					detect= false;
					G4double MinLat = CosZenithHistoMinLat[n][j];
					G4double MaxLat = CosZenithHistoMaxLat[n][j];
					if ( latitude >= MinLat && latitude <= MaxLat){   
						G4double MinLon = CosZenithHistoMinLon[n][j];
						G4double MaxLon = CosZenithHistoMaxLon[n][j];
	         				if (longitude >= MinLon && longitude <= MaxLon) detect =true;
					}
				}
				if (detect ){	
					G4double w=weight;
					if  (CosZenDivByCosTh[n][j]) 
							w = weight_DivideByCosTh;
	        			theAnalysisManager->FillHistogram1D(cos_zen,w, CosZenithHisto[n][j]);
					
				}
			}	
						
		}
            }
	if (DetectAzimuth){
		//G4cout<<"Test detection7"<<std::endl;
         	for (unsigned  int j=0; j < AzimuthPDGCode[n].size();j++){
	    		G4int code =AzimuthPDGCode[n][j];
	     		if (code == PDGCode || code == -1){
				G4bool detect= true;
				if (theAnalysisManager->CheckIfSphericalGeometry()){
					detect= false;
					G4double MinLat = AzimuthHistoMinLat[n][j];
					G4double MaxLat = AzimuthHistoMaxLat[n][j];
					if ( latitude >= MinLat && latitude <= MaxLat){   
						G4double MinLon = AzimuthHistoMinLon[n][j];
						G4double MaxLon = AzimuthHistoMaxLon[n][j];
	         				if (longitude >= MinLon && longitude <= MaxLon)detect =true;
					}
				}
				if (detect){		
					G4double w=weight;
					if  (AzDivByCosTh[n][j]) 
								w = weight_DivideByCosTh;
	        			theAnalysisManager->FillHistogram1D(azimuth/degree,weight, AzimuthHisto[n][j]);
					
				}
			}	
						
		}
   	}
	if (DetectDownwardFluxVsPos && zenith<90.*degree){
		//G4cout<<"Test detection8"<<std::endl;
		for (unsigned  int j=0; j < DownwardFluxVsPosPDGCode[n].size();j++){
	    		//G4cout<<"Flux test1"<<std::endl;
			G4int code = DownwardFluxVsPosPDGCode[n][j];
			//G4cout<<"Flux test111"<<std::endl;
			if (code == PDGCode ) {
			        //G4cout<<"Flux test1111"<<std::endl;
				G4double MinE = DownwardFluxVsPosMinE[n][j];
				G4double MaxE = DownwardFluxVsPosMaxE[n][j];
				//G4cout<<"Flux test11111"<<std::endl;
				if ( energy >= MinE && energy <= MaxE){
					G4double w=weight;
					//G4cout<<"Flux test22"<<std::endl;
					if  (DownFluxVsPosDivByCosTh[n][j]) 
						w = weight_DivideByCosTh;
					//G4cout<<"Flux test22"<<std::endl;  
					theAnalysisManager->FillHistogram2D(longitude,
							latitude,
							w, DownwardFluxVsPosHisto[n][j]);
	    				//G4cout<<"Flux test222"<<std::endl;
				}
			}
		}
		// G4cout<<"Flux test1"<<std::endl;
	 }
	if (DetectUpwardFluxVsPos && zenith>90.*degree){
		//G4cout<<"Test detection9"<<std::endl;
         	for (unsigned  int j=0; j < UpwardFluxVsPosPDGCode[n].size();j++){
	    		G4int code = UpwardFluxVsPosPDGCode[n][j];
			if (code == PDGCode ) {
				G4double MinE = UpwardFluxVsPosMinE[n][j];
				G4double MaxE = UpwardFluxVsPosMaxE[n][j];
				if ( energy >= MinE && energy <= MaxE){
					G4double w=weight;
					if  (DownFluxVsPosDivByCosTh[n][j]) 
							w = weight_DivideByCosTh;  
						 theAnalysisManager->FillHistogram2D(longitude,
							  	 latitude,
								 w, UpwardFluxVsPosHisto[n][j]);
				}
	    		
			}
		}
		//G4cout<<"Flux test3"<<std::endl;
	 }
	 if (DetectDownwardFluxVsPosFlat && zenith<90.*degree){
	 	//G4cout<<"Test detection10"<<std::endl;
		for (unsigned  int j=0; j < DownwardFluxVsPosFlatPDGCode[n].size();j++){
	    		//G4cout<<"Test detection101"<<std::endl;
			G4int code = DownwardFluxVsPosFlatPDGCode[n][j];
			if (code == PDGCode ) {
				G4double MinE = DownwardFluxVsPosFlatMinE[n][j];
				G4double MaxE = DownwardFluxVsPosFlatMaxE[n][j];
				//G4cout<<"Test detection102"<<std::endl;
				if ( energy >= MinE && energy <= MaxE){
					G4double w=weight;
					//G4cout<<"Test detection103"<<std::endl;
					if  (DownFluxVsPosFlatDivByCosTh[n][j])
						w = weight_DivideByCosTh;  
					//G4cout<<"Test detection104"<<std::endl;
					theAnalysisManager->FillHistogram2D(posX,
							posY,
							w, DownwardFluxVsPosFlatHisto[n][j]);
	    			}
			}
		}
	 }
	if (DetectUpwardFluxVsPosFlat && zenith>90.*degree){
		//G4cout<<"Test detection11"<<std::endl;
         	for (unsigned  int j=0; j < UpwardFluxVsPosFlatPDGCode[n].size();j++){
	    		//G4cout<<"Test detection12"<<std::endl;
			G4int code = UpwardFluxVsPosFlatPDGCode[n][j];
			//G4cout<<"Test detection13"<<std::endl;
			if (code == PDGCode ) {
				G4double MinE = UpwardFluxVsPosFlatMinE[n][j];
				G4double MaxE = UpwardFluxVsPosFlatMaxE[n][j];
				if ( energy >= MinE && energy <= MaxE){
					G4double w=weight;
					if  (UpFluxVsPosFlatDivByCosTh[n][j]) 
							w = weight_DivideByCosTh;  
						 theAnalysisManager->FillHistogram2D(posX,
							  	 posY,
								 w, UpwardFluxVsPosFlatHisto[n][j]);
				}
	    		
			}
		}
	 }
  }	 
}
//////////////////////////////////////////////////////////////////////////			 	 
//
void PLANETOCOSFluxDetectionAnalyser::RegisterLowestAltitudeNeeded(G4int , 
   				   				 G4double , 
				   				 G4double )
{ // G4cout<<"Test detection"<<std::endl;
  /*for (unsigned int i= 0;i<LowestAltitudeNeededHisto.size();i++){
  	
	if (PDGcode ==  LowestAltitudeNeededPDGCode[i]) {
		
		theAnalysisManager->FillHistogram2D(Energy/MeV,
						    alt/km,1., 
						    LowestAltitudeNeededHisto[i]);
	}
  
  }*/
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSFluxDetectionAnalyser::InitialiseHistograms(G4int nb_detectors)
{ //G4cout<<"Init "<<nb_detectors<<std::endl;
  DownwardFluxVsEkinHisto.clear();
  DownwardFluxVsEkinPDGCode.clear();
  DownwardFluxVsEkinMinLat.clear();
  DownwardFluxVsEkinMaxLat.clear();
  DownwardFluxVsEkinMinLon.clear();
  DownwardFluxVsEkinMaxLon.clear();
  DownFluxVsEkinDivByCosTh.clear();
 
  UpwardFluxVsEkinHisto.clear();
  UpwardFluxVsEkinPDGCode.clear();
  UpwardFluxVsEkinMinLat.clear();
  UpwardFluxVsEkinMaxLat.clear();
  UpwardFluxVsEkinMinLon.clear();
  UpwardFluxVsEkinMaxLon.clear();
  UpFluxVsEkinDivByCosTh.clear();
  
  EnergyVsCosZenithHisto.clear();
  EnergyVsCosZenithPDGCode.clear();
  EnergyVsCosZenithMinLat.clear();
  EnergyVsCosZenithMaxLat.clear();
  EnergyVsCosZenithMinLon.clear();
  EnergyVsCosZenithMaxLon.clear();
  EnergyVsCosZenithDivByCosTh.clear();
  
 
  CosZenithHisto.clear();
  CosZenithPDGCode.clear();
  CosZenithHistoMinLat.clear();
  CosZenithHistoMaxLat.clear();
  CosZenithHistoMinLon.clear();
  CosZenithHistoMaxLon.clear();
  CosZenDivByCosTh.clear();
 
  AzimuthHisto.clear();
  AzimuthPDGCode.clear();
  AzimuthHistoMinLat.clear();
  AzimuthHistoMaxLat.clear();
  AzimuthHistoMinLon.clear();
  AzimuthHistoMaxLon.clear();
  AzDivByCosTh.clear();
 
  DownwardFluxVsPosHisto.clear();
  DownwardFluxVsPosPDGCode.clear();
  DownwardFluxVsPosMinE.clear();
  DownwardFluxVsPosMaxE.clear();
  DownFluxVsPosDivByCosTh.clear();
 
  UpwardFluxVsPosHisto.clear();
  UpwardFluxVsPosPDGCode.clear();
  UpwardFluxVsPosMinE.clear();
  UpwardFluxVsPosMaxE.clear();
  UpFluxVsPosDivByCosTh.clear();
  
  
  DownwardFluxVsPosFlatHisto.clear();
  DownwardFluxVsPosFlatPDGCode.clear();
  DownwardFluxVsPosFlatMinE.clear();
  DownwardFluxVsPosFlatMaxE.clear();
  DownFluxVsPosFlatDivByCosTh.clear();
 
  UpwardFluxVsPosFlatHisto.clear();
  UpwardFluxVsPosFlatPDGCode.clear();
  UpwardFluxVsPosFlatMinE.clear();
  UpwardFluxVsPosFlatMaxE.clear();
  UpFluxVsPosFlatDivByCosTh.clear();
  
  HistoList.clear();
  GeoFactors.clear();
  
  
 

 for (int i=0;i<nb_detectors;i++){
    	DownwardFluxVsEkinHisto.push_back(std::vector< HISTO1D*>());
    	DownwardFluxVsEkinPDGCode.push_back(std::vector< G4int>());
    	DownwardFluxVsEkinMinLat.push_back(std::vector< G4double >());
    	DownwardFluxVsEkinMaxLat.push_back(std::vector< G4double >());
    	DownwardFluxVsEkinMinLon.push_back(std::vector< G4double >());
    	DownwardFluxVsEkinMaxLon.push_back(std::vector< G4double >());
   	 DownFluxVsEkinDivByCosTh.push_back(std::vector< G4bool >());
    
    	UpwardFluxVsEkinHisto.push_back(std::vector< HISTO1D*>());
    	UpwardFluxVsEkinPDGCode.push_back(std::vector< G4int>());
    	UpwardFluxVsEkinMinLat.push_back(std::vector< G4double >());
    	UpwardFluxVsEkinMaxLat.push_back(std::vector< G4double >());
    	UpwardFluxVsEkinMinLon.push_back(std::vector< G4double >());
    	UpwardFluxVsEkinMaxLon.push_back(std::vector< G4double >());
    	UpFluxVsEkinDivByCosTh.push_back(std::vector< G4bool >());
	
	EnergyVsCosZenithHisto.push_back(std::vector< HISTO2D*>());
    	EnergyVsCosZenithPDGCode.push_back(std::vector< G4int>());
    	EnergyVsCosZenithMinLat.push_back(std::vector< G4double >());
    	EnergyVsCosZenithMaxLat.push_back(std::vector< G4double >());
    	EnergyVsCosZenithMinLon.push_back(std::vector< G4double >());
    	EnergyVsCosZenithMaxLon.push_back(std::vector< G4double >());
    	EnergyVsCosZenithDivByCosTh.push_back(std::vector< G4bool >());
	
	
    
    	CosZenithHisto.push_back(std::vector< HISTO1D*>());
    	CosZenithPDGCode.push_back(std::vector< G4int>());
    	CosZenithHistoMinLat.push_back(std::vector< G4double >());
    	CosZenithHistoMaxLat.push_back(std::vector< G4double >());
    	CosZenithHistoMinLon.push_back(std::vector< G4double >());
    	CosZenithHistoMaxLon.push_back(std::vector< G4double >());
    	CosZenDivByCosTh.push_back(std::vector< G4bool >());
    
    	AzimuthHisto.push_back(std::vector< HISTO1D*>());
    	AzimuthPDGCode.push_back(std::vector< G4int>());
    	AzimuthHistoMinLat.push_back(std::vector< G4double >());
    	AzimuthHistoMaxLat.push_back(std::vector< G4double >());
    	AzimuthHistoMinLon.push_back(std::vector< G4double >());
    	AzimuthHistoMaxLon.push_back(std::vector< G4double >());
    	AzDivByCosTh.push_back(std::vector< G4bool >());
    
    	DownwardFluxVsPosHisto.push_back(std::vector< HISTO2D*>());
    	DownwardFluxVsPosPDGCode.push_back(std::vector< G4int>());
    	DownwardFluxVsPosMinE.push_back(std::vector< G4double >());
    	DownwardFluxVsPosMaxE.push_back(std::vector< G4double >());
    	DownFluxVsPosDivByCosTh.push_back(std::vector< G4bool >());
    
    	UpwardFluxVsPosHisto.push_back(std::vector< HISTO2D*>());
   	UpwardFluxVsPosPDGCode.push_back(std::vector< G4int>());
    	UpwardFluxVsPosMinE.push_back(std::vector< G4double >());
    	UpwardFluxVsPosMaxE.push_back(std::vector< G4double >());
    	UpFluxVsPosDivByCosTh.push_back(std::vector< G4bool >());
    
    	DownwardFluxVsPosFlatHisto.push_back(std::vector< HISTO2D*>());
    	DownwardFluxVsPosFlatPDGCode.push_back(std::vector< G4int>());
    	DownwardFluxVsPosFlatMinE.push_back(std::vector< G4double >());
    	DownwardFluxVsPosFlatMaxE.push_back(std::vector< G4double >());
    	DownFluxVsPosFlatDivByCosTh.push_back(std::vector< G4bool >());
    
    	UpwardFluxVsPosFlatHisto.push_back(std::vector< HISTO2D*>());
    	UpwardFluxVsPosFlatPDGCode.push_back(std::vector< G4int>());
    	UpwardFluxVsPosFlatMinE.push_back(std::vector< G4double >());
    	UpwardFluxVsPosFlatMaxE.push_back(std::vector< G4double >());
    	UpFluxVsPosFlatDivByCosTh.push_back(std::vector< G4bool >());
   
  } 
  nb_flux_detector = nb_detectors;
  selected_detectors.clear();
  selected_detectors.insert(selected_detectors.end(),
                            (unsigned int) nb_detectors,
			    false);
   SelectAllDetectors();
   /*CreateDownwardFluxVsEkinHisto("proton",
                              "1",
			      10,100.,10000., "LIN");	*/		    
}
//////////////////////////////////////////////////////////////////////////			 	 
//
void PLANETOCOSFluxDetectionAnalyser::CreateDownwardFluxVsEkinHisto(G4String aParticleName,
                              G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType)
{ //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
  
  
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
		G4String dir = "/FLUX/DET"+str_i+"/"+aParticleName;
		//G4String dir = "/TEST/"+aParticleName;
		
         	DownwardFluxVsEkinPDGCode[i].push_back(PDGCode);
	   	DownFluxVsEkinDivByCosTh[i].push_back(DivideByCosTh);
           	G4String title= "Downward flux vs Ekin of "+aParticleName;
		if (DivideByCosTh) title = "Downward flux (1/costh) vs Ekin of "+aParticleName;
           	G4double GeoFactor =1.;
		if (theAnalysisManager->CheckIfSphericalGeometry()){  
			GeoFactor = (std::sin(MaxLat*degree)-std::sin(MinLat*degree))*(MaxLong-MinLong)/720.;
		}
		GeoFactors.push_back(1/GeoFactor);	
		/*HISTO1D* theHisto =theAnalysisManager->Create1DHisto(
						       label,dir,title,"Ekin[MeV]",
	                                               nE, E1/MeV, E2/MeV,ScaleType);
		*/
		DownwardFluxVsEkinHisto[i].push_back(theAnalysisManager->Create1DHisto(
						       label,dir,title,"Ekin[MeV]",
	                                               nE, E1/MeV, E2/MeV,ScaleType));
	   	
		//G4cout<<theHisto<<std::endl;
	   	HistoList.push_back(DownwardFluxVsEkinHisto[i].back());
		DownwardFluxVsEkinMaxLat[i].push_back(MaxLat);
	   	DownwardFluxVsEkinMinLat[i].push_back(MinLat);
	   	DownwardFluxVsEkinMaxLon[i].push_back(MaxLong);
	  	DownwardFluxVsEkinMinLon[i].push_back(MinLong);
#ifndef USE_ANALYSIS_ROOT
	  	std::stringstream astream1;
	   	G4String alt,depth,max_lat,max_lon,min_lat,min_lon;
	   	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	          	<<(*DetectorDepths)[i]*cm2/g<<'\t'<<MaxLat<<'\t'
		  	<<MinLat<<'\t'<<MaxLong<<'\t'<<MinLong;
	   	astream1>>alt>>depth>>max_lat>>min_lat>>max_lon>>min_lon;
		DownwardFluxVsEkinHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	   	DownwardFluxVsEkinHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth);
	   	DownwardFluxVsEkinHisto[i].back()->annotation().addItem("MinLat",min_lat);
	   	DownwardFluxVsEkinHisto[i].back()->annotation().addItem("MaxLat",max_lat);
	   	DownwardFluxVsEkinHisto[i].back()->annotation().addItem("MinLon",min_lon);
	   	DownwardFluxVsEkinHisto[i].back()->annotation().addItem("MaxLon",max_lon);
#endif	   	

	   	DetectDownwardFluxVsEkin = true;
	}	
 }
	  

 if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
 }
 SetDetectionFlux();
}
///////////////////////////////////////////////////////////////////////////////
//		 
void PLANETOCOSFluxDetectionAnalyser::CreateUpwardFluxVsEkinHisto(G4String aParticleName,
                              G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType)
{
  
  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
  
  
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
		G4String dir = "/FLUX/DET"+str_i+"/"+aParticleName;
		//G4String dir = "/TEST/"+aParticleName;
		
         	UpwardFluxVsEkinPDGCode[i].push_back(PDGCode);
	   	UpFluxVsEkinDivByCosTh[i].push_back(DivideByCosTh);
           	G4String title= "Upward flux vs Ekin of "+aParticleName;
		if (DivideByCosTh) title = "Upward flux (1/costh) vs Ekin of "+aParticleName;
           	G4double GeoFactor =1.;
		if (theAnalysisManager->CheckIfSphericalGeometry()){  
			GeoFactor = (std::sin(MaxLat*degree)-std::sin(MinLat*degree))*(MaxLong-MinLong)/720.;
		}
		GeoFactors.push_back(1./GeoFactor);
		UpwardFluxVsEkinHisto[i].push_back(theAnalysisManager->Create1DHisto(
						       label,dir,title,"Ekin[MeV]",
	                                               nE, E1/MeV, E2/MeV,ScaleType));
	   	
	   	HistoList.push_back(UpwardFluxVsEkinHisto[i].back());
		UpwardFluxVsEkinMaxLat[i].push_back(MaxLat);
	   	UpwardFluxVsEkinMinLat[i].push_back(MinLat);
	   	UpwardFluxVsEkinMaxLon[i].push_back(MaxLong);
	  	UpwardFluxVsEkinMinLon[i].push_back(MinLong);
#ifndef USE_ANALYSIS_ROOT
	  	std::stringstream astream1;
	   	G4String alt,depth,max_lat,max_lon,min_lat,min_lon;
	   	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	          	<<(*DetectorDepths)[i]*cm2/g<<'\t'<<MaxLat<<'\t'
		  	<<MinLat<<'\t'<<MaxLong<<'\t'<<MinLong;
	   	astream1>>alt>>depth>>max_lat>>min_lat>>max_lon>>min_lon;
		UpwardFluxVsEkinHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	   	UpwardFluxVsEkinHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth);
	   	UpwardFluxVsEkinHisto[i].back()->annotation().addItem("MinLat",min_lat);
	   	UpwardFluxVsEkinHisto[i].back()->annotation().addItem("MaxLat",max_lat);
	   	UpwardFluxVsEkinHisto[i].back()->annotation().addItem("MinLon",min_lon);
	   	UpwardFluxVsEkinHisto[i].back()->annotation().addItem("MaxLon",max_lon);
#endif	   	

	   	DetectUpwardFluxVsEkin = true;
	}	
 }
	  

 if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
 }
 SetDetectionFlux();	

}
///////////////////////////////////////////////////////////////////////////////
//		 
void PLANETOCOSFluxDetectionAnalyser::CreateEnergyVsCosZenithHisto(G4String aParticleName,
                              G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType, 
			      G4int nCosZenith, G4double cosZen1, G4double cosZen2)
{
  
  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
  
  
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
		G4String dir = "/FLUX/DET"+str_i+"/"+aParticleName;
         	EnergyVsCosZenithPDGCode[i].push_back(PDGCode);
	   	EnergyVsCosZenithDivByCosTh[i].push_back(DivideByCosTh);
           	G4String title= "Ekin vs cos_theta distribution of "+aParticleName;
		if (DivideByCosTh) title = "Ekin vs cos_theta (w = 1/costh)  of "+aParticleName;
           	G4double GeoFactor =1.;
		if (theAnalysisManager->CheckIfSphericalGeometry()){  
			GeoFactor = (std::sin(MaxLat*degree)-std::sin(MinLat*degree))*(MaxLong-MinLong)/720.;
		}
		GeoFactors.push_back(1/GeoFactor);
		EnergyVsCosZenithHisto[i].push_back(theAnalysisManager->Create2DHisto(
						       label,dir,title,"Ekin[MeV]","CosTheta",
	                                               nE, E1/MeV, E2/MeV,ScaleType, 
			             		       nCosZenith, cosZen1, cosZen2));
	   	HistoList.push_back(EnergyVsCosZenithHisto[i].back());
		EnergyVsCosZenithMaxLat[i].push_back(MaxLat);
	   	EnergyVsCosZenithMinLat[i].push_back(MinLat);
	   	EnergyVsCosZenithMaxLon[i].push_back(MaxLong);
	  	EnergyVsCosZenithMinLon[i].push_back(MinLong);
#ifndef USE_ANALYSIS_ROOT
	  	std::stringstream astream1;
	   	G4String alt,depth,max_lat,max_lon,min_lat,min_lon;
	   	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	          	<<(*DetectorDepths)[i]*cm2/g<<'\t'<<MaxLat<<'\t'
		  	<<MinLat<<'\t'<<MaxLong<<'\t'<<MinLong;
	   	astream1>>alt>>depth>>max_lat>>min_lat>>max_lon>>min_lon;
		EnergyVsCosZenithHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	   	EnergyVsCosZenithHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth);
	   	EnergyVsCosZenithHisto[i].back()->annotation().addItem("MinLat",min_lat);
	   	EnergyVsCosZenithHisto[i].back()->annotation().addItem("MaxLat",max_lat);
	   	EnergyVsCosZenithHisto[i].back()->annotation().addItem("MinLon",min_lon);
	   	EnergyVsCosZenithHisto[i].back()->annotation().addItem("MaxLon",max_lon);
#endif	   	

	   	DetectEnergyVsCosZenith = true;
	}	
 }
	  

 if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
 }
 SetDetectionFlux();	

}
/////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSFluxDetectionAnalyser::CreateCosZenithHisto(G4String aParticleName,
                         G4String label, 
			 G4int nZenith, G4double cosZen1, G4double cosZen2)
{ 
  
  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;

  for (unsigned int i=0; i<selected_detectors.size() ;i++){
      	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
         	G4String dir = "/FLUX/DET"+str_i+"/"+aParticleName;
	 	CosZenithPDGCode[i].push_back(PDGCode);
	 	CosZenDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Flux vs cos_th of "+aParticleName;
		if (DivideByCosTh)  title= "Flux (1/costh) vs cos_th  of "+aParticleName;
         	G4double GeoFactor =1.;
		if (theAnalysisManager->CheckIfSphericalGeometry()){  
			GeoFactor = (std::sin(MaxLat*degree)-std::sin(MinLat*degree))*(MaxLong-MinLong)/720.;
		}
		GeoFactors.push_back(1/GeoFactor);
		CosZenithHisto[i].push_back(theAnalysisManager->Create1DHisto(
						label,dir,title,"cos_zenith",
	                                        nZenith, cosZen1, cosZen2, "lin"));
	
	 	HistoList.push_back(CosZenithHisto[i].back());
		CosZenithHistoMaxLat[i].push_back(MaxLat);
	 	CosZenithHistoMinLat[i].push_back(MinLat);
	 	CosZenithHistoMaxLon[i].push_back(MaxLong);
	 	CosZenithHistoMinLon[i].push_back(MinLong);
#ifndef USE_ANALYSIS_ROOT	 
	 	std::stringstream astream1;
	 	G4String alt,depth,max_lat,max_lon,min_lat,min_lon;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	        	<<(*DetectorDepths)[i]*cm2/g<<'\t'<<MaxLat<<'\t'
		  	<<MinLat<<'\t'<<MaxLong<<'\t'<<MinLong;
         	astream1>>alt>>depth>>max_lat>>min_lat>>max_lon>>min_lon;
	 	CosZenithHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	CosZenithHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth);
	 	CosZenithHisto[i].back()->annotation().addItem("MinLat",min_lat);
	 	CosZenithHisto[i].back()->annotation().addItem("MaxLat",max_lat);
	 	CosZenithHisto[i].back()->annotation().addItem("MinLon",min_lon);
	 	CosZenithHisto[i].back()->annotation().addItem("MaxLon",max_lon);
#endif
	 	DetectCosZenith = true;
	}
  
  }			 	
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	
}
////////////////////////////////////////////////////////////////////////////////
//	 
void PLANETOCOSFluxDetectionAnalyser::CreateAzimuthHisto(G4String aParticleName,
                         G4String label, 
			 G4int nAzim, G4double Azim1, G4double Azim2)
{ 
  
  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;


  for (unsigned int i=0; i<selected_detectors.size() ;i++){
     	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
		G4String dir = "/FLUX/DET"+str_i+"/"+aParticleName;
         	AzimuthPDGCode[i].push_back(PDGCode);
	 	AzDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Flux vs azimuth  of "+aParticleName;
		if (DivideByCosTh) title= "Flux (1/costh) vs azimuth of "+aParticleName; 
         	G4double GeoFactor =1.;
		if (theAnalysisManager->CheckIfSphericalGeometry()){  
			GeoFactor = (std::sin(MaxLat*degree)-std::sin(MinLat*degree))*(MaxLong-MinLong)/720.;
		}
		GeoFactors.push_back(1./GeoFactor);
		AzimuthHisto[i].push_back(theAnalysisManager->Create1DHisto(label,dir,title,"Azimuth[deg]",
	                                        nAzim, Azim1, Azim2, "lin"));
	 	
	 	HistoList.push_back(AzimuthHisto[i].back());
		AzimuthHistoMaxLat[i].push_back(MaxLat);
	 	AzimuthHistoMinLat[i].push_back(MinLat);
	 	AzimuthHistoMaxLon[i].push_back(MaxLong);
	 	AzimuthHistoMinLon[i].push_back(MinLong);
#ifndef USE_ANALYSIS_ROOT	 	
		std::stringstream astream1;
	 	G4String alt,depth,max_lat,max_lon,min_lat,min_lon;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	          <<(*DetectorDepths)[i]*cm2/g<<'\t'<<MaxLat<<'\t'
		  <<MinLat<<'\t'<<MaxLong<<'\t'<<MinLong;
	 	astream1>>alt>>depth>>max_lat>>min_lat>>max_lon>>min_lon;
	 	AzimuthHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	AzimuthHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth);
	 	AzimuthHisto[i].back()->annotation().addItem("MinLat",min_lat);
	 	AzimuthHisto[i].back()->annotation().addItem("MaxLat",max_lat);
	 	AzimuthHisto[i].back()->annotation().addItem("MinLon",min_lon);
	 	AzimuthHisto[i].back()->annotation().addItem("MaxLon",max_lon);  
#endif
	 	DetectAzimuth = true;
	}
  }	
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	 	

}	
////////////////////////////////////////////////////////////////////////////////
//	 
void PLANETOCOSFluxDetectionAnalyser::CreateDownwardFluxVsPosHisto(G4String aParticleName,G4String label,
                              G4int nlon,G4double lon_min,G4double lon_max,
			      G4int nlat,G4double lat_min,G4double lat_max)
{ //Check the geometry
  //------------------
  if (!theAnalysisManager->CheckIfSphericalGeometryWithErrorMessage()) return; 

  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
  
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
         	G4String dir ="/FLUX/DET"+str_i+"/"+aParticleName;
         	DownwardFluxVsPosPDGCode[i].push_back(PDGCode);
	 	DownFluxVsPosDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Downward flux vs position of "+aParticleName;
		if (DivideByCosTh) title= "Downward flux (1/costh) vs position of "+aParticleName;
		G4double GeoFactor = (std::sin(lat_max*degree)-std::sin(lat_min*degree))*(lon_max-lon_min)/720./nlon/nlat;
		GeoFactors.push_back(1./GeoFactor);
		DownwardFluxVsPosHisto[i].push_back(theAnalysisManager->Create2DHisto(
								label,dir,title,"Longitude[deg]","Latitude[deg]",
	                                        		nlon, lon_min, lon_max, "lin",
								nlat, lat_min, lat_max,"SIN"));
	 
	 	HistoList.push_back(DownwardFluxVsPosHisto[i].back());
		DownwardFluxVsPosMaxE[i].push_back(EkinMax);
	 	DownwardFluxVsPosMinE[i].push_back(EkinMin);
	 
#ifndef USE_ANALYSIS_ROOT	 
	 	std::stringstream astream1;
	 	G4String alt,depth,max_e,min_e;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	         	<<(*DetectorDepths)[i]*cm2/g<<'\t'
		 	<<EkinMax/MeV<<'\t'<<EkinMin/MeV;
	 	astream1>>alt>>depth>>max_e>>min_e;
	 
	 	DownwardFluxVsPosHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	DownwardFluxVsPosHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth); 
	 	DownwardFluxVsPosHisto[i].back()->annotation().addItem("EkinMin [MeV]",min_e);
	 	DownwardFluxVsPosHisto[i].back()->annotation().addItem("EkinMax [MeV]",max_e);
#endif	   
	 	DetectDownwardFluxVsPos = true;
	}
  }	
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	 

}
////////////////////////////////////////////////////////////////////////////////
//	 
void PLANETOCOSFluxDetectionAnalyser::CreateUpwardFluxVsPosHisto(G4String aParticleName,G4String label,
                              G4int nlon,G4double lon_min,G4double lon_max,
			      G4int nlat,G4double lat_min,G4double lat_max)
{ //Check the geometry
  //------------------
  if (!theAnalysisManager->CheckIfSphericalGeometryWithErrorMessage()) return; 

  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detectors have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
  
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
        	ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
         	astream>>str_i;
         	G4String dir ="/FLUX/DET"+str_i+"/"+aParticleName;
         	UpwardFluxVsPosPDGCode[i].push_back(PDGCode);
	 	DownFluxVsPosDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Upward flux vs position of "+aParticleName;
		if (DivideByCosTh) title= "Upward flux (1/costh) vs position of "+aParticleName;
         	G4double GeoFactor = (std::sin(lat_max*degree)-std::sin(lat_min*degree))*(lon_max-lon_min)/720./nlon/nlat;
		GeoFactors.push_back(1./GeoFactor);
		UpwardFluxVsPosHisto[i].push_back(theAnalysisManager->Create2DHisto(
								label,dir,title,"Longitude[deg]","Latitude[deg]",
	                                        		nlon, lon_min, lon_max, "lin",
								nlat, lat_min, lat_max,"SIN"));
	 
	 	HistoList.push_back(UpwardFluxVsPosHisto[i].back());
		UpwardFluxVsPosMaxE[i].push_back(EkinMax);
	 	UpwardFluxVsPosMinE[i].push_back(EkinMin);
	 
#ifndef USE_ANALYSIS_ROOT	 
	 	std::stringstream astream1;
	 	G4String alt,depth,max_e,min_e;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	         	<<(*DetectorDepths)[i]*cm2/g<<'\t'
		 	<<EkinMax/MeV<<'\t'<<EkinMin/MeV;
	 	astream1>>alt>>depth>>max_e>>min_e;
	 
	 	UpwardFluxVsPosHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	UpwardFluxVsPosHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth); 
	 	UpwardFluxVsPosHisto[i].back()->annotation().addItem("EkinMin [MeV]",min_e);
	 	UpwardFluxVsPosHisto[i].back()->annotation().addItem("EkinMax [MeV]",max_e);
#endif	   
	 	DetectUpwardFluxVsPos = true;
	}
  }	
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	 

}
////////////////////////////////////////////////////////////////////////////////
//	 
void PLANETOCOSFluxDetectionAnalyser::CreateDownwardFluxVsPosFlatHisto(G4String aParticleName,G4String label,
                              G4int nx, G4int ny)
{ //Check the geometry
  //------------------
  if (!theAnalysisManager->CheckIfFlatGeometry()) return; 
  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detector have been defined"<<std::endl;
    	return;
  }  
  G4int ndet =0;
        
  G4SolidStore* theSolidStore = G4SolidStore::GetInstance();
  G4Box* aSolidBox = dynamic_cast<G4Box*>((*theSolidStore)[0]);
  G4double xmin, ymin, xmax, ymax; 
  xmin=-aSolidBox->GetXHalfLength()/km;
  ymin=-aSolidBox->GetYHalfLength()/km;
  xmax=aSolidBox->GetXHalfLength()/km;
  ymax=aSolidBox->GetYHalfLength()/km;
	
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
		ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
		astream>>str_i;
		G4String dir ="/FLUX/DET"+str_i+"/"+aParticleName;
        	DownwardFluxVsPosFlatPDGCode[i].push_back(PDGCode); 
	 	DownFluxVsPosFlatDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Downward flux vs XY position  of "+aParticleName;
		if (DivideByCosTh) title = "Downward flux (w=1/costh) vs XY position of "+aParticleName;
		GeoFactors.push_back(nx*ny);
		DownwardFluxVsPosFlatHisto[i].push_back(theAnalysisManager->Create2DHisto(label,dir,title,"X[km]","Y[km]",
	                                        			nx, xmin, ymax, "lin",
									ny, ymin, ymax));
	 	
	 	HistoList.push_back(DownwardFluxVsPosFlatHisto[i].back());
		DownwardFluxVsPosFlatMaxE[i].push_back(EkinMax);
	 	DownwardFluxVsPosFlatMinE[i].push_back(EkinMin);
	
#ifndef USE_ANALYSIS_ROOT
	 	std::stringstream astream1;
	 	G4String alt,depth,max_e,min_e;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	         	<<(*DetectorDepths)[i]*cm2/g<<'\t'
		 	<<EkinMax/MeV<<'\t'<<EkinMin/MeV;
	 	astream1>>alt>>depth>>max_e>>min_e;
	 
	 	DownwardFluxVsPosFlatHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	DownwardFluxVsPosFlatHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth); 
	 	DownwardFluxVsPosFlatHisto[i].back()->annotation().addItem("EkinMin [MeV]",min_e);
	 	DownwardFluxVsPosFlatHisto[i].back()->annotation().addItem("EkinMax [MeV]",max_e);
#endif   
	 	DetectDownwardFluxVsPosFlat = true;
	}
  }
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	
}
////////////////////////////////////////////////////////////////////////////////
//	 
void PLANETOCOSFluxDetectionAnalyser::CreateUpwardFluxVsPosFlatHisto(G4String aParticleName,G4String label,
                              G4int nx, G4int ny)
{ //Check the geometry
  //------------------
  if (!theAnalysisManager->CheckIfFlatGeometry()) return; 

  //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
   
  //create histogram
  //---------------
  if (nb_flux_detector == 0) {
   	G4cout<<"No flux detector have been defined"<<std::endl;
    	return;
  }  
  DetectFlux =true;
  G4int ndet =0;
        
  G4SolidStore* theSolidStore = G4SolidStore::GetInstance();
  G4Box* aSolidBox = dynamic_cast<G4Box*>((*theSolidStore)[0]);
  G4double xmin, ymin, xmax, ymax; 
  xmin=-aSolidBox->GetXHalfLength()/km;
  ymin=-aSolidBox->GetYHalfLength()/km;
  xmax=aSolidBox->GetXHalfLength()/km;
  ymax=aSolidBox->GetYHalfLength()/km;
	
  for (unsigned int i=0; i<selected_detectors.size() ;i++){
  	if (selected_detectors[i]){
		ndet++;
         	std::stringstream astream;
         	G4String str_i;
         	astream<<i+1;
		astream>>str_i;
		G4String dir ="/FLUX/DET"+str_i+"/"+aParticleName;
        	UpwardFluxVsPosFlatPDGCode[i].push_back(PDGCode); 
	 	UpFluxVsPosFlatDivByCosTh[i].push_back(DivideByCosTh);
         	G4String title= "Upward flux vs XY position  of "+aParticleName;
		if (DivideByCosTh) title = "Upward flux (w=1/costh) vs XY position of "+aParticleName;
		GeoFactors.push_back(nx*ny);
		UpwardFluxVsPosFlatHisto[i].push_back(theAnalysisManager->Create2DHisto(label,dir,title,"X[km]","Y[km]",
	                                        			nx, xmin, ymax, "lin",
									ny, ymin, ymax));
	 	
	 	HistoList.push_back(UpwardFluxVsPosFlatHisto[i].back());
		UpwardFluxVsPosFlatMaxE[i].push_back(EkinMax);
	 	UpwardFluxVsPosFlatMinE[i].push_back(EkinMin);
	
#ifndef USE_ANALYSIS_ROOT
	 	std::stringstream astream1;
	 	G4String alt,depth,max_e,min_e;
	 	astream1<<(*DetectorAltitudes)[i]/km<<'\t'
	         	<<(*DetectorDepths)[i]*cm2/g<<'\t'
		 	<<EkinMax/MeV<<'\t'<<EkinMin/MeV;
	 	astream1>>alt>>depth>>max_e>>min_e;
	 
	 	UpwardFluxVsPosFlatHisto[i].back()->annotation().addItem("Altitude[km]",alt);
	 	UpwardFluxVsPosFlatHisto[i].back()->annotation().addItem("Depth[g/cm2]",depth); 
	 	UpwardFluxVsPosFlatHisto[i].back()->annotation().addItem("EkinMin [MeV]",min_e);
	 	UpwardFluxVsPosFlatHisto[i].back()->annotation().addItem("EkinMax [MeV]",max_e);
#endif   
	 	DetectUpwardFluxVsPosFlat = true;
	}
  }
  if (ndet==0) {
        G4cout<<"Before creating a flux histogram you should select";
	G4cout<<" at least one detector"<<std::endl;
	return;
  }
  SetDetectionFlux();	

}
//////////////////////////////////////////////////////////////////////////			 	 
//
G4double PLANETOCOSFluxDetectionAnalyser::FindGeoFactor(HISTOBASE* anHisto)
{
  for (unsigned int i=0;i<HistoList.size();i++){
  	if (HistoList[i] == anHisto){
		return GeoFactors[i];
	}
  }
  return 1.;
}
//////////////////////////////////////////////////////////////////////////			 	 
//
void PLANETOCOSFluxDetectionAnalyser::CreateLowestAltitudeNeededHisto(
			      G4String aParticleName,
                              G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType,
			      G4int nAlt,G4double AltitudeMin, G4double AltitudeMax)
{ //Check if particle exists
  //------------------------ 
  G4int PDGCode =theAnalysisManager->GetParticlePDGCode(aParticleName);
  if (PDGCode == -999999) return; 
 
  //create histogram
  //---------------
  DetectLowestAltitudeNeeded=true;	
  LowestAltitudeNeededPDGCode.push_back(PDGCode);
  G4String title= "Lowest altitude neede for detection  of "+aParticleName;
  LowestAltitudeNeededHisto.push_back(theAnalysisManager->Create2DHisto(
						       label,"/LOWESTALTITUDE/"+aParticleName,title,"Ekin[MeV]","lowest altitude needed", 
	                                               nE, E1/MeV, E2/MeV,ScaleType,
						       nAlt, AltitudeMin/km, AltitudeMax/km ));
}

////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSFluxDetectionAnalyser::SetDetectionFlux()
{ if (!DetectFlux){
	DetectFlux = true;
	PLANETOCOSSteppingAction* theSteppingAction 
	       = dynamic_cast<PLANETOCOSSteppingAction*>
	            (G4EventManager::GetEventManager()->GetUserSteppingAction());
	theSteppingAction->SetDetectFlux(true); 
  }  
}
