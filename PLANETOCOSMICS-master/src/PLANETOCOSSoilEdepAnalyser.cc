#include "PLANETOCOSSoilEdepAnalyser.hh"
#include "myfunctions.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinatePlanet.hh"
#include "myfunctions.hh"
#include "G4ParticleTable.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "PLANETOCOSEventAction.hh"
#include "PLANETOCOSFluxHit.hh"
#include "PLANETOCOSSoilSD.hh"
#include "PlanetManager.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PlanetUnits.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSoilEdepAnalyser::PLANETOCOSSoilEdepAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager):
theAnalysisManager(anAnalysisManager)
{ 
  		  
  //Edep detection
  //--------------  
  DetectEdep = false;
  DetectEdepVsThickness = false;
  DetectEdepVsDepth = false;
  EdepVsThicknessHisto =0;
  EdepVsDepthHisto =0; 
  
  
  //Bin content for histogram 
  // 1 is for unnormalised case
  // 2 is for normalisation to primary flux
  // 3 is for normalisation to one  pimary particle
  //-------------------
  EdepBin1 ="Edep[rad*cm2]";
  EdepBin2 ="Edep[rad/s]";
  EdepBin3 ="Edep[rad*cm2/nb_prim_part]";
} 
////////////////////////////////////////////////////////////////////////////////
//  
PLANETOCOSSoilEdepAnalyser::~PLANETOCOSSoilEdepAnalyser()
{;
}
////////////////////////////////////////////////////////////////////////////////
//  
void PLANETOCOSSoilEdepAnalyser::Analyse(PLANETOCOSEdepHitsCollection* EdepHC)
{  int n_hit = EdepHC->entries();
   //G4cout<<"test1"<<std::endl;
   for (int i=0;i<n_hit;i++){
      	G4double weight = (*EdepHC)[i]->GetWeight();
      
      	//the normalise histogram will be given in rad/s
       	G4double edep = weight*(*EdepHC)[i]->GetEdep()/erg/100.;
        
       	if (DetectEdepVsThickness){
        	G4double h1= (*EdepHC)[i]->GetAltitude1()/km;
	 	G4double h2= (*EdepHC)[i]->GetAltitude2()/km;
		double xmin =std::min(h1,h2);
	 	double xmax = std::max(h1,h2);
		
#ifndef USE_ANALYSIS_ROOT		
		int nbins=EdepVsThicknessHisto->axis().bins();
		int n1=EdepVsThicknessHisto->coordToIndex(xmin);
         	int n2=EdepVsThicknessHisto->coordToIndex(xmax);		
#else	
		int nbins=EdepVsThicknessHisto->GetXaxis()->GetNbins();
		int n1=EdepVsThicknessHisto->FindBin(xmin)-1;
    		int n2=EdepVsThicknessHisto->FindBin(xmax)-1;
#define fill Fill		
#endif	
		
	 	double up =SoilLayerDepthsInLengthUnit->back()/km;
	 	double low =SoilLayerDepthsInLengthUnit->front()/km;
         	double dx_bin=(up-low)/nbins;
	 	
		//G4cout<<"fill value"<<edep/DepthForThicknessHisto[n1]<<std::endl;
	 	if (n1==n2) 
			EdepVsThicknessHisto->fill(xmin,edep/DepthForThicknessHisto[n1]); 
         	else { 
	   		//G4cout<<"dx_bin "<<dx_bin<<std::endl;
	    		double weight_per_bin=edep*dx_bin/(xmax-xmin);
            		double x1=low+(n1+1)*dx_bin;
            		double x2=low+(n2)*dx_bin;
            		double weight_n1=edep*(x1-xmin)/(xmax-xmin);
            		double weight_n2=edep*(xmax-x2)/(xmax-xmin);
	    		EdepVsThicknessHisto->fill(xmin,weight_n1/DepthForThicknessHisto[n1]);
	    		EdepVsThicknessHisto->fill(xmax,weight_n2/DepthForThicknessHisto[n2]);
	    		for (int j=n1+1;j<n2;j++){
	      			double x=low+ (j+0.5)*dx_bin;
	       			EdepVsThicknessHisto->fill(x,weight_per_bin/DepthForThicknessHisto[j]);
	      		}	
	   	}
		
	
	 }
       
       	if (DetectEdepVsDepth){
		
        	G4double depth1= (*EdepHC)[i]->GetDepth1()*cm2/g;
	 	G4double depth2= (*EdepHC)[i]->GetDepth2()*cm2/g;
	 	myfunc::LinearDistributionInHistogram(EdepVsDepthHisto,
	                                       				depth1,
					       				depth2,
					       				edep/delta_depth);
								
		
	}
   }
   //G4cout<<"test2"<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSSoilEdepAnalyser::CreateEdepVsThickness(G4String label,
                                                    G4int n_Bins)
{ DetectEdep = true;
  
  if (!DetectEdepVsThickness) {
    	DetectEdepVsThickness = true;
    	PLANETOCOSSoilSD* sd  = dynamic_cast<PLANETOCOSSoilSD*>
                  (G4SDManager::GetSDMpointer()->FindSensitiveDetector("soilSD"));
   	if (sd) sd->SetDetectEdepVsThickness(true);
    
 }  
 
 if (EdepVsThicknessHisto) delete EdepVsThicknessHisto;
 	
 //create histogram 
 G4String title ="Deposited energy vs soil thickness";
 G4String dir ="/SOIL_EDEP/DEPTH_IN_KM";
 G4double h_min = SoilLayerDepthsInLengthUnit->front();
 G4double h_max = SoilLayerDepthsInLengthUnit->back();
 DepthForThicknessHisto.clear();
 G4double dh=(h_max-h_min)/n_Bins;
 G4double old_depth=0.;
 for (int i=0; i<n_Bins;i++){
  	G4double h =dh*(i+1);
   	G4double depth = myfunc::LinearInterpolation(*SoilLayerDepthsInLengthUnit,
                                                *SoilLayerDepthsInDepthUnit,h);
   	G4double ddepth=depth-old_depth;
   	old_depth =depth;						
  	DepthForThicknessHisto.push_back(ddepth*cm2/g);
 }
 
 EdepVsThicknessHisto =
              theAnalysisManager->Create1DHisto(label, dir,title,"depth[km]",n_Bins,h_min/km,h_max/km,"linear");
}
////////////////////////////////////////////////////////////////////////////////
//		    
void PLANETOCOSSoilEdepAnalyser::CreateEdepVsDepth(G4String label,G4int n_Bins)
{ 
  DetectEdep = true;
  if (!DetectEdepVsDepth){ 
    	DetectEdepVsDepth = true;
    	PLANETOCOSSoilSD* sd  = dynamic_cast<PLANETOCOSSoilSD*>
                  	(G4SDManager::GetSDMpointer()->FindSensitiveDetector("soilSD"));
    	sd->SetDetectEdepVsDepth(true);
  }  

  //remove histogram if exist  
  if (EdepVsDepthHisto) delete EdepVsDepthHisto;
  
  //create histogram 
 
  G4String title ="Deposited energy vs depth";
  G4String dir ="/SOIL_EDEP/DEPTH_IN_G_PER_CM2";
  G4double depth_min = SoilLayerDepthsInDepthUnit->front();
  G4double depth_max = SoilLayerDepthsInDepthUnit->back();
  delta_depth = (depth_max-depth_min)*cm2/g/n_Bins;
  G4double depth1 = double( int(depth_min*cm2/g));
  G4double depth2 = double (int(depth_max*cm2/g) + 0.0001);
  EdepVsDepthHisto= theAnalysisManager->Create1DHisto(label,dir,title,"depth[g/cm2]",
	               				      n_Bins,
		       				      depth1,
		       				      depth2,
						      "linear");

  
   
   
  
  						      
}
