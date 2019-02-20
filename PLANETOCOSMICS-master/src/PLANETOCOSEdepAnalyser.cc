#include "PLANETOCOSEdepAnalyser.hh"
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
#include "PLANETOCOSSD.hh"
#include "PlanetManager.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PlanetUnits.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEdepAnalyser::PLANETOCOSEdepAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager):
theAnalysisManager(anAnalysisManager)
{ 
  		  
  //Edep detection
  //--------------  
  DetectEdep = false;
  DetectEdepvsAltitude = false;
  DetectEdepvsDepth = false;
  EdepVsAltitudeHisto =0;
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
PLANETOCOSEdepAnalyser::~PLANETOCOSEdepAnalyser()
{;
}
////////////////////////////////////////////////////////////////////////////////
//  
void PLANETOCOSEdepAnalyser::Analyse(PLANETOCOSEdepHitsCollection* EdepHC)
{  int n_hit = EdepHC->entries();
  
   for (int i=0;i<n_hit;i++){
      	G4double weight = (*EdepHC)[i]->GetWeight();
      
      	//the normalise histogram will be given in rad/s
       	G4double edep = weight*(*EdepHC)[i]->GetEdep()/erg/100.;
        
       	if (DetectEdepvsAltitude){
        	G4double alt1= (*EdepHC)[i]->GetAltitude1()/km;
	 	G4double alt2= (*EdepHC)[i]->GetAltitude2()/km;
		double xmin =std::min(alt1,alt2);
	 	double xmax = std::max(alt1,alt2);
		
#ifndef USE_ANALYSIS_ROOT		
		int nbins=EdepVsAltitudeHisto->axis().bins();
		int n1=EdepVsAltitudeHisto->coordToIndex(xmin);
         	int n2=EdepVsAltitudeHisto->coordToIndex(xmax);		
#else	
		int nbins=EdepVsAltitudeHisto->GetXaxis()->GetNbins();
		int n1=EdepVsAltitudeHisto->FindBin(xmin)-1;
    		int n2=EdepVsAltitudeHisto->FindBin(xmax)-1;
#define fill Fill		
#endif	
		
	 	double low =AtmosphericLayerAltitudes->back()/km;
	 	double up =AtmosphericLayerAltitudes->front()/km;
         	double dx_bin=(up-low)/nbins;
	 	
		//G4cout<<"fill value"<<edep/DepthForAltitudeHisto[n1]<<std::endl;
	 	if (n1==n2) 
			EdepVsAltitudeHisto->fill(xmin,edep/DepthForAltitudeHisto[n1]); 
         	else { 
	   		//G4cout<<"dx_bin "<<dx_bin<<std::endl;
	    		double weight_per_bin=edep*dx_bin/(xmax-xmin);
            		double x1=low+(n1+1)*dx_bin;
            		double x2=low+(n2)*dx_bin;
            		double weight_n1=edep*(x1-xmin)/(xmax-xmin);
            		double weight_n2=edep*(xmax-x2)/(xmax-xmin);
	    		EdepVsAltitudeHisto->fill(xmin,weight_n1/DepthForAltitudeHisto[n1]);
	    		EdepVsAltitudeHisto->fill(xmax,weight_n2/DepthForAltitudeHisto[n2]);
	    		for (int j=n1+1;j<n2;j++){
	      			double x=low+ (j+0.5)*dx_bin;
	       			EdepVsAltitudeHisto->fill(x,weight_per_bin/DepthForAltitudeHisto[j]);
	      		}	
	   	}
		
	
	 }
       
       	if (DetectEdepvsDepth){
		
		
        	G4double depth1= (*EdepHC)[i]->GetDepth1()*cm2/g;
	 	G4double depth2= (*EdepHC)[i]->GetDepth2()*cm2/g;
	 	myfunc::LinearDistributionInHistogram(EdepVsDepthHisto,
	                                       				depth1,
					       				depth2,
					       				edep/delta_depth);
		
#ifdef  TEST_EDEP
	int ParticleTypeCode  = (*EdepHC)[i]->GetParticleTypeCode();
	if ((ParticleTypeCode) >-1) { 
		myfunc::LinearDistributionInHistogram(TestEdepVsDepthHisto[ParticleTypeCode],
	                                       				depth1,
					       				depth2,
					       				edep/delta_depth);
		
	}								
#endif									
		
	}
   }
}
////////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSEdepAnalyser::CreateEdepVsAltitude(G4String label,
                                                    G4int n_Alt)
{ DetectEdep = true;
  
  if (!DetectEdepvsAltitude) {
    	DetectEdepvsAltitude = true;
    	PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                  (G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
   	if (sd) sd->SetDetectEdepvsAltitude(true);
    
 }  
 
 if (EdepVsAltitudeHisto) delete EdepVsAltitudeHisto;
 	
 //create histogram 
 G4String title ="Deposited energy vs altitude";
 G4String dir ="/EDEP/ALTITUDE";
 G4double Alt_max = AtmosphericLayerAltitudes->front();
 G4double Alt_min = AtmosphericLayerAltitudes->back();
 DepthForAltitudeHisto.clear();
 G4double dalt=(Alt_max-Alt_min)/n_Alt;
 G4double depth_alt0= AtmosphericLayerDepths->back();
 
 for (int i=0; i<n_Alt;i++){
  	G4double altitude =Alt_min+(i+1)*dalt;
   	G4double depth = myfunc::LinearInterpolation(*AtmosphericLayerAltitudes,
                                                *AtmosphericLayerDepths,altitude);
   	G4double ddepth= depth_alt0-depth;
   	depth_alt0 =depth;						
  	DepthForAltitudeHisto.push_back(ddepth*cm2/g);
 }
 
 EdepVsAltitudeHisto =
              theAnalysisManager->Create1DHisto(label, dir,title,"altitude[km]",n_Alt,Alt_min/km,Alt_max/km,"linear");
}
////////////////////////////////////////////////////////////////////////////////
//		    
void PLANETOCOSEdepAnalyser::CreateEdepVsDepth(G4String label,G4int n_depth)
{ 
  DetectEdep = true;
  if (!DetectEdepvsDepth){ 
    	DetectEdepvsDepth = true;
    	PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                  	(G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
    	sd->SetDetectEdepvsDepth(true);
  }  

  //remove histogram if exist  
  if (EdepVsDepthHisto) delete EdepVsDepthHisto;
  
  //create histogram 
 
  G4String title ="Deposited energy vs depth";
  G4String dir ="/EDEP/DEPTH";
  G4double depth_min = AtmosphericLayerDepths->front();
  G4double depth_max = AtmosphericLayerDepths->back();
  delta_depth = (depth_max-depth_min)*cm2/g/n_depth;
  G4double depth1 = double( int(depth_min*cm2/g));
  G4double depth2 = double (int(depth_max*cm2/g) + 1);
  EdepVsDepthHisto= theAnalysisManager->Create1DHisto(label,dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear");
#ifdef  TEST_EDEP    
   
   TestEdepVsDepthHisto.clear();
   
   dir ="/EDEP/DEPTH/gamma/";
   title="Deposited energy by non cut gamma  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut gamma  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   						      
   
   dir ="/EDEP/DEPTH/e-/";
   title="Deposited energy by non cut e-  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut e-  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   dir ="/EDEP/DEPTH/e+/";
   title="Deposited energy by non cut e+  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut e+  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   dir ="/EDEP/DEPTH/mu/";
   title="Deposited energy by non cut mu  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut mu  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   dir ="/EDEP/DEPTH/tau/";
   title="Deposited energy by non cut tau  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut tau  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   
   dir ="/EDEP/DEPTH/proton/";
   title="Deposited energy by non cut proton  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut proton  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   dir ="/EDEP/DEPTH/neutron/";
   title="Deposited energy by non cut neutron  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut neutron  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
   dir ="/EDEP/DEPTH/other_baryon/";
   title="Deposited energy by non cut other_baryon  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut other_baryon  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
						      
  
  
   dir ="/EDEP/DEPTH/ions/";
   title="Deposited energy by non cut ions  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut ions  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
    
						      
  
  dir ="/EDEP/DEPTH/mesons/";
   title="Deposited energy by non cut mesons  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("1",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));
   title="Deposited energy by cut mesons  vs depth";
   TestEdepVsDepthHisto.push_back(theAnalysisManager->Create1DHisto("2",dir,title,"depth[g/cm2]",
	               				      n_depth,
		       				      depth1,
		       				      depth2,
						      "linear"));					      
   						      
						      
  						      
						      

#endif  
  
  
   
   
  
  						      
}
