#include "PLANETOCOSPrimaryFluxAnalyser.hh"
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

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryFluxAnalyser::PLANETOCOSPrimaryFluxAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager):
theAnalysisManager(anAnalysisManager)
{//Primary Detection
  //--------------------
  Initialise();
  
  //Bin content for histogram 
  // 1 is for unnormalised case
  // 2 is for normalisation to primary flux
  // 3 is for normalisation to one  pimary particle
  //-------------------
 
  DirDifFluxBin1 ="Flux[#]";
  DirDifFluxBin2 ="Flux[#/cm2/s]";
  DirDifFluxBin3 ="Flux[#/nb_primary_part]";
  
  
  OmniDifFluxBin1 ="Flux[#]";
  OmniDifFluxBin2 ="Flux[#/cm2/s]";
  OmniDifFluxBin3 ="Flux[#/nb_primary_part]";
  
  IntegFluxBin1 ="Flux[#]";
  IntegFluxBin2 ="Flux[#/cm2/s]";
  IntegFluxBin3 ="Flux[#/nb_primary_part]";
       
}
////////////////////////////////////////////////////////////////////////////////
//  
PLANETOCOSPrimaryFluxAnalyser::~PLANETOCOSPrimaryFluxAnalyser() 
{;
}
////////////////////////////////////////////////////////////////////////////////
//  
void PLANETOCOSPrimaryFluxAnalyser::Analyse(PLANETOCOSPrimaryHitsCollection* PrimaryHC) 
{ //G4cout<<"PrimTest1"<<std::endl;
  int n_hit = PrimaryHC->entries();
  //G4cout<<"PrimTest2"<<std::endl;
  for (int i=0;i<n_hit;i++){
  
      	G4double weight = (*PrimaryHC)[i]->GetWeight();
       	G4double energy = (*PrimaryHC)[i]->GetEnergy();
       	G4double zenith = (*PrimaryHC)[i]->GetTheta();
       	G4double cos_zen = std::cos(zenith); 
	//G4cout<<"cos_zen "<<cos_zen<<std::endl;
       	G4double azimuth = (*PrimaryHC)[i]->GetPhi();
       	G4double PDGCode = (*PrimaryHC)[i]->GetPDGCode();
  	
        //G4cout<<"PrimTest3"<<std::endl;
       	if (DetectEnergyVsZenithPrimary){
		for (unsigned  int j=0; j <EnergyVsZenithPrimaryPDGCode.size();j++){
	    		G4int code =EnergyVsZenithPrimaryPDGCode[j];
	     		if (code == PDGCode || code == -1)
	                 	theAnalysisManager->FillHistogram2D(energy/MeV,cos_zen,weight,
			                 	                    EnergyVsZenithPrimaryHisto[j]);
	     
            	}
	 }
	//G4cout<<"PrimTest4"<<std::endl;
        if (DetectZenithVsAzimuthPrimary){
         	for (unsigned  int j=0; j <ZenithVsAzimuthPrimaryPDGCode.size();j++){
	   		G4int code =ZenithVsAzimuthPrimaryPDGCode[j];
	     		if (code == PDGCode || code == -1)
	                 	theAnalysisManager->FillHistogram2D(azimuth/degree,cos_zen,weight,
		                 	      		ZenithVsAzimuthPrimaryHisto[j]);
	    	}
	 }
	//G4cout<<"PrimTest5"<<std::endl;
        if (DetectOmniFluxPrimary){
         	for (unsigned  int j=0; j < OmniFluxPrimaryPDGCode.size();j++){
	    		G4int code =OmniFluxPrimaryPDGCode[j];
	   
	     		if (code == PDGCode || code == -1)
	       			theAnalysisManager->FillHistogram1D(energy/MeV, weight, OmniFluxPrimaryHisto[j]);
	     
	    	}
	}
	//G4cout<<"PrimTest6"<<std::endl;
	if (DetectZenithPrimary){
         	for (unsigned  int j=0; j <ZenithPrimaryPDGCode.size();j++){
	    		//G4cout<<"PrimTest61"<<std::endl;
			G4int code =ZenithPrimaryPDGCode[j];
			//G4cout<<"PrimTest62"<<std::endl;
			//G4cout<<ZenithPrimaryHisto[j]<<std::endl;
	     		if (code == PDGCode || code == -1)
	        		theAnalysisManager->FillHistogram1D(cos_zen,weight, ZenithPrimaryHisto[j]);		
            		//G4cout<<"PrimTest63"<<std::endl;
		}
	} 
	//G4cout<<"PrimTest7"<<std::endl;
	if (DetectAzimuthPrimary){
         	for (unsigned  int j=0; j < AzimuthPrimaryPDGCode.size();j++){
	    		G4int code =AzimuthPrimaryPDGCode[j];
	     		if (code == PDGCode || code == -1)
	        		theAnalysisManager->FillHistogram1D(azimuth/degree,weight,AzimuthPrimaryHisto[j]);
	    	}
	} 
	//G4cout<<"PrimTest8"<<std::endl;
 }
}
////////////////////////////////////////////////////////////////////////////////  
//			 
void PLANETOCOSPrimaryFluxAnalyser::CreateEnergyVsZenithPrimaryHisto(G4String aParticleName,
                                 G4String label,
			         G4int nE,G4double E1,G4double E2, G4String ScaleType, 
			         G4int nZenith, G4double cosZen1, G4double cosZen2)
{ 
  G4int  PDGCode = GetPDGCodeAndInitialise(aParticleName); 
  if (PDGCode == -999999) return;
  
  EnergyVsZenithPrimaryPDGCode.push_back(PDGCode);
  G4String dir_name = "/PRIMARY/"+aParticleName;
  G4String histo_name= "Primary flux of "+aParticleName+" cos_zen vs ekin ";
  EnergyVsZenithPrimaryHisto.push_back(theAnalysisManager->Create2DHisto(label,dir_name,histo_name,"Energy[MeV]","cos_theta",
	                                                                 nE,E1/MeV,E2/MeV,ScaleType, 
 		                                                         nZenith, cosZen1, cosZen2));
  DetectEnergyVsZenithPrimary = true;
 
}  
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryFluxAnalyser::CreateZenithVsAzimuthPrimaryHisto(G4String aParticleName,
                                  G4String label, 
			          G4int nZenith, G4double cosZen1, G4double cosZen2,
				  G4int nAzim, G4double Azim1, G4double Azim2)
{ G4int  PDGCode = GetPDGCodeAndInitialise(aParticleName); 
  if (PDGCode == -999999) return;
  G4String dir_name = "/PRIMARY/"+aParticleName;
  ZenithVsAzimuthPrimaryPDGCode.push_back(PDGCode);
  G4String histo_name= "Primary flux of "+aParticleName+" cos_zen vs azim";
  ZenithVsAzimuthPrimaryHisto.push_back(theAnalysisManager->Create2DHisto(label,dir_name,histo_name,"phi [deg]","cos_theta",  
	                                              nAzim,Azim1,Azim2,"lin", 
			                              nZenith, cosZen1, cosZen2));
   
  DetectZenithVsAzimuthPrimary = true;
}
////////////////////////////////////////////////////////////////////////////////
//	
void PLANETOCOSPrimaryFluxAnalyser::OmnifluxFluxPrimaryHisto(G4String aParticleName,
                              G4String label,
			      G4int nE,G4double E1,G4double E2, G4String ScaleType) 			 
{ G4int  PDGCode = GetPDGCodeAndInitialise(aParticleName); 
  if (PDGCode == -999999) return;
  G4String dir_name = "/PRIMARY/"+aParticleName;
  OmniFluxPrimaryPDGCode.push_back(PDGCode);
  G4String histo_name= "Primary flux of "+aParticleName;
  OmniFluxPrimaryHisto.push_back(theAnalysisManager->Create1DHisto(label,dir_name,histo_name,"Energy[MeV]",
	                                      nE,E1/MeV,E2/MeV,ScaleType));
  DetectOmniFluxPrimary = true;
}
////////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSPrimaryFluxAnalyser::CreateZenithPrimaryHisto(G4String aParticleName,
                         G4String label,
			 G4int nZenith, G4double cosZen1, G4double cosZen2)
{ G4int  PDGCode = GetPDGCodeAndInitialise(aParticleName); 
  if (PDGCode == -999999) return;
  G4String dir_name = "/PRIMARY/"+aParticleName;  
  ZenithPrimaryPDGCode.push_back(PDGCode);
  G4String histo_name= "Primary zenith angular distribution of "+aParticleName;
  ZenithPrimaryHisto.push_back(theAnalysisManager->Create1DHisto(label,dir_name,histo_name,"cos_theta",
	                                  nZenith,cosZen1,cosZen2,"lin"));
  DetectZenithPrimary = true;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryFluxAnalyser::CreateAzimuthPrimaryHisto(G4String aParticleName,
                         G4String label, 
			 G4int nAzim, G4double Azim1, G4double Azim2)
{ G4int  PDGCode = GetPDGCodeAndInitialise(aParticleName); 
  if (PDGCode == -999999) return;
  G4String dir_name = "/PRIMARY/"+aParticleName;  
  AzimuthPrimaryPDGCode.push_back(PDGCode);
  G4String histo_name= "Primary azimuth angular distribution of "+aParticleName;
  AzimuthPrimaryHisto.push_back(theAnalysisManager->Create1DHisto(label,dir_name,histo_name,"phi[deg]",
	                                      			  nAzim,Azim1,Azim2,"lin"));
  DetectAzimuthPrimary = true;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int PLANETOCOSPrimaryFluxAnalyser::GetPDGCodeAndInitialise(G4String aParticleName){
  
  G4int PDGCode;
  PDGCode =-1;
  if (aParticleName != "all") {
   	G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    	G4ParticleDefinition* aPartDefinition = 
                              theParticleTable->FindParticle(aParticleName);
    	if 	(aPartDefinition) PDGCode = aPartDefinition->GetPDGEncoding();	
    	else {
     		G4cout<<"The particle "<<aParticleName<<" is not defined"<<std::endl;
      		G4cout<<"No primary  histogram will be  created"<<std::endl;
      		PDGCode=-999999;
		return PDGCode;
     	}	      
  }

  //Initialise if needed
  //---------------
  if (!DetectPrimary) {
  	DetectPrimary=true;
    	PLANETOCOSPrimaryGeneratorAction* thePrimaryAction = 
				(PLANETOCOSPrimaryGeneratorAction*)
	       					(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    	thePrimaryAction->SetDetectPrimary(true);
  
  }
  return PDGCode;
}   
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryFluxAnalyser::Initialise()
{ EnergyVsZenithPrimaryHisto.clear();
  EnergyVsZenithPrimaryPDGCode.clear();
  DetectEnergyVsZenithPrimary = false;
 
  ZenithVsAzimuthPrimaryHisto.clear();
  ZenithVsAzimuthPrimaryPDGCode.clear();
  DetectZenithVsAzimuthPrimary = false;
  
  OmniFluxPrimaryHisto.clear();
  OmniFluxPrimaryPDGCode.clear();
  DetectOmniFluxPrimary = false;
  
  ZenithPrimaryHisto.clear();
  ZenithPrimaryPDGCode.clear();
  DetectZenithPrimary = false;
  
  AzimuthPrimaryHisto.clear();
  AzimuthPrimaryPDGCode.clear();
  DetectAzimuthPrimary = false;
}
