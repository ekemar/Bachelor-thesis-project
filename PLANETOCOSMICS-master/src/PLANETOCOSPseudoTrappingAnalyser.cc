#include "PLANETOCOSPseudoTrappingAnalyser.hh"
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
PLANETOCOSPseudoTrappingAnalyser::PLANETOCOSPseudoTrappingAnalyser(PLANETOCOSAnalysisManager* anAnalysisManager):
theAnalysisManager(anAnalysisManager)
{ //Post track hit 
  //--------------
  
  DetectPseudoTrapping = false;
  DetectNbOfPlanetTurnVsStartEkin = false;
  DetectNbOfPlanetTurnVsLifeTime = false;
  DetectNbOfEquatorCrossingVsStartEkin = false;
  DetectNbOfEquatorCrossingVsLifeTime = false;
  DetectNbOfEquatorCrossingVsNbOfPlanetTurn = false;
  DetectLifeTimeVsStartEkin = false; 
  		  
  
} 
////////////////////////////////////////////////////////////////////////////////
//  
PLANETOCOSPseudoTrappingAnalyser::~PLANETOCOSPseudoTrappingAnalyser()
{;
}
////////////////////////////////////////////////////////////////////////////////
//  
void PLANETOCOSPseudoTrappingAnalyser::Analyse(PLANETOCOSPostTrackHitsCollection* PostTrackHC)
{ int n_hit = PostTrackHC->entries();
		
  for (int i=0;i<n_hit;i++){
	G4double weight = (*PostTrackHC)[i]->GetWeight();
	G4double start_energy = (*PostTrackHC)[i]->GetStartEnergy()/MeV;
	G4double life_time = (*PostTrackHC)[i]->GetLifeTime()/s;
	/*G4cout<<life_time<<std::endl;
	G4cout<<start_energy<<std::endl;*/
	G4double nb_of_planet_turn = (*PostTrackHC)[i]->GetNbOfPlanetTurn();
	G4double nb_of_equator_crossing = (*PostTrackHC)[i]->GetNbOfEquatorCrossing();	
	G4int PDGCode = (*PostTrackHC)[i]->GetPDGCode();
	if (DetectNbOfPlanetTurnVsStartEkin){
		for (unsigned  int j=0; j <NbOfPlanetTurnVsStartEkinPDGCode.size();j++){
	   			G4int code =NbOfPlanetTurnVsStartEkinPDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D(start_energy, nb_of_planet_turn,weight,
		                        		NbOfPlanetTurnVsStartEkinHisto[j]);
	    	}		
			
			
	}
	if (DetectNbOfPlanetTurnVsLifeTime){
		for (unsigned  int j=0; j <NbOfPlanetTurnVsLifeTimePDGCode.size();j++){
	   			G4int code =NbOfPlanetTurnVsLifeTimePDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D( life_time, nb_of_planet_turn,weight,
		                        		NbOfPlanetTurnVsLifeTimeHisto[j]);
	    	}		
			
			
	}
	if (DetectNbOfEquatorCrossingVsStartEkin){
		for (unsigned  int j=0; j <NbOfEquatorCrossingVsStartEkinPDGCode.size();j++){
	   			G4int code =NbOfEquatorCrossingVsStartEkinPDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D(start_energy, nb_of_equator_crossing,weight,
		                        		NbOfEquatorCrossingVsStartEkinHisto[j]);
	    	}		
			
			
	}
	if (DetectNbOfEquatorCrossingVsLifeTime){
		for (unsigned  int j=0; j <NbOfEquatorCrossingVsLifeTimePDGCode.size();j++){
	   			G4int code =NbOfEquatorCrossingVsLifeTimePDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D( life_time, nb_of_equator_crossing,weight,
		                        		NbOfEquatorCrossingVsLifeTimeHisto[j]);
	    	}		
			
			
	}
	if (DetectNbOfEquatorCrossingVsNbOfPlanetTurn){
		for (unsigned  int j=0; j <NbOfEquatorCrossingVsNbOfPlanetTurnPDGCode.size();j++){
	   			G4int code =NbOfEquatorCrossingVsNbOfPlanetTurnPDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D( nb_of_planet_turn, nb_of_equator_crossing,weight,
		                        		NbOfEquatorCrossingVsNbOfPlanetTurnHisto[j]);
	    	}		
			
			
	}
	if (DetectLifeTimeVsStartEkin){
		for (unsigned  int j=0; j <LifeTimeVsStartEkinPDGCode.size();j++){
	   			G4int code =LifeTimeVsStartEkinPDGCode[j];
	     			if (code == PDGCode || code == -1)
	                 		theAnalysisManager->FillHistogram2D( start_energy, life_time,weight,
		                        		LifeTimeVsStartEkinHisto[j]);
	    	}		
			
			
	}
				 
		
  } 
}
////////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSPseudoTrappingAnalyser::CreateNbOfPlanetTurnVsStartEkinHisto(G4String aParticleName,G4String label,
			                     G4int nE,G4double E1,G4double E2, 
					     G4String ScaleType,
			                     G4int nBinTurn,G4double nTurnMax )
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
   
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName;
 G4String title="Nb of planet turn before stop vs start Ekin for "+aParticleName;
 NbOfPlanetTurnVsStartEkinPDGCode.push_back(PDGCode);
 NbOfPlanetTurnVsStartEkinHisto.push_back(
	       theAnalysisManager->Create2DHisto(label,dir,title,"Ekin [MeV]","Nb of Planet Turn ",
                                                 nE, E1/MeV, E2/MeV, ScaleType,
			                         nBinTurn, 0., nTurnMax)
					         );
 DetectPseudoTrapping = true;
 DetectNbOfPlanetTurnVsStartEkin = true;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////			    
void PLANETOCOSPseudoTrappingAnalyser::CreateNbOfPlanetTurnVsLifeTimeHisto(G4String aParticleName,G4String label,
			                    G4int nTime,G4double time1,G4double time2, 
					    G4String ScaleType,
			                    G4int nBinTurn,G4double nTurnMax )
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName; 
 NbOfPlanetTurnVsLifeTimePDGCode.push_back(PDGCode);
 G4String title= "Nb of planet turn before stop vs lifetime for ";
 title += aParticleName;
 NbOfPlanetTurnVsLifeTimeHisto.push_back(
	theAnalysisManager->Create2DHisto(label,dir,title,"Time[s]","Nb of Planet Turn ",
                                          nTime, time1/s, time2/s, ScaleType,
			                  nBinTurn, 0., nTurnMax)
					  );
 DetectPseudoTrapping = true;
 DetectNbOfPlanetTurnVsLifeTime = true;
 
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPseudoTrappingAnalyser::CreateNbOfEquatorCrossingVsStartEkinHisto(G4String aParticleName,G4String label,
			                     G4int nE,G4double E1,G4double E2, 
					     G4String ScaleType,
			                     G4int nCrossingMax )
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName; 
 NbOfEquatorCrossingVsStartEkinPDGCode.push_back(PDGCode);
 G4String title= "Nb of equator crossing before stop vs start ekin for ";
 title += aParticleName;
	NbOfEquatorCrossingVsStartEkinHisto.push_back(
		theAnalysisManager->Create2DHisto(label,dir,title,"Ekin [MeV]","Nb of equator crossing",
                                                 nE, E1/MeV, E2/MeV, ScaleType,
			                         nCrossingMax+1, -0.5,
						 nCrossingMax+0.5)
					         );
 DetectPseudoTrapping = true;
 DetectNbOfEquatorCrossingVsStartEkin = true;

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////			    
void PLANETOCOSPseudoTrappingAnalyser::CreateNbOfEquatorCrossingVsLifeTimeHisto(G4String aParticleName,G4String label,
			                    G4int nTime,G4double time1,G4double time2, 
					    G4String ScaleType,
			                    G4int nCrossingMax )
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName; 
 NbOfEquatorCrossingVsLifeTimePDGCode.push_back(PDGCode);
 G4String title= "Nb of equator crossing  before stop vs lifetime for ";
 title += aParticleName;
 NbOfEquatorCrossingVsLifeTimeHisto.push_back(
	       theAnalysisManager->Create2DHisto(label,dir,title,"Time[s]","Nb of equator crossing ",
                                                 nTime, time1/s, time2/s, ScaleType,
			                         nCrossingMax+1, -0.5,
						 nCrossingMax+0.5));
 DetectPseudoTrapping = true;
 DetectNbOfEquatorCrossingVsLifeTime = true;
 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////			    
void PLANETOCOSPseudoTrappingAnalyser::CreateNbOfEquatorCrossingVsNbOfPlanetTurnHisto(G4String aParticleName,G4String label,
			                    G4int nBinTurn,G4double nTurnMax,
			                    G4int nCrossingMax )
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName; 
 NbOfEquatorCrossingVsNbOfPlanetTurnPDGCode.push_back(PDGCode);
 G4String title= "Nb of equator crossing vs nb of planet turn before stop";
 title += aParticleName;
 NbOfEquatorCrossingVsNbOfPlanetTurnHisto.push_back(
	      theAnalysisManager->Create2DHisto(label,dir,title,"Nb of Planet Turn","Nb of equator crossing ",
                                                 nBinTurn, 0., nTurnMax,"lin",
			                         nCrossingMax+1, -0.5,
						 nCrossingMax+0.5)
					         );
 DetectPseudoTrapping = true;
 DetectNbOfEquatorCrossingVsNbOfPlanetTurn = true;
 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////			    
void PLANETOCOSPseudoTrappingAnalyser::CreateLifeTimeVsStartEkinHisto(G4String aParticleName,G4String label,
			                G4int nE,G4double E1,G4double E2, 
					G4String ScaleType1,
			                G4int nTime, G4double time1,G4double time2, 
					G4String ScaleType2)
{//Check if particle exists
 //------------------------ 
 G4int PDGCode = theAnalysisManager->GetParticlePDGCode(aParticleName);
 if (PDGCode == -999999  ) return;
 
 //check the type of geometry
 //--------------------------
 if (!theAnalysisManager->CheckIfSphericalGeometry()) return;
 
 //create histogram
 //----------------
 G4String dir ="/PSEUDOTRAPPING/"+aParticleName; 
 LifeTimeVsStartEkinPDGCode.push_back(PDGCode);
 G4String title= "Lifetime vs start Ekin for ";
 title += aParticleName;
 G4cout<<time1/s<<std::endl;
 G4cout<<time2/s<<std::endl;
 LifeTimeVsStartEkinHisto.push_back(
	       theAnalysisManager->Create2DHisto(label,dir,title,"Ekin[MeV]","LifeTime[s]",
                                                 nE, E1, E2, ScaleType1,
			                         nTime, time1/s, time2/s, ScaleType2));
 DetectPseudoTrapping = true;
 DetectLifeTimeVsStartEkin = true;
}
