#include "MagneticShieldingToolAnalysisManager.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "PLANETOCOSTrackInformation.hh"
#include "SpaceCoordinatePlanet.hh"
#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "PlanetUnits.hh"
#include "PlanetManager.hh"


MagneticShieldingToolAnalysisManager* MagneticShieldingToolAnalysisManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
MagneticShieldingToolAnalysisManager::MagneticShieldingToolAnalysisManager()  
{
}
////////////////////////////////////////////////////////////////////////////////
//
MagneticShieldingToolAnalysisManager::~MagneticShieldingToolAnalysisManager() 
{
}
////////////////////////////////////////////////////////////////////////////////
//
MagneticShieldingToolAnalysisManager* MagneticShieldingToolAnalysisManager::GetInstance()
{
  if (instance == 0) instance = new MagneticShieldingToolAnalysisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::OpenAsymptoticDirectionFile(G4String fileName)
{ // Asymptotic direction file
  theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"Rigidity"<<'\t'<<"Filter"<<'\t'<<"Asympt. Lat."<<'\t'
                    <<"Asympt. Long."<<'\t'
                    <<"Position Xpla"<<'\t'<<"Ypla"<<'\t'<<"Zpla"
		    <<std::endl;
		    
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue)
{   //Rigidity 
    
    G4double rigidity =LastMomentaOnCharge.mag()/1000.;
    
    //Asymptotic direction
    
    G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
    G4double Aslat=90.-(AsympDir.theta()/degree);
    G4double Aslong=AsympDir.phi()/degree;
    
    
    //LastPosition
    
    G4ThreeVector LastPLAPosition;
    LastPLAPosition = LastPosition/PlanetManager::GetInstance()->GetRplanet();
    
    
    //Nb of turn and equatorial crossing
    //-----------------------------------
    
    /*const G4Track* aTrack = G4EventManager::GetEventManager()
    					   ->GetTrackingManager()
					   ->GetTrack();
        				
    PLANETOCOSTrackInformation* aTrackInfo 
  			= dynamic_cast<PLANETOCOSTrackInformation*>
				            (aTrack->GetUserInformation());
    double NbOfPlanetTurn = aTrackInfo->GetNbOfPlanetTurn();
    double NbOfEquatorCrossing = aTrackInfo->GetNbOfEquatorCrossing();*/
    //write on asymptotic direction file 
    
	
    //theAsciiFile<<std::setiosflags(std::ios::scientific);;
    theAsciiFile<<std::setprecision(4);				   
    theAsciiFile<<rigidity<<'\t'<<'\t';
    theAsciiFile<<FilterValue<<'\t';
    theAsciiFile<<Aslat<<'\t'<<'\t'<<Aslong<<'\t'<<'\t'<<'\t';
    theAsciiFile<<LastPLAPosition.x()<<'\t'
                   <<LastPLAPosition.y()<<'\t'
		   <<LastPLAPosition.z()<<'\t'
		//   <<NbOfPlanetTurn<<'\t'
		//   <<NbOfEquatorCrossing<<'\t'
		   <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::OpenAsymptoticDirVsDirFile(G4String coord_sys,G4double rigidity,G4String fileName)
{ // Asymptotic direction file
  theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"The coordinate system for defining arrival position is "<< coord_sys<<G4endl;
  theAsciiFile<<"rigidity: "<<rigidity /GV <<" GV " <<G4endl;
  theAsciiFile<<"zenith"<<'\t'<<"azimuth"<<'\t'<<"filter"<<'\t'<<"Asympt. Lat. PLA"<<'\t'
             <<"Asympt. Long. PLA"<<'\t'
             <<"Position Xgeo"<<'\t'<<"Ygeo"<<'\t'<<"Zgeo"<<G4endl;
		    
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::RegisterAsymptoticDirVsDir(G4double zenith, G4double azimuth,
                                                             G4ThreeVector LastPosition,
                                                             G4ThreeVector LastMomentaOnCharge,
			                                     G4int  FilterValue)
{  
    
    //Asymptotic direction
    
    G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
    G4double Aslat=90.-(AsympDir.theta()/degree);
    G4double Aslong=AsympDir.phi()/degree;
    
    
    //LastPosition
    
    G4ThreeVector LastPLAPosition;
    LastPLAPosition = LastPosition/PlanetManager::GetInstance()->GetRplanet();
    
    //write on asymptotic direction file 
    
    theAsciiFile.precision(2);
    theAsciiFile.setf(std::ios::fixed);					   
    theAsciiFile<<zenith/degree<<'\t'<<azimuth/degree<<'\t'<<'\t';
    theAsciiFile<<FilterValue<<'\t';
    theAsciiFile<<Aslat<<'\t'<<'\t'<<Aslong<<'\t'<<'\t'<<'\t';
    theAsciiFile<<LastPLAPosition.x()<<'\t'
                   <<LastPLAPosition.y()<<'\t'
		   <<LastPLAPosition.z()<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::CloseAsymptoticDirectionFile
                                    (G4double Rc,G4double Rm,G4double Rs)				
{//write rigidity cutoff
 
  theAsciiFile<<"Rl "<<Rs<<'\t';
  theAsciiFile<<"Rc "<<Rc<<'\t';
  theAsciiFile<<"Ru "<<Rm<<G4endl;
  
 //close file
  theAsciiFile.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::
            OpenCutoffVsPositionFile(G4String fileName,G4String CoordSys,
	                             G4double Altitude, G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  //theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Altitude"<<'\t'<<Altitude/km<<G4endl;
  theAsciiFile<<"Zenith"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<'\t'<<"Ilat2"<<'\t'<<'\t'<<"Ilat3"<<'\t'<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::OpenCutoffVsSpenvisPositionGridFile(G4String fileName,
							    G4String CoordSys,
	                             			    G4double zenith,
				     			    G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  //theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Zenith [deg]:"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth [deg]"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Altitude [km]"<<'\t'<<"Latitude [deg]"<<'\t'<<"Longitude[deg]"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<"Ilat2"<<'\t'<<"Ilat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::OpenCutoffVsSpenvisTrajectoryFile(G4String fileName,G4String CoordSys,
	                             G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  //theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Zenith [deg]:"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth [deg]:"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Mod Julian Day"<<'\t'<<"Altitude [km]"<<'\t'<<"Latitude [deg]"<<'\t'<<"Longitude[deg]"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"ILat1"<<'\t'<<"ILat2"<<'\t'<<"ILat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::RegisterCutoffVsPosition(G4double Rc,G4double Rm,G4double Rs,
                                  G4double latitude, G4double longitude, std::vector<G4double> ILat)
{ 
 
 theAsciiFile<<latitude/degree<<'\t'<<longitude/degree<<'\t';
 if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<ILat[0]/degree<<'\t'<<ILat[1]/degree<<'\t'<<ILat[2]/degree<<'\t'<<'\t';
 theAsciiFile << Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
				       
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::
            OpenCutoffVsPositionOnLShellFile(G4String fileName,G4String CoordSys,
	                             G4double L, G4double zenith,
				     G4double azimuth)
{ 
  theAsciiFile.open(fileName, std::ios::out);
  //theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4); 
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"L:"<<'\t'<<L<<G4endl;
  theAsciiFile<<"Zenith:"<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth:"<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<"Ilat2"<<'\t'<<"Ilat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::
            OpenCutoffVsDirectionFile(G4String fileName,G4String CoordSys, bool Append)
{
  if(!Append)//PVD
	{//PVD
	theAsciiFile.open(fileName, std::ios::out);
	////theAsciiFile<<std::setiosflags(std::ios::scientific);;
	theAsciiFile<<std::setprecision(4);
	theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
	
	PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
	(PLANETOCOSPrimaryGeneratorAction*)
		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	
	if (CoordSys=="PLAG"){
		G4double PLAGalt;
		G4double PLAGlat;
		G4double PLAGlong;
		thePrimaryGenerator->GetPLAGPosition(PLAGalt,PLAGlat,PLAGlong);
		theAsciiFile<<"Altitude:"<<'\t'<<PLAGalt/km<<" km"<<G4endl;
		theAsciiFile<<"Latitude:"<<'\t'<<PLAGlat/degree<<" deg"<<G4endl;
		theAsciiFile<<"Longitude:"<<'\t'<<PLAGlong/degree<<" deg"<<G4endl;  
	}
	else{ 
		G4ThreeVector PLAposition = thePrimaryGenerator->GetPLAPosition();
		G4ThreeVector position=SpaceCoordinatePlanet::GetInstance()
					->Transform(PLAposition,"PLA",CoordSys)/PlanetManager::GetInstance()->GetRplanet();
			
		theAsciiFile<<"Position:"<<'\t'<<position.x()<<'\t'
					<<position.y()<<'\t'
					<<position.z()<<G4endl;
	
	}
	theAsciiFile<<"Zenith"<<'\t'<<"Azimuth"<<'\t'<<"Ru"<<'\t'
			<<"Rc"<<'\t'<<"Rl"<<G4endl; 
	}
  else
	{
	theAsciiFile.open(fileName, std::ios::out | std::ios::app);
	////theAsciiFile<<std::setiosflags(std::ios::scientific);;
	theAsciiFile<<std::setprecision(4);
	}
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::RegisterCutoffVsDirection
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double zenith, G4double azimuth)
{ theAsciiFile<<zenith/degree<<'\t'<<'\t'<<azimuth/degree<<'\t'<<
                                            Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;		    
					    
 					       
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::
            OpenCutoffVsTimeFile(G4String fileName)
{ theAsciiFile.open(fileName, std::ios::out);
  //theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
 	(PLANETOCOSPrimaryGeneratorAction*)
          	G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
 
  G4ThreeVector PLAposition = thePrimaryGenerator->GetPLAPosition()/PlanetManager::GetInstance()->GetRplanet();
  G4ThreeVector PLAdirection = thePrimaryGenerator->GetPLADirection();
  
  
  // G4ThreeVector position=SpaceCoordinatePlanet::GetInstance()
  //                                 ->Transform(PLAposition,"PLA","CoordSys")/PlanetManager::GetInstance()->GetRplanet();
                   
  theAsciiFile<<"PLA Position:"<<'\t'<<PLAposition.x()<<'\t'
                                     <<PLAposition.y()<<'\t'
				     <<PLAposition.z()<<G4endl;
  theAsciiFile<<"PLA Direction:"<<'\t'<<PLAdirection.x()<<'\t'
                                     <<PLAdirection.y()<<'\t'
				   <<PLAdirection.z()<<G4endl;
   
  theAsciiFile<<"Time"<<'\t'<<"Ru"<<'\t'
                    <<"Rc"<<'\t'<<"Rl"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::RegisterCutoffVsTime
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double time)
{ theAsciiFile<<time/s<<'\t'<<Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MagneticShieldingToolAnalysisManager::CloseAsciiFile()
{ theAsciiFile.close();
}
