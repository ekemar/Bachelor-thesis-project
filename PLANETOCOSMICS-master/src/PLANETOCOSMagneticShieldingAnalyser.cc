#include "PLANETOCOSMagneticShieldingAnalyser.hh"
#include "G4RunManager.hh"
#include "PlanetManager.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "SpaceCoordinatePlanet.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSMagneticShieldingAnalyser::PLANETOCOSMagneticShieldingAnalyser()
{
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSMagneticShieldingAnalyser::~PLANETOCOSMagneticShieldingAnalyser()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::OpenAsymptoticDirectionFile(G4String fileName)
{ // Asymptotic direction file
  theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"Rigidity"<<'\t'<<"Filter"<<'\t'<<"Asympt. Lat."<<'\t'
              <<"Asympt. Long."<<'\t'
              <<"Position Xpla"<<'\t'<<"Ypla"<<'\t'<<"Zpla"<<G4endl;
		    
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue)
{ //Rigidity 
    
  G4double rigidity =LastMomentaOnCharge.mag()/1000.;
    
  //Asymptotic direction
    
  G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
  G4double Aslat=90.-(AsympDir.theta()/degree);
  G4double Aslong=AsympDir.phi()/degree;
    
    
  //LastPosition
    
  G4ThreeVector LastPLAPosition;
  G4double rplanet=PlanetManager::GetInstance()->GetRplanet();
  LastPLAPosition = LastPosition/rplanet;
    
  //write on asymptotic direction file 
    
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);					   
  theAsciiFile<<rigidity<<'\t'<<'\t';
  theAsciiFile<<FilterValue<<'\t';
  theAsciiFile<<Aslat<<'\t'<<'\t'<<Aslong<<'\t'<<'\t'<<'\t';
  theAsciiFile<<LastPLAPosition.x()<<'\t'
                   <<LastPLAPosition.y()<<'\t'
		   <<LastPLAPosition.z()<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::CloseAsymptoticDirectionFile
                                    (G4double Rc,G4double Rm,G4double Rs)				
{//write rigidity cutoff
 
   theAsciiFile<<"Rs "<<Rs<<'\t';
   theAsciiFile<<"Rc "<<Rc<<'\t';
   theAsciiFile<<"Rm "<<Rm<<G4endl;
  
 //close file
  theAsciiFile.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::
            OpenCutoffVsPositionFile(G4String fileName,G4String CoordSys,
	                             G4double Altitude, G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Altitude"<<'\t'<<Altitude/km<<G4endl;
  theAsciiFile<<"Zenith"<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth"<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'<<"Rm"<<'\t'
                    <<"Rc"<<'\t'<<"Rs"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::RegisterCutoffVsPosition(G4double Rc,G4double Rm,G4double Rs,
                                     G4double latitude, G4double longitude)
{ theAsciiFile<<latitude/degree<<'\t'<<'\t'<<longitude/degree<<'\t'<<
                                            Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
 					       
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::
            OpenCutoffVsPositionOnLShellFile(G4String fileName,G4String CoordSys,
	                             G4double L, G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"L"<<'\t'<<L<<G4endl;
  theAsciiFile<<"Zenith"<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth"<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'<<"Rm"<<'\t'
                    <<"Rc"<<'\t'<<"Rs"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::
            OpenCutoffVsDirectionFile(G4String fileName,G4String CoordSys)
{ theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
 
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
 		(PLANETOCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  G4double rplanet=PlanetManager::GetInstance()->GetRplanet();
  if (CoordSys=="PLAG"){
  	G4double PLAGalt;
   	G4double PLAGlat;
   	G4double PLAGlong;
   	thePrimaryGenerator->GetPLAGPosition(PLAGalt,PLAGlat,PLAGlong);
   	theAsciiFile<<"Altitude:"<<'\t'<<PLAGalt/km<<" km"<<G4endl;
   	theAsciiFile<<"Latitude:"<<'\t'<<PLAGlat/degree<<" deg"<<G4endl;
   	theAsciiFile<<"Longitude:"<<'\t'<<PLAGlong/degree<<" deg"<<G4endl;  
  }
  else {
  	G4ThreeVector PLAposition = thePrimaryGenerator->GetPLAPosition();
    	G4ThreeVector position=SpaceCoordinatePlanet::GetInstance()
                                  	->Transform(PLAposition,"PLA","CoordSys")/rplanet;
                   
    	theAsciiFile<<"Position:"<<'\t'<<position.x()<<'\t'
                                   <<position.y()<<'\t'
				   <<position.z()<<G4endl;
    
  }
  theAsciiFile<<"Zenith"<<'\t'<<"Azimuth"<<'\t'<<"Rm"<<'\t'
                    <<"Rc"<<'\t'<<"Rs"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::RegisterCutoffVsDirection
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double zenith, G4double azimuth)
{theAsciiFile<<zenith/degree<<'\t'<<'\t'<<azimuth/degree<<'\t'<<
                                            Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
 					       
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::
            OpenCutoffVsTimeFile(G4String fileName)
{ theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
 
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
 	(PLANETOCOSPrimaryGeneratorAction*)
          	G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  G4double rplanet=PlanetManager::GetInstance()->GetRplanet();
  G4ThreeVector PLAposition = thePrimaryGenerator->GetPLAPosition()/rplanet;
  G4ThreeVector PLAdirection = thePrimaryGenerator->GetPLADirection();
  
  G4ThreeVector position=SpaceCoordinatePlanet::GetInstance()
                                  ->Transform(PLAposition,"PLA","CoordSys")/rplanet;
                   
  theAsciiFile<<"PLA Position:"<<'\t'<<PLAposition.x()<<'\t'
                                   <<PLAposition.y()<<'\t'
				   <<PLAposition.z()<<G4endl;
  theAsciiFile<<"PLA Direction:"<<'\t'<<PLAdirection.x()<<'\t'
                                   <<PLAdirection.y()<<'\t'
				   <<PLAdirection.z()<<G4endl;
   
  theAsciiFile<<"Time"<<'\t'<<"Rm"<<'\t'
                    <<"Rc"<<'\t'<<"Rs"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::RegisterCutoffVsTime
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double time)
{ theAsciiFile<<time/s<<'\t'<<Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSMagneticShieldingAnalyser::CloseAsciiFile()
{ theAsciiFile.close();
}
