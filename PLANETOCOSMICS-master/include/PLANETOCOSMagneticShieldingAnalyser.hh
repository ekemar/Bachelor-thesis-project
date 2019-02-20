#ifndef PLANETOCOSMagneticShieldingAnalyser_HH
#define PLANETOCOSMagneticShieldingAnalyser_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include"vector"
#include"globals.hh"
#include"fstream"
#include"G4ThreeVector.hh"
class PLANETOCOSMagneticShieldingAnalyser
{  
public:
   
   PLANETOCOSMagneticShieldingAnalyser();
   ~PLANETOCOSMagneticShieldingAnalyser();
   
   
   //Asymptotic direction
   //--------------------
   void OpenAsymptoticDirectionFile(G4String fileName);
   void RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue);
   void CloseAsymptoticDirectionFile(G4double Rc,G4double Rm,G4double Rs);
   
   
   
   //Rigidity ForDiffrent position
   //-----------------------------
   void OpenCutoffVsPositionFile(G4String fileName,G4String CoordSys,
	                             G4double Altitude, G4double zenith,
				     G4double azimuth);
   void RegisterCutoffVsPosition(G4double Rc,G4double Rm,G4double Rs,
                                     G4double latitude, G4double longitude);
   
   //Rigidity For Different position on a magnetic dipole shell
   void OpenCutoffVsPositionOnLShellFile(G4String fileName,G4String CoordSys,
	                             G4double L, G4double zenith,
				     G4double azimuth);
  
    //Rigidity For Different directions
   void OpenCutoffVsDirectionFile(G4String fileName,G4String CoordSys);
   void RegisterCutoffVsDirection(G4double Rc,G4double Rm,G4double Rs,
                                     G4double zenith, G4double azimuth);
   
   
   
    //Rigidity For Different times
   void OpenCutoffVsTimeFile(G4String fileName);
   void RegisterCutoffVsTime(G4double Rc,G4double Rm,G4double Rs,
                                     G4double time);
   
   
  
   void CloseAsciiFile();

private:
   
   std::ofstream theAsciiFile;
};
#endif
