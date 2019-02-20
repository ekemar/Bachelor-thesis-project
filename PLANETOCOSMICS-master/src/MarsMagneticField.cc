#include "MarsMagneticField.hh"
#include "PlanetMagneticField.hh"
#include "MarsMagneticFieldMessenger.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4MagIntegratorStepper.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "SpaceCoordinatePlanet.hh"
#include "PlanetUnits.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4UImanager.hh"
#include"G4RunManager.hh"
#include"BlineTool.hh"




////////////////////////////////////////////////////////////////////////////////
//
MarsMagneticField::MarsMagneticField(): PlanetMagneticField("Mars")
{ 
  ListOfInternalFieldModels.push_back("CAIN90");
  ListOfInternalFieldModels.push_back("PURUCKER");
  ListOfInternalFieldModels.push_back("CAIN50");
  
  
  HasAGlobalField =false;
  //messenger
  theFieldMessenger = new MarsMagneticFieldMessenger(this); 
  
  SetInternalField("CAIN90");
  Purucker_cut_r2 = 1500.*km*1500*km;
     

  External=false;
  Internal=true;
  
  
  
  
}
////////////////////////////////////////////////////////////////////////
//
MarsMagneticField::~MarsMagneticField()
{ 
}


////////////////////////////////////////////////////////////////////////////////
//
void MarsMagneticField::SetInternalField(  G4String aString)
{
  if (!SetMotherInternalField(aString)){
  	Internal=true;
   	if (aString == "CAIN90") {
		// G4cout<<"CAIN90"<<std::endl;
  		if (getenv("MARS_CRUSTAL")){
  			G4String name_file = getenv("MARS_CRUSTAL")
		                        	+G4String("/Cain_n90Coefficients.txt");
			ReadGaussCoefficients(name_file);
			PInternalField=&MarsMagneticField::GetCAIN;
		
  		}
		else {
			G4cout<<"You should define the environment variable MARS_CRUSTAL"<<std::endl;
			Internal=false;
		}
  	
  	}
  	else if (aString == "CAIN50"){
  		if (getenv("MARS_CRUSTAL")){
  			G4String name_file = getenv("MARS_CRUSTAL")+
		                          G4String("/Cainn50Coefficients.txt");
			ReadGaussCoefficients(name_file);
			PInternalField=&MarsMagneticField::GetCAIN;
		
  		}
		else {
			G4cout<<"You should define the environment variable MARS_CRUSTAL"<<std::endl;
			Internal=false;
		}
  	
  	} 
  	else if (aString == "PURUCKER"){
  		ReadPuruckerData();
  		PInternalField=&MarsMagneticField::GetPURUCKER;
	}
	else Internal=false;
	
  }
  else {
  	PInternalField=PMotherInternalField;
  }
    
}

///////////////////////////////////////////////////////////////////////////////
//
void MarsMagneticField::SetExternalField(  G4String )
{External=false;
 SetMagnetopauseModel("SPHERE");
  
}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MarsMagneticField::GetInternalField(G4ThreeVector aVector) const
{ return (this->*PInternalField)(aVector);  
}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MarsMagneticField::GetExternalField(G4ThreeVector aVector) const
{ return (this->*PExternalField)(aVector); 
}

////////////////////////////////////////////////////////////////////////////////
//
void  MarsMagneticField::SetMagnetopauseModel( G4String aString) 
{ Magnetopause=SetMotherMagnetopauseModel(aString); 
  POutsideMagnetosphere=PMotherOutsideMagnetosphere; 	    
}
////////////////////////////////////////////////////////////////////////////////
//
void  MarsMagneticField::ComputeBfieldParametersAtReferenceTime()
{;
}
////////////////////////////////////////////////////////////////////////////////
//
bool  MarsMagneticField::OutsideMagnetosphere( G4ThreeVector pos) const
{return  (this->*POutsideMagnetosphere)(pos);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MarsMagneticField::GetPURUCKER(G4ThreeVector pos) const
{G4ThreeVector Bfield = G4ThreeVector(0.,0.,0.); 
 for (unsigned int i=0; i<Purucker_dip_Pos.size()-1;i++){
	G4ThreeVector diff_r = pos - Purucker_dip_Pos[i];
	G4double r2 = diff_r.mag2();
	if (r2 <= Purucker_cut_r2){  
		G4ThreeVector M_vec_unit= Purucker_dip_Pos[i]/Purucker_dip_Pos[i].mag();
		G4double r = std::sqrt(r2);
		G4double r3 = r*r2;
		G4double f1 = Purucker_dip_M[i]/r3;
		G4double f2 = f1*3.*diff_r.dot(M_vec_unit)/r2;
		Bfield += f1*M_vec_unit -f2* diff_r;
		
	}	
 
	
 }
 
 
 return Bfield;
 
 
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MarsMagneticField::GetCAIN(G4ThreeVector pos) const
{//The position should be scaled because the Spherical coefficients of
 //CAIN model are defined by considering Rmars = 3390.*km 
 G4ThreeVector scaled_pos = pos*rplanet/ (3390.*km);
 return  GetGAUSS(scaled_pos);
 
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MarsMagneticField::ReadPuruckerData() 
{ G4String name_file = getenv("MARS_CRUSTAL")
		                        +G4String("/Purucker_radial_dipole.txt");
		

  std::fstream File_Input(name_file,  std::ios::in);
 
 //clear vectors
  Purucker_dip_Pos.clear();
  Purucker_dip_M.clear();
  G4double a,b,c;
  
  do{
  	File_Input>>a>>b>>c;
	G4double cos_theta = std::cos((90.-b)*degree); 
	G4double sin_theta = std::sin((90.-b)*degree);
	G4double cos_phi = std::cos(a*degree);
	G4double sin_phi = std::sin(a*degree);
	Purucker_dip_Pos.push_back( 3393.5*km*G4ThreeVector(sin_theta*cos_phi,
							 sin_theta*sin_phi,
							 cos_theta));
	G4double fac = nanotesla*100.0*59.2*59.2*1.89*1.89*40.0*km*km*km;
	Purucker_dip_M.push_back(c*fac);
 }while (!File_Input.eof());	
}

