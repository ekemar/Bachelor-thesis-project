#include "MercuryMagneticField.hh"
#include "PlanetMagneticField.hh"
#include "MercuryMagneticFieldMessenger.hh"
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
MercuryMagneticField::MercuryMagneticField(): PlanetMagneticField("Mercury")
{ 
  ListOfMagnetopauseModels.push_back("TSY96MERCURY");
  //messenger
  theFieldMessenger = new MercuryMagneticFieldMessenger(this); 
  HasAGlobalField =true;
  
  
  SetInternalField("DIPOLE");
  SetMagnetopauseModel("SPHERE");
  DipoleB0=300.*nT;
  SetDipoleAxis(0.,0.);
 
  External=false;
  Internal=true;
  
  
  
  
}
////////////////////////////////////////////////////////////////////////
//
MercuryMagneticField::~MercuryMagneticField()
{ 
}


////////////////////////////////////////////////////////////////////////////////
//
void MercuryMagneticField::SetInternalField(  G4String aString)
{ Internal= SetMotherInternalField(aString);
  PInternalField=PMotherInternalField;
}

///////////////////////////////////////////////////////////////////////////////
//
void MercuryMagneticField::SetExternalField(  G4String )
{ External=false;
  SetMagnetopauseModel("SPHERE");
  
}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MercuryMagneticField::GetInternalField(G4ThreeVector aVector) const
{ return (this->*PInternalField)(aVector);  
}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MercuryMagneticField::GetExternalField(G4ThreeVector aVector) const
{ return (this->*PExternalField)(aVector); 
}

////////////////////////////////////////////////////////////////////////////////
//
void  MercuryMagneticField::SetMagnetopauseModel( G4String aString) 
{ Magnetopause=SetMotherMagnetopauseModel(aString);  
  if (!Magnetopause){
        Magnetopause =true;
  	if (aString == "TSY96MERCURY")
     		POutsideMagnetosphere
                  	=&MercuryMagneticField::TSY96MercuryOutsideMagnetosphere; 		  
	else	Magnetopause =false;
 		
  }
  else 	POutsideMagnetosphere = PMotherOutsideMagnetosphere;  
}
////////////////////////////////////////////////////////////////////////////////
//
void  MercuryMagneticField::ComputeBfieldParametersAtReferenceTime()
{;
}
////////////////////////////////////////////////////////////////////////////////
//
bool  MercuryMagneticField::OutsideMagnetosphere( G4ThreeVector pos) const
{ if (Magnetopause) return  (this->*POutsideMagnetosphere)(pos);
  else return true;
}
////////////////////////////////////////////////////////////////////////////////
//
bool MercuryMagneticField::TSY96MercuryOutsideMagnetosphere(G4ThreeVector pos) const
{SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

 G4ThreeVector pos_pla=pos;
 // if (ConsiderDipoleShift) pos_pla-= DipoleSchift;  
   //PLAtoPSM
 G4ThreeVector pos_psm =theConvertor->Transform(pos_pla,"PLA","PSM");
  

 
 G4double xpsm=pos_psm.x()/rplanet;
 G4double ypsm=pos_psm.y()/rplanet;
 G4double zpsm=pos_psm.z()/rplanet;

 G4double XAPPA = 7.;
 G4double AM0, S0, X00, DSIG;
 AM0 = 70.;
 S0 = 1.08;
 X00 = 5.48;
 DSIG = 0.005;
 G4double X0=X00/XAPPA;
 G4double AM=AM0/XAPPA;
 G4double RHO2= ypsm*ypsm + zpsm*zpsm;
 G4double ASQ=AM * AM;
 G4double XMXM=AM+xpsm-X0;
 if (XMXM < 0.) XMXM=0.; 
 G4double AXX0 =XMXM * XMXM;
 G4double ARO=ASQ+RHO2 ;
 G4double SIGMA=std::sqrt((ARO+AXX0+ 
                            std::sqrt((ARO+AXX0)*(ARO+AXX0)-4.*ASQ*AXX0))
			                                            /(2.*ASQ));
								    
								    
								    
 /* if (!(SIGMA < S0+DSIG ) ) G4cout<<"Outside Tsy96 Magneto"<<std::endl; 
  else 	G4cout<<"Inside Tsy96 Magneto"<<pos_psm.mag()/rm<<std::endl; */					    
 return  (!(SIGMA < S0+DSIG ));
 
}
