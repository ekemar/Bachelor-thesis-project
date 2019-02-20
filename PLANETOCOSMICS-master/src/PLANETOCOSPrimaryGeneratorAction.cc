#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSPrimaryMessenger.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
//#include "MYGeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinatePlanet.hh"
#include "PlanetUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "PlanetMagneticField.hh"
#include "PLANETOCOSPrimaryHit.hh"
#include "PLANETOCOSSD.hh"
#include "PlanetManager.hh"
#include "G4SDManager.hh"
#include "myfunctions.hh"

#include "isotropic.hh"	//PVD
#include <TRandom3.h>		//PVD
#include <TH1.h>		//PVD
#include <TF1.h>		//PVD
#include <TGraph.h>		//PVD





#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Geantino.hh"
#include "G4PrimaryVertex.hh"


#include "G4RunManager.hh"




#include "PLANETOCOSGeometryConstruction.hh"
//#include "PLANETOCOSAnalysisManager.hh"
#include "G4UnitsTable.hh"

#ifdef USE_ROOT_SOURCE
#include "ROOTPartSourceGenerator.hh"
#endif



////////////////////////////////////////////////////////////////////////////////
//




PLANETOCOSPrimaryGeneratorAction::PLANETOCOSPrimaryGeneratorAction()
{ myParticleSource= new MYGeneralParticleSource();
  myMessenger = new PLANETOCOSPrimaryMessenger(this);
  
  pGeneratePrimaries= 
         &PLANETOCOSPrimaryGeneratorAction::GenerateStandardPrimaries;
  SetDefaultRigidityVector();
  


  // first position and direction 
  G4bool test;
  test = SetPositionAndDirection("PLAG", 20.*km, 7.98*degree, 46.55*degree, 0.,0.);
  myParticleSource->SetParticleDefinition(G4Proton::Proton());
   
  verbosity=0;
   
  DetectPrimary =false;
  nb_of_primaries_for_scaling_flux =0;
  
  //cutoff parameter			     
  consider_cutoff = false;
  cutoff_hasbeen_already_defined = false;
  cutoff_rigidity = -100.;
  min_cutoff=0.;
  max_cutoff=0.;

#ifdef USE_ROOT_SOURCE
	UseRootSource = false;
	theROOTSourceGenerator = ROOTPartSourceGenerator::GetInstance();
	var1_type ="NONE";
	var2_type ="NONE";
	
	
	
#endif 
 
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryGeneratorAction::~PLANETOCOSPrimaryGeneratorAction()
{ delete myParticleSource;}
/////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{  (this->*pGeneratePrimaries)(anEvent);
   
   
   if (verbosity >0) 
          PrintPrimaryVertexInformation(anEvent->GetPrimaryVertex());
}
/////////////////////////////////////////////////////////////////////////////
void PLANETOCOSPrimaryGeneratorAction::GenerateStandardPrimaries(G4Event* anEvent)
{  //temporary event
  G4Event* tempEvent = new G4Event();
  G4double energy, rigidity, azimuth, zenith;
  G4ThreeVector position; 
#ifdef USE_ROOT_SOURCE 
   ComputePrimariesFromROOTHisto(tempEvent,energy, rigidity,position,zenith,azimuth); 
#else
  ComputePrimaries(tempEvent,energy, rigidity,position,zenith,azimuth);
#endif  
  nb_of_primaries_for_scaling_flux++;
  G4double rc=ComputeCutoff(zenith,azimuth,position);
  
  if (consider_cutoff){
   	while (rigidity <rc){
        	ComputePrimaries(tempEvent, energy, rigidity,position,zenith,azimuth);
	 	rc = ComputeCutoff(zenith,azimuth,position);
	 	nb_of_primaries_for_scaling_flux++;
	}
	
  }
  //Put first primary on event
  G4ThreeVector pos = tempEvent->GetPrimaryVertex()->GetPosition();
  G4double t0=tempEvent->GetPrimaryVertex()->GetT0();
  G4PrimaryVertex* thePrimaryVertex =
  		new G4PrimaryVertex(pos.x(),pos.y(),pos.z(),t0);
  G4PrimaryParticle* thePrimaryparticle = new G4PrimaryParticle();
  G4PrimaryParticle* tempPrimary = tempEvent->GetPrimaryVertex()->GetPrimary();
  thePrimaryparticle->SetG4code(tempPrimary->GetG4code());
  G4ThreeVector P =  tempPrimary->GetMomentum();
  
  
  
  thePrimaryparticle->SetMomentum(P.x(),P.y(),P.z());
  thePrimaryVertex->SetPrimary(thePrimaryparticle);
  anEvent->AddPrimaryVertex(thePrimaryVertex);
  delete tempEvent;
  
  //primary hit that is used by the analysis manager for purpose of
  //histograms normalisation and also for registering of primaries
  
  G4PrimaryVertex* aVertex = anEvent->GetPrimaryVertex();
  G4double weight = aVertex->GetWeight();
  G4int PDGCode = aVertex->GetPrimary()->GetPDGcode();
  
  //PVD: get latorX and longorY from position like in PLANETOCOSSteppingAction
  
  G4double lat_or_y, long_or_x;
  
  const PLANETOCOSGeometryConstruction* theGeometry = dynamic_cast<const PLANETOCOSGeometryConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction()); 
  G4String geometry = theGeometry->GetGeometryType();  
  
  if (geometry == "SPHERICAL") 
  	{
      	lat_or_y = 90.*degree -position.theta();
        long_or_x = position.phi();
        if(long_or_x>180.*degree) long_or_x = long_or_x-360*degree; 
        }
  else 
  	{
       	lat_or_y =position.y();
       	long_or_x =position.x();
	}  
  
  // this test is due to a problem appaearing with g4.8.1
  if (weight <= 0) weight =1.;

  PLANETOCOSPrimaryHit* aHit =
	            new  PLANETOCOSPrimaryHit(weight, PDGCode,energy,
                                                     azimuth, zenith,
						     lat_or_y, long_or_x);
  PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                      (G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
  sd->SetPrimaryHit(aHit);
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::
      GeneratePrimariesForComputingCutOffRigidity(G4Event* anEvent)
{ 
  Rigidity=rigidity_values[rigidity_index];
  G4double P = Rigidity* myParticleSource
                           ->GetParticleDefinition()->GetPDGCharge();
  G4double E0 = myParticleSource->GetParticleDefinition()->GetPDGMass()/GeV;
  G4double Etot=sqrt(E0*E0+P*P)*GeV;
#ifdef USE_OLD_GPS  
  myParticleSource->SetEnergyDisType("Mono");
  myParticleSource->SetMonoEnergy(Etot-E0*GeV);
#else
  myParticleSource->GetCurrentSource()
  		  ->GetEneDist()
		  ->SetEnergyDisType("Mono");
  myParticleSource->GetCurrentSource()
  		  ->GetEneDist()
		  ->SetMonoEnergy(Etot-E0*GeV);
#endif  
  
  rigidity_index++;
  myParticleSource->GeneratePrimaryVertex(anEvent);
}
////////////////////////////////////////////////////////////////////////////
G4bool PLANETOCOSPrimaryGeneratorAction
       ::SetPositionAndDirection(const G4String CoordSys, 
                                 const G4double anAlt,const G4double aLong,
				 const G4double aLat,const G4double aZenith,
				 const G4double  anAzimuth) 
{ G4bool test = SetPosition(CoordSys,anAlt,aLong,aLat);
  if (test)  test = SetDirection(CoordSys,aZenith,anAzimuth);
  return test;  
}
//////////////////////////////////////////////////////////////////////////////
void PLANETOCOSPrimaryGeneratorAction::SetPositionAndDirection(const G4String CoordSys, 
                                 const G4ThreeVector position,
				 const G4ThreeVector direction) 
{ G4bool test = SetPosition(CoordSys,position);
  if (test) test = SetDirection(CoordSys,direction);  
}
////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSPrimaryGeneratorAction::SetPosition(const G4String CoordSys, 
                     const G4double anAlt,const G4double aLong,
                     const G4double aLat)
{ 
 SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();
 G4double rplanet=PlanetManager::GetInstance()->GetRplanet();
 if (CoordSys == "PLAG" || CoordSys == "GEODETIC"){
 	PLAGaltitude=anAlt;
   	PLAGlongitude=aLong;
   	PLAGlatitude=aLat;
   	theConvertor->ComputePLAPositionFromPLAG(anAlt,aLat,aLong,PLAPosition);
 }
 else{
 	G4ThreeVector position=G4ThreeVector(1.,0.,0.);
  	position.setPhi(aLong);
  	position.setTheta(90.*degree-aLat);
  	position.setMag(rplanet+anAlt);
  	PLAPosition=theConvertor->Transform(position,CoordSys,"PLA");
  	theConvertor->ComputePLAGCoordinatesFromPLAPosition
                     (PLAPosition,PLAGaltitude, PLAGlongitude, PLAGlatitude);
 }
 G4ThreeVector xrot = G4ThreeVector(1.,0.,0.);
 xrot.setRThetaPhi(1.,90*degree,aLong+90*degree);
 G4ThreeVector zrot = PLAPosition/PLAPosition.mag();
 G4ThreeVector yrot = zrot.cross(xrot);
#ifdef USE_OLD_GPS
 myParticleSource->SetCentreCoords(PLAPosition);
 myParticleSource->DefineAngRefAxes("angref1",xrot);
 myParticleSource->DefineAngRefAxes("angref2",yrot);
#else
 myParticleSource->GetCurrentSource()
  		 ->GetPosDist()
		 ->SetCentreCoords(PLAPosition);
 myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref1",xrot);
 myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref2",yrot);		  
#endif 
 return true; 		     
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSPrimaryGeneratorAction::SetPosition(const G4String CoordSys, 
                                               const G4ThreeVector aPosition)

{ SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();
 	    
  PLAPosition=theConvertor->Transform(aPosition,CoordSys,"PLA");
 
  theConvertor->ComputePLAGCoordinatesFromPLAPosition
                     (PLAPosition,PLAGaltitude, PLAGlongitude, PLAGlatitude);
  
  G4ThreeVector xrot = G4ThreeVector(1.,0.,0.);
  xrot.setRThetaPhi(1.,90*degree,PLAGlongitude+90*degree);
  G4ThreeVector zrot = PLAPosition/PLAPosition.mag();
  G4ThreeVector yrot = zrot.cross(xrot);
#ifdef USE_OLD_GPS
 myParticleSource->SetCentreCoords(PLAPosition);
 myParticleSource->DefineAngRefAxes("angref1",xrot);
 myParticleSource->DefineAngRefAxes("angref2",yrot);
#else
 myParticleSource->GetCurrentSource()
  		 ->GetPosDist()
		 ->SetCentreCoords(PLAPosition);
 myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref1",xrot);
 myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref2",yrot);		  
#endif
 
  return true; 		     
}
///////////////////////////////////////////////////////////////////////////////
//
/*G4bool PLANETOCOSPrimaryGeneratorAction::SetPositionOnDipoleMagneticShell
                     (const G4String Reference, const G4double L, const G4double
		     latitude, const G4double  longitude)
{
 if (Reference != "PM" && Reference != "PSMG")
         {G4cout<<Reference<<" is not a good reference"<<G4endl;
	  return false;}
 
 G4String system;
 system=Reference;
 
 SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();

 G4double coslat=std::cos(latitude);
 G4double r=L*rm*coslat*coslat;
 
 G4ThreeVector Position=r*G4ThreeVector(1.,0.,0.);
 Position.setPhi(longitude); 
 Position.setTheta(90.*degree -latitude);
  

 PLAPosition = theConvertor->Transform(Position,system,"PLA");
			  

// dipole Shift 

 G4ThreeVector DipoleShift = ((PLANETOCOSMagneticField*)
   G4TransportationManager::GetTransportationManager()
      ->GetFieldManager()->GetDetectorField())->GetDipoleShift();
 
 PLAPosition += DipoleShift;

 theConvertor->ComputePLAGCoordinatesFromPLAPosition
                     (PLAPosition,PLAGaltitude, PLAGlongitude, PLAGlatitude);
#ifdef USE_OLD_GPS		     
  myParticleSource->SetCentreCoords(PLAPosition);
#else
  myParticleSource->GetCurrentSource()
  		  ->GetPosDist()
		  ->SetCentreCoords(PLAPosition);
#endif 
 return true; 
 
}
*/		     
////////////////////////////////////////////////////////////////////
G4bool PLANETOCOSPrimaryGeneratorAction::SetDirection(const G4String CoordSys, 
                      const G4double aZenith,const G4double anAzimuth)
{ SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();
 
   if (CoordSys == "PLAG" || CoordSys == "GEODETIC"){
  	 theConvertor->ComputePLADirectionAndPositionFromPLAG
	            (PLAGaltitude, PLAGlatitude, PLAGlongitude,
		     aZenith, anAzimuth,
		     PLAPosition, PLADirection);
   }
   else {
  	G4ThreeVector position=theConvertor->Transform(PLAPosition,"PLA",CoordSys);
  	G4ThreeVector direction =G4ThreeVector(0.,0.,1.);
  	direction.setTheta(aZenith);
  	direction.setPhi(180.*degree-anAzimuth);
  	direction=-direction.rotateY(position.theta()).rotateZ(position.phi());
  	PLADirection=theConvertor->Transform(direction,CoordSys,"PLA");
   }
#ifdef USE_OLD_GPS 
   
   myParticleSource->SetParticleMomentumDirection(PLADirection);
#else
   myParticleSource->GetCurrentSource()
		 ->GetAngDist()
		 ->SetAngDistType("planar");
   myParticleSource->GetCurrentSource()
  		  ->GetAngDist()
		  ->SetParticleMomentumDirection(PLADirection);
#endif   
 

   return true; 		     
}
//position and direction definitiojn for flat geometry
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::SetPosition(const G4ThreeVector aPosition)
{ const PLANETOCOSGeometryConstruction* theGeometry =
	 dynamic_cast<const PLANETOCOSGeometryConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
  G4double zpos_at_zero_alt = theGeometry->GetZpositionAtZeroAltitude();

  G4ThreeVector position = aPosition+G4ThreeVector(0.,0.,zpos_at_zero_alt);
  //G4cout<<position/km<<std::endl;
  G4ThreeVector xrot = G4ThreeVector(1.,0.,0.);
  G4ThreeVector yrot = G4ThreeVector(0.,1.,0.);
#ifdef USE_OLD_GPS		     
  myParticleSource->SetCentreCoords(position);
#else
  myParticleSource->GetCurrentSource()
  		  ->GetPosDist()
		  ->SetCentreCoords(position);
  myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref1",xrot);
  myParticleSource->GetCurrentSource()
			->GetAngDist()
			->DefineAngRefAxes("angref2",yrot);
#endif
  
}
///////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::SetDirection(const G4ThreeVector aDirection)
{
#ifdef USE_OLD_GPS 
   myParticleSource->SetParticleMomentumDirection(aDirection);
#else
   myParticleSource->GetCurrentSource()
		 ->GetAngDist()
		 ->SetAngDistType("planar");
   myParticleSource->GetCurrentSource()
  		  ->GetAngDist()
		  ->SetParticleMomentumDirection(aDirection);
#endif 
 

}
///////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::SetDirection(const G4double aZenith, const G4double anAzimuth)
{ G4ThreeVector direction =G4ThreeVector(0.,0.,1.);
  direction.setTheta(aZenith);
  direction.setPhi(180.*degree-anAzimuth);
  direction=-direction;
#ifdef USE_OLD_GPS 
   myParticleSource->SetParticleMomentumDirection(direction);
#else
   myParticleSource->GetCurrentSource()
  		  ->GetAngDist()
		  ->SetParticleMomentumDirection(direction);
#endif  

}


///////////////////////////////////////////////////////////////////////////////
//
G4bool PLANETOCOSPrimaryGeneratorAction::SetDirection(const G4String CoordSys, 
                                               const G4ThreeVector aDirection)
{SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();
 PLADirection=theConvertor->Transform(aDirection,CoordSys,"PLA");
#ifdef USE_OLD_GPS 
  myParticleSource->SetParticleMomentumDirection(PLADirection);
#else
  myParticleSource->GetCurrentSource()
  		  ->GetAngDist()
		  ->SetParticleMomentumDirection(PLADirection);
#endif  
 return true;
}	
///////////////////////////////////////////////////////////////////////////////
//
G4double  PLANETOCOSPrimaryGeneratorAction::GetSourceAltitude()
{

#ifdef USE_OLD_GPS 
 G4ThreeVector pos = myParticleSource->GetCentreCoords();
#else
 G4ThreeVector pos = myParticleSource->GetCurrentSource()
 				     ->GetPosDist()
				     ->GetCentreCoords();
#endif 
 const PLANETOCOSGeometryConstruction* theGeometry =
	 dynamic_cast<const PLANETOCOSGeometryConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
  
 G4String geometry = theGeometry->GetGeometryType();
 G4double zposition_at_zero_altitude 
                              = theGeometry->GetZpositionAtZeroAltitude();
 G4double altitude;
 if (geometry =="EUCLIDEAN"){
  	altitude =pos.z()-zposition_at_zero_altitude;  
 }
 else {
  	altitude =pos.mag()-zposition_at_zero_altitude;
 }
 return altitude;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::SetRigidity(G4double aRigidity)
{ G4double E0 = myParticleSource->GetParticleDefinition()->GetPDGMass();
  G4double P = aRigidity* myParticleSource
                           ->GetParticleDefinition()->GetPDGCharge();
  G4double Etot=sqrt(E0*E0+P*P);
#ifdef USE_OLD_GPS 
  myParticleSource->SetEnergyDisType("Mono");
  myParticleSource->SetMonoEnergy(Etot-E0);
#else
  myParticleSource->GetCurrentSource()
 		 ->GetEneDist()
		 ->SetEnergyDisType("Mono");
  myParticleSource->GetCurrentSource()
 		 ->GetEneDist()
		 ->SetMonoEnergy(Etot-E0);
#endif 		      
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::AddValuesToRigidityVector
                                     (G4int nvalues,G4double val1, G4double step)
{G4double val2=val1+G4double(nvalues)*step;
 if (val1 >=0 && val2 >=0) 
   {G4double val_max= std::max(val1,val2);
    G4double val_min= std::min(val1,val2);
    step=std::abs(step);
    if (rigidity_values.empty())
     {for (G4int i=0; i<nvalues;i++)
               rigidity_values.push_back(val_max - step * i);
     }
    else if (val_min > rigidity_values[0])
     {for (G4int i=0; i<nvalues;i++)
               rigidity_values.insert(rigidity_values.begin(),val_max - step * i);
     }
    else if  (val_max < rigidity_values.back())
     {for (G4int i=0; i<nvalues;i++)
               rigidity_values.push_back(val_max - step * i);
     }
    else 
     {G4cout<<"error when adding values to rigidity vector"<<G4endl; 
     }
   }
 else 
   {G4cout<<"error when adding values to rigidity vector"<<G4endl;
   }   
}
void PLANETOCOSPrimaryGeneratorAction::SetDefaultRigidityVector()
{rigidity_values.clear();
 AddValuesToRigidityVector(100,20,-0.1);
 AddValuesToRigidityVector(900,10,-0.01);
}


void PLANETOCOSPrimaryGeneratorAction::SelectTypeOfPrimaries(G4String aString)
{if (aString == "RigidityFilter") 
      pGeneratePrimaries= 
           &PLANETOCOSPrimaryGeneratorAction::GeneratePrimariesForComputingCutOffRigidity;
 else pGeneratePrimaries= 
           &PLANETOCOSPrimaryGeneratorAction::GenerateStandardPrimaries;
}
void PLANETOCOSPrimaryGeneratorAction
                 ::PrintPrimaryVertexInformation(const G4PrimaryVertex* aPrimaryVertex)
{G4double rplanet=PlanetManager::GetInstance()->GetRplanet();
 G4ThreeVector pla_pos = aPrimaryVertex->GetPosition()/rplanet;
 const G4PrimaryParticle* aPrimary = aPrimaryVertex->GetPrimary();
 
 G4ThreeVector pla_dir = aPrimary->GetMomentum();
 G4double P = pla_dir.mag();
 pla_dir/=P;
 G4double charge = myParticleSource
                           ->GetParticleDefinition()->GetPDGCharge();
 G4double rigidity = std::abs(P/charge/GV);
 G4String particle_name =aPrimary->GetG4code()->GetParticleName();
 G4double E0=aPrimary->GetG4code()->GetPDGMass();
 G4double E= std::sqrt(E0*E0 + P*P);
 G4double Ekin= (E-E0)/GeV;
 
 G4cout.precision(5);
 //G4cout<<setiosflags(0x1000);
 G4cout<<"New primary "<<std::endl;
 
 G4cout<<"Particle, Energy, Rigidity : "<<particle_name<<", "
                                        <<Ekin<<", "   
                                        <<rigidity<<std::endl;
					
 G4cout<<"PLA position : X "<<pla_pos.x()<<", Y "
                            <<pla_pos.y()<<", Z "   
                            <<pla_pos.z()<<", theta "
			    <<pla_pos.theta()/degree<<", phi "
			    <<pla_pos.phi()/degree<<std::endl; 
 G4cout<<"PLA direction : X "<<pla_dir.x()<<", Y "
                             <<pla_dir.y()<<", Z "   
                             <<pla_dir.z()<<", theta "
			     <<pla_dir.theta()/degree<<", phi "
			     <<pla_dir.phi()/degree<<std::endl; 
 if (verbosity>1){
 	SpaceCoordinatePlanet* theConvertor=
                            SpaceCoordinatePlanet::GetInstance();
  	G4ThreeVector psm_pos = theConvertor->Transform(pla_pos,"PLA","PSM");
  	G4ThreeVector pmag_pos = theConvertor->Transform(pla_pos,"PLA","PMAG");
  	G4ThreeVector psmag_pos = theConvertor->Transform(pla_pos,"PLA","PSMAG");
  
  	G4ThreeVector psm_dir = theConvertor->Transform(pla_dir,"PLA","PSM");
  	G4ThreeVector pmag_dir = theConvertor->Transform(pla_dir,"PLA","PMAG");
  	G4ThreeVector psmag_dir = theConvertor->Transform(pla_dir,"PLA","PSMAG");
  	G4double altitude,longitude,latitude;
  	theConvertor->ComputePLAGCoordinatesFromPLAPosition(pla_pos*rplanet,
                                                       altitude, 
						       longitude, 
						       latitude); 
						       
  	altitude/=km;
	
	G4cout<<"PSM position : X "<<psm_pos.x()<<", Y "
                            <<psm_pos.y()<<", Z "   
                            <<psm_pos.z()<<", theta "
			    <<psm_pos.theta()/degree<<", phi "
			    <<psm_pos.phi()/degree<<std::endl; 
        G4cout<<"PSM direction : X "<<psm_dir.x()<<", Y "
                             <<psm_dir.y()<<", Z "   
                             <<psm_dir.z()<<", theta "
			     <<psm_dir.theta()/degree<<", phi "
			     <<psm_dir.phi()/degree<<std::endl; 	  
  
 	
          
	G4cout<<"PMAG position : X "<<pmag_pos.x()<<", Y "
                            <<pmag_pos.y()<<", Z "   
                            <<pmag_pos.z()<<", theta "
			    <<pmag_pos.theta()/degree<<", phi "
			    <<pmag_pos.phi()/degree<<std::endl; 
        G4cout<<"PMAG direction : X "<<pmag_dir.x()<<", Y "
                             <<pmag_dir.y()<<", Z "   
                             <<pmag_dir.z()<<", theta "
			     <<pmag_dir.theta()/degree<<", phi "
			     <<pmag_dir.phi()/degree<<std::endl; 
	
	G4cout<<"PSMAG position : X "<<psmag_pos.x()<<", Y "
                            <<psmag_pos.y()<<", Z "   
                            <<psmag_pos.z()<<", theta "
			    <<psmag_pos.theta()/degree<<", phi "
			    <<psmag_pos.phi()/degree<<std::endl; 
        G4cout<<"PSMAG direction : X "<<psmag_dir.x()<<", Y "
                             <<psmag_dir.y()<<", Z "   
                             <<psmag_dir.z()<<", theta "
			     <<psmag_dir.theta()/degree<<", phi "
			     <<psmag_dir.phi()/degree<<std::endl; 
	
	G4cout<<"PLAG psoition alt, lat, long : "<<altitude<<", "
                           <<latitude/degree<<", "   
                           <<longitude/degree<<std::endl;		   		   
			   		   
	
			   
	}		   
 //G4cout<<setiosflags(0x800);                     		    
 
 }

///////////////////////////////////////////
void PLANETOCOSPrimaryGeneratorAction::PrintBfieldAtPrimary()
{
#ifdef USE_OLD_GPS
 G4ThreeVector pos = myParticleSource->GetCentreCoords();
#else
 G4ThreeVector pos = myParticleSource->GetCurrentSource()
 				     ->GetPosDist()
				     ->GetCentreCoords();
#endif 
 const PlanetMagneticField* theField =
                    reinterpret_cast<const PlanetMagneticField*>
	                 (G4TransportationManager::GetTransportationManager()->
		                                GetFieldManager()->GetDetectorField()); 
 //G4cout<<theField<<std::endl;
 if (theField){
 	theField->PrintBfield(pos);
 }
 else {
 	G4cout<<"The field is not switch on or does not exist!!!"<<std::endl;
 }	 	
}



////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::ComputePrimaries
                                       (G4Event* anEvent, G4double& energy,
                                        G4double& rigidity, G4ThreeVector& position,
			                G4double& zenith, G4double& azimuth )
{ //check if anEvent has already a vertex if yes remove it
  
  G4PrimaryVertex* aPrimaryVertex;
 /* aPrimaryVertex=anEvent->GetPrimaryVertex();
  if (aPrimaryVertex) delete aPrimaryVertex;*/
  
  
  //compute a new vertex
  myParticleSource->GeneratePrimaryVertex(anEvent);
  aPrimaryVertex=anEvent->GetPrimaryVertex();
  G4PrimaryVertex* lastVertex = aPrimaryVertex->GetNext();
  
  if (lastVertex) {
  		G4PrimaryVertex* vertex = lastVertex;
		while (vertex) {
			lastVertex=vertex;
			vertex=lastVertex->GetNext();
		}
		G4PrimaryParticle* lastPrimary = lastVertex->GetPrimary();
		G4PrimaryParticle* aPrimary = aPrimaryVertex->GetPrimary();
		G4ThreeVector P = lastPrimary->GetMomentum();
		aPrimary->SetMomentum(P.x(),P.y(),P.z());
		G4ThreeVector pos = lastVertex->GetPosition();
		aPrimaryVertex->SetPosition(pos.x(),pos.y(),pos.z());
		aPrimaryVertex->SetT0(lastVertex->GetT0());
 }
 
 /*if (IsAPitchAngleDistribution) {
  	//Check if Bfield exist
	G4PrimaryParticle* aPrimary = aPrimaryVertex->GetPrimary();
	G4ThreeVector pos = aPrimaryVertex->GetPosition();
	const G4Field* theField = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->GetDetectorField();
	if (theField) {
		G4double vec_pos[7],B[3];
		vec_pos[0]=pos.x();
		vec_pos[1]=pos.y();
		vec_pos[2]=pos.z();
		theField->GetFieldValue(vec_pos,B);
		G4ThreeVector Bf = G4ThreeVector(B[0],B[1],B[2]);
		//pitch angle distribution should be  toward the atmosphere
		if ( myDetector->GetGeometryType() == "SPHERICAL"){
			if (Bf.dot(pos) <0) Bf=-Bf;
		} 
		else if (Bf.dot(G4ThreeVector(0.,0.,1.))<0.) Bf=-Bf;
		G4ThreeVector P = aPrimary->GetMomentum();
	        P.rotateY(Bf.theta()).rotateZ(Bf.phi());
		aPrimary->SetMomentum(P.x(),P.y(),P.z());
	}
	else {
		G4cout<<" The magnetic field has not been switch on!"<<std::endl; 
		G4cout<<"The pitch angle distribution will be unselected."<<std::endl;
		IsAPitchAngleDistribution = false;
	}		      
  	
  }*/
  
  
  //G4cout<<"Primary 1"<<std::endl;
  G4PrimaryParticle* aPrimary = aPrimaryVertex->GetPrimary();
   
  position=aPrimaryVertex->GetPosition();
   
  G4ThreeVector direction=aPrimary->GetMomentum();
  G4double P=direction.mag();
   
  G4double charge =aPrimary->GetG4code()->GetPDGCharge();
  rigidity = P/charge;
  
  // G4cout<<rigidity<<'\t'<<energy<<std::endl;
  G4double E=aPrimary->GetMass();
  G4double E0=aPrimary->GetG4code()->GetPDGMass();	
  E=std::sqrt(P*P+E0*E0); 
  energy = E-E0;
   
  G4ThreeVector direction1 =-direction;
  
  const PLANETOCOSGeometryConstruction* theGeometry =
	 dynamic_cast<const PLANETOCOSGeometryConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
  
  G4String geometry = theGeometry->GetGeometryType();
  if (geometry == "SPHERICAL") direction1=direction1.rotateZ(-position.phi())
				     			.rotateY(-position.theta());
  zenith =  direction1.theta();
  azimuth = 180.*degree -direction1.phi();
  // G4cout<<"Primary 2"<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::ReadCutOffVsDirection(G4String file_name)
{std::fstream File_Input(file_name,  std::ios::in);

 
 
  //read the file
  char ch;
  G4String first_word;
  std::stringstream* aline = new std::stringstream();
 
  G4int i=0;
  std::vector<double> zenith;
  std::vector<double> azim;
  std::vector<double> Rc;
  G4double a,b,c,d,e;
  while (!File_Input.eof()){ 
 	i++;
  	delete aline;
  	aline = new std::stringstream();
  	File_Input.get(ch);
  	*aline<<ch;
  	while (ch != '\n' && !File_Input.eof()){
    		File_Input.get(ch);
	 	*aline<<ch;
    	} 
  	if (i>5 && !File_Input.eof()) {
		*aline>>a>>b>>c>>d>>e;
    		zenith.push_back(a*degree);
    		azim.push_back(b*degree);
    		Rc.push_back(d*GV);
   		// G4cout<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<'\t'<<e<<std::endl;
   	}
  }
 
 
 //fill the vectors
 
  cutoff_theta.clear();
  cutoff_phi.clear();
  cutoff_rigidities.clear();
  cutoff_phi.push_back(azim[0]);
  for (unsigned int i=1; i<azim.size(); i++){
  	if (azim[i] != azim[0]) cutoff_phi.push_back(azim[i]);
  	else i=azim.size();
  } 
  for (unsigned int i=0; i<azim.size(); i+=cutoff_phi.size()){
  	cutoff_theta.push_back(zenith[i]);
  } 
  int index=0; 
  for (unsigned int i=0; i<cutoff_theta.size(); i++){
  	cutoff_rigidities.push_back(std::vector<G4double>());
   	cutoff_rigidities[i].clear();
   	for (unsigned int j=0; j<cutoff_phi.size(); j++){
    		// G4cout<<cutoff_theta[i]/degree<<'\t';
     		// G4cout<<cutoff_phi[j]/degree<<'\t';
     		cutoff_rigidities[i].push_back(Rc[index]);
     		// G4cout<<cutoff_rigidities[i][j]/GV<<std::endl;
     		index++;
    	}
  }
  
 cutoff_hasbeen_already_defined = true;
 consider_cutoff=true;
 cutoff_type= "DIRECTION";  
 

 return;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double PLANETOCOSPrimaryGeneratorAction::ComputeCutoff(G4double zenith_angle,
						        G4double azimuth_angle,
							G4ThreeVector )
{
  if (cutoff_type == "DIRECTION")
 	return myfunc::BilinearInterpolation(cutoff_theta,cutoff_phi,
                                           cutoff_rigidities, 
					   zenith_angle, azimuth_angle);
 
  return cutoff_rigidity;
}
///////////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::SetCutoffRigidity(G4double aVal) 
{ cutoff_rigidity = aVal;
  cutoff_hasbeen_already_defined = true;
  consider_cutoff=true;
  cutoff_type= "FIXED";
}
///////////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::SetConsiderCutoff(G4bool aBool) 
{if (aBool && !cutoff_hasbeen_already_defined){
   	G4cout<<"Cutoff rigidity  vs direction "<<std::endl;
    	G4cout<<" or constant for all direction should be defined first"<<std::endl; 
    	return;
 }
 else {
   	consider_cutoff = aBool;
 }  
}  
////////////////////////////////////////////////////////////////////////////
G4double PLANETOCOSPrimaryGeneratorAction
          ::ModulatedDifferentialGalacticFlux(G4ParticleDefinition* aParticle,
                                              G4double energy,
				              G4double phi_mod)
{//this gives one model of the omnidirectionnal modulated galactic flux
 // the model is based on the field force approximation developed in
 // Gleeson and Axford, Solar modulation of galactic cosmic rays, Astrophys.
 // Journal,154, p.1011, 1968
 // the flux is given in function of energy and not in function 
 // of energy per nuc




  //general terms
  //--------------
  G4double nuc =std::abs(aParticle->GetBaryonNumber());
  G4double Ze = std::abs(aParticle->GetPDGCharge());
  
  G4double E0= aParticle->GetPDGMass()/MeV;
  //kinetic energy at Earth
  G4double Ekin_1AU = energy/MeV;
  //kinetic energy of the particle in the local interstellar medium 
  G4double Ekin_LIS = Ekin_1AU+phi_mod*Ze;
  
 
  G4double JE=Ekin_1AU*(Ekin_1AU+2.*E0)/(Ekin_LIS)
                                      /(Ekin_LIS+2.*E0);
 
  //model depedent term
  G4double LIS=0.; //Local interstellar spectrum. 
                  //This component is also model dependent. 
		  //We take here the model derived by Garcia Munoz
 
  
  G4double Ekin_LIS_nuc = Ekin_LIS/nuc;
  if (aParticle == G4Proton::Proton()){ 
     	G4double Cp=1.244e06/cm2/s/MeV;
      	G4double x=Ekin_LIS_nuc+780.*std::exp(-2.5e-04*Ekin_LIS_nuc);
      	LIS =Cp*std::pow(x,-2.65);
  }
      
  else if (aParticle == G4Alpha::Alpha()){
     	G4double Cp=3.1415926*4.*1.8e08/m2/s/MeV;
  
      	G4double x=Ekin_LIS_nuc+660.*std::exp(-1.4e-04*Ekin_LIS_nuc);
      	LIS =Cp*std::pow(x,-2.77)/nuc;
  }
 
  JE*=LIS;
  return JE;      
}
////////////////////////////////////////////////////////////////////////////////
// 
void  PLANETOCOSPrimaryGeneratorAction
          ::SelectModulatedGalacticFlux(G4ParticleDefinition* aParticle,
	                                G4double Enuc_min,
				        G4double Enuc_max,
					G4double phi_mod)							 
{
  G4int nuc = std::abs(aParticle->GetBaryonNumber());
  G4double Emin=Enuc_min*nuc;
  G4double Emax=Enuc_max*nuc;
  
  if (Emin>=Emax){
   G4cout<<"Your selected maximum energy per nucleon is below the cutoff rigidity"<<std::endl;
   G4cout<<"The spectrum definition procedure will be interrupted"<<std::endl;
   return;
  }
 
 
 
  G4int nE=int(std::log10(Emax/Emin)*10);
  if (nE>254) nE=254; 
  if (nE<3) nE=3;
  std::vector<G4double> energy_vec;
  std::vector<G4double> omni_flux_vec;
  G4double energy;
  G4double omni_flux;
  if (nE >=0 && (aParticle ==G4Proton::Proton() ||
                                       aParticle == G4Alpha::Alpha())){
    	nE++;
   	//Set spectrum parameters
#ifdef USE_OLD_GPS	
    	myParticleSource->SetEnergyDisType("Arb");
    	myParticleSource->ReSetHist("arb"); //arbitry point wise spectrum
    	myParticleSource->SetParticleDefinition(aParticle);
#else
	/*
	myParticleSource->ClearAll();
	myParticleSource->AddaSource(1.);
	*/
	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->SetEnergyDisType("Arb");
    	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->ReSetHist("arb"); //arbitry point wise spectrum
    	myParticleSource->SetParticleDefinition(aParticle);
	
	
#endif
   	// myParticleSource->InputDifferentialSpectra(true); //differential spectrum
    	for (int i=0;i<nE;i++){  // define the spectrum point per point
     		energy = Emin*std::pow(Emax/Emin,double(i)/(nE-1));
      		omni_flux = ModulatedDifferentialGalacticFlux(aParticle,
                                                    		energy,
					            		phi_mod);
      		energy_vec.push_back(energy);
      		omni_flux_vec.push_back(omni_flux);
      		//G4cout<<energy<<'\t'<<omni_flux<<std::endl;
#ifdef USE_OLD_GPS		
      		myParticleSource->ArbEnergyHisto(G4ThreeVector(energy,
                                                    	       omni_flux,
						    	       1.));
#else
		myParticleSource->GetCurrentSource()
				->GetEneDist()
				->ArbEnergyHisto(G4ThreeVector(energy,
                                                    	       omni_flux,
						    	       1.));	
#endif 					    						    
     	}
	G4cout<<"Test3"<<std::endl;
#ifdef USE_OLD_GPS	
    	myParticleSource->ArbInterpolate("Log");
	myParticleSource->SetAngDistType("cos");
#else
	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->ArbInterpolate("Log");
	myParticleSource->GetCurrentSource()
			->GetAngDist()
			->SetAngDistType("cos");
#endif    
    	//integration of the flux for later normalisation
    	//-----------------------------------------------
    	G4double energy_threshold;
    	energy_threshold = energy_vec[0];
    	G4double e2 =energy_vec.back();
    
    	
    	G4double integral_flux
              = myfunc::IntegrationOfY_pow(energy_vec, omni_flux_vec,
                                                        energy_vec[0],
							e2-0.0001*MeV);
    
    	integral_flux /=4.;//incoming flux at the top of the atmosphere a cos law is
                       // considered
		       
	G4cout<<"Primary flux [#/cm2/s]: "<<integral_flux*cm2*s<<std::endl;	       
    	
    	scaling_flux = integral_flux;
    	PLANETOCOSAnalysisManager* theAnalysisManager =
                                      PLANETOCOSAnalysisManager::GetInstance();
    	theAnalysisManager->SetPrimaryIntegralFlux(integral_flux);
                    
  }
  else { 
  	G4cout<<"The galactic cosmic ray flux you have selected will not be defined"<<std::endl;
  } 
  energy_vec.clear();
  omni_flux_vec.clear(); 

}
////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::
         SelectModulatedGalacticFluxAtSolMin(G4ParticleDefinition* aParticle,
	                                     G4double Enuc_min,
				             G4double Enuc_max)
{G4cout<<"SolMin"<<std::endl; 
 SelectModulatedGalacticFlux(aParticle,Enuc_min,Enuc_max,400.);                             
}
////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::
             SelectModulatedGalacticFluxAtSolMax(G4ParticleDefinition* aParticle,
	                                         G4double Enuc_min,
				                 G4double Enuc_max)
{ SelectModulatedGalacticFlux(aParticle,Enuc_min,Enuc_max,1000.);                             
}
////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::
             SelectMeanModulatedGalacticFlux(G4ParticleDefinition* aParticle,
	                                 G4double Enuc_min,
				         G4double Enuc_max)
{ SelectModulatedGalacticFlux(aParticle,Enuc_min,Enuc_max,550.);                             
} 
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::ReadUserSpectrum(G4String file_name)
{ 
  std::fstream File_Input(file_name,  std::ios::in);
  G4double energy_unit = MeV;
  G4double flux_unit = 1./m2/s/sr/MeV; 
  G4ParticleDefinition* primary_particle =0;
  G4String energy_type;
  G4String flux_type;
  G4String interpolation_type="Spline";
 

  //read the first data line
  //if file start by \comments
  //the first line start after \data or \definition
 
  char ch;
  G4String first_word;
  std::stringstream* aline = new std::stringstream();
 
  ch='\n';
  // avoid blank lines
  while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
 

  while (ch != '\n') {
 	*aline << ch;
  	File_Input.get(ch);
  }
  *aline>>first_word;
 
  // read the section of comments if exist
  //-----------------------------------------
 
  if (first_word=="\\comments"){
  	while (first_word != "\\data" && 
	       first_word != "\\definition" && !File_Input.eof()){
     		delete aline;
      		aline = new std::stringstream();
      		File_Input.get(ch);
      		*aline<<ch;
      		while (ch != '\n' && !File_Input.eof()){
        		File_Input.get(ch);
	 		*aline<<ch;
        	}
      	*aline>>first_word;
  	}
     
     
   	if (File_Input.eof()){ 
   		std::cout<<"No spectrum data have been read"<<std::endl;
       		return;
     	}
 }
 
 //read the section of definition if exist
 //----------------------------------------
  const PLANETOCOSGeometryConstruction* theGeometry =
	 dynamic_cast<const PLANETOCOSGeometryConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());	
  G4double altitude = theGeometry->GetAtmosphereHmax()*1.000001;
  if (first_word=="\\definition"){
  	while (first_word != "\\data" && !File_Input.eof()){
     		delete aline;
      		aline = new std::stringstream();
      		File_Input.get(ch);
      		while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
      		*aline<<ch;
      		while (ch != '\n' && !File_Input.eof()){
        		File_Input.get(ch);
	 		*aline<<ch;
        	}
      		*aline>>first_word;
      		if (first_word.find("\\energy_unit") == 0 ||
          	    first_word.find("\\flux_unit") == 0 ||
	  	    first_word.find("\\particle") == 0 || 
	  	    first_word.find("\\interpolation") == 0 ||
		    first_word.find("\\altitude") == 0){
	 
	   		str_size last= first_word.find("{");
	    		G4String tag_name = first_word;
	    		tag_name = tag_name.remove(last);
	    		tag_name = tag_name.remove(0,1);
            		G4String arg_val = first_word.remove(0,last+1);
	    		if (arg_val.find("}")  == arg_val.size()-1){
	     			arg_val = arg_val.remove(arg_val.size()-1);
	      			G4cout<<"you have selected the following "<<tag_name<<" : ";
	      			G4cout<<arg_val<<std::endl;
	      			if (tag_name =="energy_unit"){
	        			G4double aVal =G4UnitDefinition::GetValueOf(arg_val);
	         			G4String category =G4UnitDefinition::GetCategory(arg_val);
		 			energy_unit = aVal;
		 			if (category == "Energy"){
		    				energy_type="Energy";
		    			}
		 			else if (category == "Electric potential"){
		    				energy_type="Rigidity";
		    			} 
		 			else if (category == "Energy per nucleon"){
		    				energy_type="E_per_nuc";
		    			}    
		 			else{ 
		    				G4cout<<"your selected  energy unit is not defined"<<std::endl;
		     				return;
		    			}    
				}
	      			else if (tag_name == "flux_unit"){
	       				G4double aVal =G4UnitDefinition::GetValueOf(arg_val);
	        			G4String category =G4UnitDefinition::GetCategory(arg_val);
	        			flux_unit = aVal;
					if (category ==	"Integral dir flux")
		                           			flux_type = "int_dir";
					else if  (category ==	"Integral omni flux") 
		                            			flux_type = "int_omni";
					else if  (category ==	"Differential dir flux") 
		                            			flux_type = "dif_dir";
					else if  (category ==	"Differential omni flux")
					    			flux_type = "dif_omni";
					else { 
		 				G4cout<<"your selected  flux unit is not defined"<<std::endl;
		  				return;
					}			    
	      			 }
	     		 	else if (tag_name == "particle"){ 	
	        			primary_particle =  
		     			G4ParticleTable::GetParticleTable()->FindParticle(arg_val);
		  			if (!primary_particle){
		     				G4cout<<"your selected particle is not defined"<<std::endl;
		      				return;
		     			} 
		 		}
	      			else if (tag_name == "interpolation"){ 	
	          			if (arg_val == "Log" ||
		       			    arg_val == "Exp" ||
		                            arg_val == "Spline" ||
		                            arg_val == "Lin"){
	 					interpolation_type=arg_val;
		     			} 
		 		}
				else if (tag_name == "altitude"){ 	
	          			std::stringstream astream;
					astream <<arg_val;
					astream >>altitude;
					altitude*=km;
		 		}	   
             		}
            	}
	    
    		if (File_Input.eof()){ 
     			std::cout<<"No spectrum data have been read"<<std::endl;
       			return;
     		}  
    	}
  }
     
  std::vector<double> energy;
  std::vector<double> flux;       
  G4int i=0; 
  G4double a,b,c;

 if (!primary_particle){
   	G4cout<<"You should define a primary particle"<<std::endl;
    	G4cout<<"The procedure is interrupted"<<std::endl;
    	return;
 }
 if (first_word=="\\data"){ 
  	do {
		delete aline;
      		aline = new std::stringstream();
      		ch='\n';
      		// avoid blank lines
      		while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
   		while (ch != '\n' && !File_Input.eof()){
       			*aline<<ch;
        		File_Input.get(ch);
       		}      
      		if (!File_Input.eof()){
        		if (flux_type.contains("int")){ //integral flux
	   			*aline>>a>>b>>c;
	     			if (i==0) {
					//flux.push_back(0.);
					//energy.push_back(a/1.0000000000001);
					flux.push_back(1e-50);
					energy.push_back(a);
				}	
	     			energy.push_back(b);
	     			flux.push_back(c);
	     			if (energy[i+1]<= energy[i]){ 
	      				G4cout<<"Your energy vector is not increasing monotically"<<std::endl;
	       				G4cout<<"The procedure is interrupted"<<std::endl;
	       				return; 
	      			}
	   		}
	 		else {
	   			*aline>>a>>b;
	    			energy.push_back(a);
	    			flux.push_back(b);
	    			if (i>0 && energy[i]<= energy[i-1]){
	     				G4cout<<"Your energy vector is not increasing monotically"<<std::endl;
	      				G4cout<<"The procedure is interrupted"<<std::endl;
	      				return;
	     			}
	   		}   

		}
      		i++;   
    	}
   	while (!File_Input.eof());
  }
  
 
   
  
//Define the spectrum with the general particle source  
// arrange the energy vector and flux vector
  G4double energy_factor = energy_unit;
  G4double flux_factor = flux_unit;
  G4double nuc =std::abs(primary_particle->GetBaryonNumber());
  G4double E0= primary_particle->GetPDGMass();
  G4double q= primary_particle->GetPDGCharge();
  if (energy_type == "E_per_nuc" && nuc>0){
        energy_factor *=nuc;
        if (flux_type.contains("dif")) flux_factor /=nuc;
  }
  
  if  (flux_type.contains("dir")) flux_factor *=4.*pi;
  
          
  for(unsigned int i=0;i< energy.size();i++){
     	energy[i]*=energy_factor;
        if (energy_type == "Rigidity"){
	  	G4double P = energy[i]*q;
	   	energy[i]=sqrt(E0*E0+ P*P) -E0;
	}
	flux[i] *= flux_factor;               
  }
          
   
  // define the flux to the source
   
  if (flux_type.contains("int")){
#ifdef USE_OLD_GPS
        myParticleSource->SetParticleDefinition(primary_particle);
     	myParticleSource->SetEnergyDisType("User");
      	myParticleSource->ReSetHist("energy");
#else   
        /*myParticleSource->ClearAll();
	myParticleSource->AddaSource(1.); */
	myParticleSource->SetParticleDefinition(primary_particle);
	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->SetEnergyDisType("User");
      	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->ReSetHist("energy");
#endif	
	G4double integral_flux =0.;
 
      	for (unsigned int i=0;i<energy.size();i++){
#ifdef USE_OLD_GPS 
        	myParticleSource->UserEnergyHisto(G4ThreeVector(energy[i],
	                                                 	flux[i],
							  	0));
#else
		myParticleSource->GetCurrentSource()
				->GetEneDist()
				->UserEnergyHisto(G4ThreeVector(energy[i],
	                                                 	flux[i],
							  	0));
#endif								
		integral_flux +=flux[i];					                    
	}
	
	integral_flux /=4.;//incoming flux at the top of the atmosphere 
	//G4cout<<integral_flux<<std::endl;
	PLANETOCOSAnalysisManager* theAnalysisManager =
                                      PLANETOCOSAnalysisManager::GetInstance();
        theAnalysisManager->SetPrimaryIntegralFlux(integral_flux); 
	
  }
  else {
	//redefine the limit according to cutoff rigidity min and max
       	G4double Emin=energy.front();
       	G4double Emax=energy.back();
	if (consider_cutoff) {
        	G4double q =primary_particle->GetPDGCharge();
         	G4double E0 =primary_particle->GetPDGMass();
         	G4double Pmin = min_cutoff*q;
         	Emin=std::max(std::sqrt(E0*E0+Pmin*Pmin)- E0,Emin);
	}
       	if (Emin>=Emax) {
          	G4cout<<"Your slected maximum energy per is below what is allowed by";
	   	G4cout<<"the cutoff rigidity"<<std::endl;
           	G4cout<<"The spectrum definition procedure will be interrupted"<<std::endl;
           	return;
   	}
       	if (Emin != energy.front()) {
        	G4double flux_min = myfunc::LogInterpolation(energy, flux,Emin);
	 	G4bool out_of_range;
	 	G4int index = myfunc::locate( energy,Emin,out_of_range);
	 	for (int i=0;i<=index;i++){
	   		energy.erase(energy.begin());
	     		flux.erase(flux.begin());
	   	}
		// energy.erase(energy.begin(),&energy[index+1]);
	 	energy.insert(energy.begin(),Emin);
		// flux.erase(flux.begin(),&flux[index+1]);
	 	flux.insert(flux.begin(),flux_min);
	}	 
      
     
   	//integrate the spectrum for normalisation
       	G4double energy_threshold;
        energy_threshold = energy[0];
        G4double e2 =energy.back();
    
        G4double integral_flux
              = myfunc::IntegrationOfY_pow(energy, flux, energy[0],e2-0.00000001*MeV);
 
        integral_flux /=4.;//incoming flux at the top of the atmosphere  
  	
    	PLANETOCOSAnalysisManager* theAnalysisManager =
                                      PLANETOCOSAnalysisManager::GetInstance();
        theAnalysisManager->SetPrimaryIntegralFlux(integral_flux);
      
        
        scaling_flux = integral_flux;
     
       //define the spectrum as arbitrary point wise spectrum
#ifdef USE_OLD_GPS
	myParticleSource->SetParticleDefinition(primary_particle);         
	myParticleSource->SetEnergyDisType("Arb");
        myParticleSource->ReSetHist("arb");
        myParticleSource->InputDifferentialSpectra(true);
#else   
        /*myParticleSource->ClearAll();
	myParticleSource->AddaSource(1.); */
	myParticleSource->SetParticleDefinition(primary_particle);  
	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->SetEnergyDisType("Arb");
        myParticleSource->GetCurrentSource()
			->GetEneDist()
			->ReSetHist("arb");
        myParticleSource->GetCurrentSource()
			->GetEneDist()
			->InputDifferentialSpectra(true);
#endif        
      
      
      
        for (unsigned int i=0;i<energy.size();i++) {
         	//G4cout<<i<<'\t'<<energy[i]<<'\t'<<flux[i]<<std::endl;
#ifdef USE_OLD_GPS 		
	  	myParticleSource->ArbEnergyHisto(G4ThreeVector(energy[i],
	                                                	flux[i],
							  	0));
#else
	  	myParticleSource->GetCurrentSource()
				->GetEneDist()
				->ArbEnergyHisto(G4ThreeVector(energy[i],
	                                                	flux[i],
							  	0));
#endif								
	}
      	
       	//G4cout<<interpolation_type<<std::endl;
#ifdef USE_OLD_GPS 	
       	myParticleSource->ArbInterpolate(interpolation_type);
#else
	myParticleSource->GetCurrentSource()
			->GetEneDist()
			->ArbInterpolate(interpolation_type);
#endif
 }
#ifdef USE_OLD_GPS 
 myParticleSource->SetAngDistType("cos");
#else
 myParticleSource->GetCurrentSource()
		 ->GetAngDist()
		 ->SetAngDistType("cos");
#endif
 G4String geometry = theGeometry->GetGeometryType();
 if (geometry == "SPHERICAL") SetPosition("PLA",altitude,0.,0.);
 else   SetPosition(G4ThreeVector(0.,0.,altitude));
}

///////////////////////////////////////////////////////////////////////////////

void PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution(
	const G4String particle_distribution,
	const G4String coord_sys,
	const G4String length_unit,
	const G4String angle_unit,
	const G4double planet_radius,
	const G4double atmo_height,
	const G4double det_height,
	const G4double low_theta, 
	const G4double up_theta, 
	const G4double low_phi, 
	const G4double up_phi,
	const G4int eventNumber
)
{
G4cout<<"PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution: Generate "<<eventNumber<<" particles with an isotropic distribution in "<<det_height<<length_unit<<" which will be started in "<<atmo_height<<length_unit<<" altitude with a range of latitude: ("<<low_theta<<","<< up_theta<<") "<<angle_unit<<" and longitude: ("<<low_phi<<","<< up_phi<<") "<<angle_unit<<" according to "<<particle_distribution<<std::endl;

//Read distribution

ifstream in(particle_distribution);

//If file exists, read fluxes

bool goodFile = in.good();

char type[100];

std::vector<double> energy;
std::vector<double> prob;

if(goodFile)
	{
	//read type

	in>>type;

	while(1) 
		{
		G4double buf_energy, buf_prob;
		in >> buf_energy >> buf_prob;
		if (!in.eof()) 
			{	
			energy.push_back(buf_energy);
			prob.push_back(buf_prob);
			}			
		else break;
		}
	}
else G4cout<<"PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution: No file with fluxes found!\n"<<std::endl;
in.close();

//-------------------------------------------------------

G4ParticleDefinition* primary_particle = G4ParticleTable::GetParticleTable()->FindParticle(type);
if (!primary_particle) G4cout<<"PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution: your selected particle is not defined"<<std::endl;

//-------------------------------------------------------

// Calculate starting positions of simulation in such a way that the final ditribution in det_height is isotropic.

G4double altitude = atmo_height * G4UnitDefinition::GetValueOf(length_unit);

bool goodAngle = true;

TRandom3 * random = new TRandom3();
random->SetSeed(0);

std::vector<double> randomLat, randomLong, randomAlpha, randomBeta;

if(angle_unit == "degree")
	{
	for(int E = 0; E < eventNumber; E++)
		{
		std::vector<double> lat_long_alpha_beta = isotropic::GetRandomStartingPosition(planet_radius, det_height, atmo_height, low_theta, up_theta, low_phi, up_phi, random);

		randomLat.push_back(lat_long_alpha_beta[0]);
		randomLong.push_back(lat_long_alpha_beta[1]);
		randomAlpha.push_back(lat_long_alpha_beta[2]);
		randomBeta.push_back(lat_long_alpha_beta[3]);
		}
	}
else 
	{
	G4cout<<"PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution: Angle unit is not degree!"<<std::endl;
	goodAngle = false;
	}
delete random;			random = 0;

//-------------------------------------------------------

if(goodFile && primary_particle && goodAngle)
	{
	//create histogram for random distribution
	
	TGraph randomGraph;
	for(int i = 0; i < int(energy.size()); i++) randomGraph.SetPoint(i,energy[i],prob[i]);
	
	int bins = 100;
	
	int random_bins = bins;
	double random_left_int = energy[0];
	double random_right_int = energy[energy.size()-1];

	double random_left = random_left_int;
	double random_right = random_right_int;
	
	double log_energy[bins+1];
	for(int i = 0; i <= bins; i++) log_energy[i] = pow(10,log10(random_left)+(log10(random_right) - log10(random_left))/double(bins)*double(i));

	TH1D * randomHisto = new TH1D("primary_randomHisto","",bins, log_energy);
	//TH1D * randomHisto = new TH1D("primary_randomHisto","",bins, random_left,random_right);
	for(int i = 1; i <= random_bins; i++) 
		{
		double temp_energy = randomHisto->GetBinCenter(i);
		double temp_prob = randomGraph.Eval(temp_energy);
		if(temp_prob < 0) temp_prob = 0;
		randomHisto->SetBinContent(i,temp_prob);
		}

	//Set source
	
	myParticleSource->SetParticleDefinition(primary_particle);
	
	time_t start = time(0); 
	
	srand(time(NULL));
	
	for(int E = 0; E < eventNumber; E++)
		{
		// Set random energy
		int randombin = 1 + rand() % bins;
		double particle_energy = gRandom->Uniform(randomHisto->GetBinLowEdge(randombin),randomHisto->GetBinLowEdge(randombin)+randomHisto->GetBinWidth(randombin)) * 1e3;

		myParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(particle_energy);

		randomLat[E] *= G4UnitDefinition::GetValueOf(angle_unit);
		randomLong[E] *= G4UnitDefinition::GetValueOf(angle_unit);

		randomAlpha[E] *= G4UnitDefinition::GetValueOf(angle_unit); 
		randomBeta[E] *= G4UnitDefinition::GetValueOf(angle_unit);

		SetPosition(coord_sys,altitude,randomLong[E],randomLat[E]);
		SetDirection(coord_sys,randomAlpha[E],randomBeta[E]);

		G4RunManager::GetRunManager()->BeamOn(1);
		}
	time_t stop = time(0);
	G4cout<<"process time: "<<stop-start<<std::endl;
	}
else G4cout<<"PLANETOCOSPrimaryGeneratorAction::RandomIsotropicDistribution: No events will be processed!"<<std::endl;
}

#ifdef USE_ROOT_SOURCE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
	
void  PLANETOCOSPrimaryGeneratorAction::SetVar1Type(G4String aString)
{ if (aString=="energy" || aString=="theta" || aString=="cos_theta"   ) {
  	var1_type=aString;
  } 
  else 	var1_type="NONE"; 
;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void  PLANETOCOSPrimaryGeneratorAction::SetVar2Type(G4String aString)
{ if (aString=="energy" || aString=="theta" || aString=="cos_theta" ) {
  	var2_type=aString;
  } 
  else 	var2_type="NONE"; 
;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
	
void PLANETOCOSPrimaryGeneratorAction::ReadROOTFirstSourceHisto(G4String file_name,G4String path)
{
 theROOTSourceGenerator->ReadFirstSourceHisto(file_name, path);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
	
void PLANETOCOSPrimaryGeneratorAction::ReadROOTSecondSourceHisto(G4String file_name,G4String path)
{
 theROOTSourceGenerator->ReadSecondSourceHisto(file_name, path);
}


////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSPrimaryGeneratorAction::ComputePrimariesFromROOTHisto
                                       (G4Event* anEvent, G4double& energy,
                                        G4double& rigidity, G4ThreeVector& position,
			                G4double& zenith, G4double& azimuth )
{ 
  
  if (UseRootSource){
  	std::vector<double> res = theROOTSourceGenerator->GeneratePrimary();
	
 	if (res.size()>=1){
  		if (var1_type=="energy"){
			myParticleSource->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			myParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(res[0]);
		
		}
		else if (var1_type=="theta" || var1_type=="cos_theta"){
			myParticleSource->GetCurrentSource() ->GetAngDist()->SetAngDistType("iso");
			G4double theta =res[0];
			if (var1_type=="cos_theta") theta = std::acos(res[0]);
			myParticleSource->GetCurrentSource()->GetAngDist()->SetMinTheta(0.9999999*theta);
			myParticleSource->GetCurrentSource()->GetAngDist()->SetMaxTheta(theta);
		}
		if (res.size()==2 && var1_type != var2_type){
			if (var2_type=="energy"){
				myParticleSource->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
				myParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(res[1]);
		
			}
			else if (var2_type=="theta" || var2_type=="cos_theta"){
				myParticleSource->GetCurrentSource() ->GetAngDist()->SetAngDistType("iso");
   				G4double theta =res[1];
				if (var2_type=="cos_theta") theta = std::acos(res[1]);
				myParticleSource->GetCurrentSource()->GetAngDist()->SetMinTheta(0.9999999*theta);
				myParticleSource->GetCurrentSource()->GetAngDist()->SetMaxTheta(theta);
			}
		  		
  		}
  	}
 }	
  
  ComputePrimaries (anEvent,energy,rigidity, position, zenith, azimuth );
  
}
#endif 
