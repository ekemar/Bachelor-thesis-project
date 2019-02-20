#include "PLANETOCOSSpenvisManager.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinatePlanet.hh"
#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "PlanetUnits.hh"

#include "PlanetMagneticField.hh"
#include "PlanetManager.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4GeneralParticleSource.hh"

#ifdef  USE_UNILIB
#include"unilib_c_pub.h"
#endif

PLANETOCOSSpenvisManager* PLANETOCOSSpenvisManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSpenvisManager::PLANETOCOSSpenvisManager()  
{ theSpenvisCSVCollection = 0;
  theSpenvisCSVFiles.clear();
  RegisterParticleTrajectory = false;
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSpenvisManager::~PLANETOCOSSpenvisManager() 
{ if (theSpenvisCSVCollection) delete  theSpenvisCSVCollection;
  
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSpenvisManager* PLANETOCOSSpenvisManager::GetInstance()
{
  if (instance == 0) instance = new PLANETOCOSSpenvisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::WriteData(std::vector< double> values,SpenvisCSV* aBlock)
{ //G4cout<<values.size()<<std::endl;
  if (!aBlock) theSpenvisCSVFiles[0].AddDataRow(values); 
  else aBlock->AddDataRow(values);
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::WriteBfieldInformation(SpenvisCSV* aBlock)
{ SpenvisCSV* theBlock;
  if (!aBlock) theBlock = &(theSpenvisCSVFiles[0]);
  else theBlock =aBlock;
  PlanetMagneticField* theField = PlanetManager::GetInstance()->GetMagneticField();

 
#ifdef  USE_UNILIB
  if ( theField->GetUseUnilibModel() && theField->Getul_ifail() >= 0){
 		theBlock->AddMetaVariableStr("BFIELD_MODEL", "UNILIB");
	
 
  }
  else {
#endif
	G4String Bintern = theField->GetIntFieldModelName();
		
	theBlock->AddMetaVariableStr("BFIELD_MODEL", "PLANETOCOS");
	//Date of reference for the magnetic field
	DateAndTime theDate = theField->GetReferenceDate();
	AddMetaVariableToCSVFile(G4String("YEAR"),
				 G4String(""),
				 theDate.year,
				 G4String("%4.0f"),theBlock);
	AddMetaVariableToCSVFile(G4String("MONTH"),
				 G4String(""),
				 theDate.month,
				 G4String("%2.0f"),theBlock);
	AddMetaVariableToCSVFile(G4String("DAY"),
				 G4String(""),
				 theDate.day,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("HOUR"),
				 G4String(""),
				 theDate.hour,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("MINUTE"),
			 	 G4String(""),
				 theDate.min,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("SECOND"),
					 G4String(""),
					 theDate.sec,
				         G4String("%2.0f"),
				         theBlock);
	AddMetaVariableToCSVFile(G4String("MSECOND"),
			         G4String(""),
			         theDate.msec,
				 G4String("%4.0f"),
				         theBlock);
		
		
	//Internal field			 		 			 			 				 
	theBlock->AddMetaVariableStr("BFIELD_INT", Bintern);
	/*if (Bintern == "IGRF"){
		AddMetaVariableToCSVFile(G4String("BFIELD_NMIGRF"),
					 G4String(""),
					 double(theField->Getnm_igrf()),
					 G4String("%2.0f"),
					 theBlock);
		}*/
	if (Bintern == "DIPOLE") {
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPB0"),
					 G4String("nT"),
					 theField->GetDipoleB0()/nanotesla,
					 G4String("%5.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPTH"),
					 G4String("degree"),
					 theField->GetDipoleTheta()/degree,
					 G4String("%4.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPPHI"),
					 G4String("degree"),
					 theField->GetDipolePhi()/degree,
					 G4String("%4.2f"),
					 theBlock);
		G4ThreeVector DipSchift =theField->GetDipoleShift();
		std::vector<double> dip_schift;
		dip_schift.push_back(DipSchift.x()/km);
		dip_schift.push_back(DipSchift.y()/km);
		dip_schift.push_back(DipSchift.z()/km);
		theBlock->AddMetaVariable(G4String("BFIELD_DIPSCHIFT"),
 					  dip_schift, 
					  "%6.2f", G4String("km"));				 
			
	}
	G4String Bextern = theField->GetExtFieldModelName();
	theBlock->AddMetaVariableStr("BFIELD_EXT", Bextern);
		
	//extern field parameter
	/*if (Bextern != "NOFIELD") {
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPTILT"),
					 G4String("degree"),
					 theField->GetDipolePS()/degree,
					 G4String("%4.2f"),
					 theBlock);
	}*/
	/*if (Bextern == "TSY89") {
		AddMetaVariableToCSVFile(G4String("BFIELD_KPIOPT"),
					 G4String(""),
					 theField->Getiopt(),
					 G4String("%1.0f"),
					 theBlock);
	}
	if (Bextern == "TSY96" || Bextern == "TSY2001" ) {
		AddMetaVariableToCSVFile(G4String("BFIELD_PDYN"),
					 G4String("Pa"),
					 theField->GetPdyn(),
					 G4String("%4.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DST"),
					 G4String("nT"),
					 theField->GetDst()/nanotesla,
					 G4String("%6.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_IMFY"),
					 G4String("nT"),
					 theField->GetImfy()/nanotesla,
					 G4String("%5.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_IMFZ"),
					 G4String("nT"),
					 theField->GetImfz()/nanotesla,
				 	 G4String("%5.2f"),
					 theBlock);
		if  (Bextern == "TSY2001") {
			AddMetaVariableToCSVFile(G4String("BFIELD_G1"),
					       	 G4String(""),
					       	 theField->GetG1(),
						 G4String("%6.2f"),
						 theBlock);
			AddMetaVariableToCSVFile(G4String("BFIELD_G2"),
					       	 G4String(""),
					       	 theField->GetG2(),
						 G4String("%6.2f"),
						 theBlock);		 
			
		}
						 
						 			 
	}*/
		
		
		
#ifdef  USE_UNILIB
 }
#endif

}

////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::SaveCSVFile(G4String name)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  for (unsigned int i=0; i<theSpenvisCSVFiles.size(); i++){
  	theSpenvisCSVCollection->AddCSVBlock(block_names[i],
						theSpenvisCSVFiles[i]);
	
  }
  theSpenvisCSVCollection->OutputCollection(name);
  for (unsigned int i=0; i<theSpenvisCSVFiles.size(); i++){
  	theSpenvisCSVCollection->DelCSVBlock(block_names[i]);
	
  }
/*  theSpenvisCSVFiles.clear();
  block_names.clear(); */
  delete theSpenvisCSVCollection;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::ClearCSVFiles()
{ theSpenvisCSVFiles.clear();
  block_names.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::InitialiseAsymptoticDirectionCSVFile()
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  
  //get a pointer to the primary generator action
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (PLANETOCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  G4double PLAGalt, PLAGlat, PLAGlong;
  // PLAGzenith, PLAGazimuth; 
  thePrimaryGenerator->GetPLAGPosition( PLAGalt,
  					 PLAGlat,
					 PLAGlong);
  					
  //thePrimaryGenerator->GetPLAGDirection( PLAGzenith, PLAGazimuth);  
  G4ThreeVector PLAPosition = thePrimaryGenerator->GetPLAPosition();
  std::vector<double> pos ;
  pos.push_back(PLAPosition.x()/km);
  pos.push_back(PLAPosition.y()/km);
  pos.push_back(PLAPosition.z()/km); 
  theSpenvisCSVFiles[0].AddMetaVariable("PLA_POSITION",pos,"%6.2f","km");
  pos.clear();
  
  G4ThreeVector PLADirection = thePrimaryGenerator->GetPLADirection();
  std::vector<double> dir ;
  dir.push_back(PLADirection.x());
  dir.push_back(PLADirection.y());
  dir.push_back(PLADirection.z());
  theSpenvisCSVFiles[0].AddMetaVariable("PLA_DIRECTION",dir,"%6.2f","");
  dir.clear();
  
  AddMetaVariableToCSVFile(G4String("PLAG_ALT"),
			   G4String("km"),
			   PLAGalt/km,
			   G4String("%6.2f"),
			   0);
  AddMetaVariableToCSVFile(G4String("PLAG_LAT"),
			   G4String("degree"),
			   PLAGlat/degree,
			   G4String("%6.2f"),
			   0);			   
  AddMetaVariableToCSVFile(G4String("PLAG_LON"),
			   G4String("degree"),
			   PLAGlong/degree,
			   G4String("%6.2f"),
			   0);
  /*AddMetaVariableToCSVFile(G4String("PLAG_ZEN"),
			   G4String("degree"),
			   PLAGzenith/degree,
			   G4String("%6.2f"),
			   0);			   
  AddMetaVariableToCSVFile(G4String("PLAG_AZIM"),
			   G4String("degree"),
			   PLAGazimuth/degree,
			   G4String("%6.2f"),
			   0);
  */
  //particle type
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
  //definition of variable
  
  theSpenvisCSVFiles[0].AddVariable("rigidity","GV",1,"rigidity","%6.3f");
  theSpenvisCSVFiles[0].AddVariable("trajectory label","",1,"trajectory label","%1.0f");
  theSpenvisCSVFiles[0].AddVariable("asymptotic latitude","degree",1,"asymptotic latitude","%4.2f");
  theSpenvisCSVFiles[0].AddVariable("asymptotic longitude","degree",1,"asymptotic longitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("last position in PLA coordinate ","Re",
  					3,"last position in PLA coordinate","%6.3f");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::RegisterAsymptoticDirection(G4ThreeVector LastPosition,
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
    LastPLAPosition = LastPosition/re;
    
    //give data
    std::vector<double > values;
    values.clear();
    values.push_back(rigidity);
    values.push_back(double(FilterValue));
    values.push_back(Aslat);
    values.push_back(Aslong);
    values.push_back(LastPLAPosition.x());
    values.push_back(LastPLAPosition.y());
    values.push_back(LastPLAPosition.z());
    WriteData(values);
    
   					    
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::InitialiseCutoffVsPosCSVFile(G4String coor_sys, 
   				     			G4double zenith,
				     			G4double azimuth)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("DIRECTION_ZEN"),
			   G4String("degree"),
			   zenith/degree,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_AZIM"),
			   G4String("degree"),
			   azimuth/degree,
			   G4String("%6.2f"),
			   0);
  			   
  
  
  //get a pointer to the primary generator action
  
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (PLANETOCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
  
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("SMJD","day",1,"Spenvis Modified julian date","%16.8f");
  theSpenvisCSVFiles[0].AddVariable("Altitude","km",1,"Altitude","%9.3f");
  theSpenvisCSVFiles[0].AddVariable("Latitude","degree",1,"Latitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Longitude","degree",1,"Longitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::InitialiseCutoffVsDirCSVFile(G4String coor_sys,
							G4double altitude,
				     			G4double latitude,
				     			G4double longitude 
   				     			)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("ALTITUDE"),
			   G4String("km"),
			   altitude/km,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("LATITUDE"),
			   G4String("degree"),
			   latitude/degree,
			   G4String("%6.2f"),
			   0);
   AddMetaVariableToCSVFile(G4String("LONGITUDE"),
			   G4String("degree"),
			   longitude/degree,
			   G4String("%6.2f"),
			   0);
  
  //get a pointer to the primary generator action
  
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (PLANETOCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
   			   			   
  			   
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("Zenith","degree",1,"Zenith angle of incoming direction","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Azimuth","degree",1,"Azimuth angle of incoming direction","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::InitialiseCutoffVsTimeCSVFile(G4String coor_sys,
							G4double altitude,
				     			G4double latitude,
				     			G4double longitude,
							G4double zenith,
							G4double azimuth 
   							)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("ALTITUDE"),
			   G4String("km"),
			   altitude/km,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("LATITUDE"),
			   G4String("degree"),
			   latitude/degree,
			   G4String("%6.2f"),
			   0);
  AddMetaVariableToCSVFile(G4String("LONGITUDE"),
			   G4String("degree"),
			   longitude/degree,
			   G4String("%6.2f"),
			   0);	
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_ZEN"),
			   G4String("degree"),
			   zenith/degree,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_AZIM"),
			   G4String("degree"),
			   azimuth/degree,
			   G4String("%6.2f"),
			   0);	
  
  			   		   			   
  			   		   
  //get a pointer to the primary generator action
  
  PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (PLANETOCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
   			   			   
  			   
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("SMJD","day",1,"Spenvis Modified julian date","%16.8f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::AddMetaVariableToCSVFile(G4String variable_name,
    				   G4String variable_unit,
				   G4double value, G4String format,
				    SpenvisCSV* aBlock)
{SpenvisCSV* theBlock;
 if (!aBlock) theBlock = &(theSpenvisCSVFiles[0]);
 else theBlock =aBlock;
 std::vector<double> values;
 values.push_back(value);
 theBlock->AddMetaVariable(variable_name,
 			   values, 
			   format, variable_unit);
 
}
/////////////////////////////////////////////////////////////////////////////
//
bool PLANETOCOSSpenvisManager::ReadSpenvisPositionGridCSVFile(G4String input_file,
   				   std::vector<double>& altitude,
   				   std::vector<double>& latitude,
				   std::vector<double>& longitude,
				   G4int nblock)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  theSpenvisCSVCollection->ReadCollection(input_file);
  G4String block_name;
  if (nblock <= 0) block_name="00001";
  else {
  	std::stringstream astream;
	astream <<nblock;
	astream >>block_name;
	if (nblock  <10) block_name="0000"+block_name;
	else if  (nblock  <100) block_name="000"+block_name;
	else if  (nblock  <1000) block_name="00"+block_name;
	else if  (nblock  <10000) block_name="0"+block_name;
	
  }
  SpenvisCSV aSpenvisCSV = theSpenvisCSVCollection->GetCSVBlock(block_name);
  delete theSpenvisCSVCollection;
  theSpenvisCSVCollection = 0;
  
  G4int nb_var = aSpenvisCSV.GetNumVariables();
  G4int nb_lines = aSpenvisCSV.GetNumDataLines();
  altitude.clear();
  latitude.clear();
  longitude.clear();
  
  if (nb_var && nb_lines >0){
  	for (int i=0; i< nb_lines;i++){
		G4double alt,lon,lat;
		alt = aSpenvisCSV.GetDataRecord("Altitude",i)[0]*km;
		lat = aSpenvisCSV.GetDataRecord("Latitude",i)[0]*degree;
		lon = aSpenvisCSV.GetDataRecord("Longitude",i)[0]*degree;
		altitude.push_back(alt);
		latitude.push_back(lat);
		longitude.push_back(lon);		
	}
  	
  		
	return true;
		
  }
  return false;
  
   
  
   
 
}				   
////////////////////////////////////////////////////////////////////////////////
//
bool PLANETOCOSSpenvisManager::ReadSpenvisTrajectoryCSVFile(G4String input_file,
   				   std::vector<double>& altitude,
   				   std::vector<double>& latitude,
				   std::vector<double>& longitude,
				   std::vector<double>& smjd,
				   G4int nblock)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  theSpenvisCSVCollection->ReadCollection(input_file);
  G4String block_name;
  if (nblock <= 0) block_name="00001";
  else {
  	std::stringstream astream;
	astream <<nblock;
	astream >>block_name;
	if (nblock  <10) block_name="0000"+block_name;
	else if  (nblock  <100) block_name="000"+block_name;
	else if  (nblock  <1000) block_name="00"+block_name;
	else if  (nblock  <10000) block_name="0"+block_name;
	
  }
  SpenvisCSV aSpenvisCSV = theSpenvisCSVCollection->GetCSVBlock(block_name);
  delete theSpenvisCSVCollection;
  theSpenvisCSVCollection = 0;
  
  G4int nb_var = aSpenvisCSV.GetNumVariables();
  G4int nb_lines = aSpenvisCSV.GetNumDataLines();
  altitude.clear();
  latitude.clear();
  longitude.clear(),
  smjd.clear();
  if (nb_var && nb_lines >0){
  	for (int i=0; i< nb_lines;i++){
		G4double alt,lon,lat,mjd;
		alt = aSpenvisCSV.GetDataRecord("Altitude",i)[0]*km;
		lat = aSpenvisCSV.GetDataRecord("Latitude",i)[0]*degree;
		lon = aSpenvisCSV.GetDataRecord("Longitude",i)[0]*degree;
		mjd = aSpenvisCSV.GetDataRecord("MJD",i)[0];
		altitude.push_back(alt);
		latitude.push_back(lat);
		longitude.push_back(lon);
		smjd.push_back(mjd);
		
	}
  	
  		
	return true;
		
  }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::InitialiseTrajectoryCSVBlock(G4String Type)
{ if (RegisterParticleTrajectory){
	G4bool BlineCase = false;
	if (Type.contains("Magnetic")) BlineCase = true;
	theParticleTrajectoryCSVBlocks.push_back(SpenvisCSV());
	SpenvisCSV* theBlock = &(theParticleTrajectoryCSVBlocks.back());
	WriteBfieldInformation(theBlock);
	theBlock->AddMetaVariableStr("TRAJECTORY_TYPE",Type);
	
	//Start position and direction  
	//get a pointer to the primary generator action
  	
	PLANETOCOSPrimaryGeneratorAction* thePrimaryGenerator =
   		(PLANETOCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  	G4double PLAGalt, PLAGlat, PLAGlong;
	// PLAGzenith, PLAGazimuth; 
  	thePrimaryGenerator->GetPLAGPosition( PLAGalt,
  					       PLAGlat,
					       PLAGlong);
  					
  //	thePrimaryGenerator->GetPLAGDirection( PLAGzenith, PLAGazimuth);  
  	G4ThreeVector PLAPosition = thePrimaryGenerator->GetPLAPosition();
  	std::vector<double> pos ;
  	pos.push_back(PLAPosition.x()/re);
  	pos.push_back(PLAPosition.y()/re);
  	pos.push_back(PLAPosition.z()/re); 
  	theBlock->AddMetaVariable("PLA_STARTPOSITION",pos,"%6.2f","re");
  	pos.clear();
  
 	G4ThreeVector PLADirection = thePrimaryGenerator->GetPLADirection();
	std::vector<double> dir ;
  	if (!BlineCase) {
  		dir.push_back(PLADirection.x());
  		dir.push_back(PLADirection.y());
  		dir.push_back(PLADirection.z());
  		theBlock->AddMetaVariable("PLA_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
	
	SpaceCoordinatePlanet* theConvertor = 
			SpaceCoordinatePlanet::GetInstance();
	G4ThreeVector Position = theConvertor->Transform(PLAPosition,"PLA","PSM");	
	G4ThreeVector Direction = theConvertor->Transform(PLADirection,"PLA","PSM");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("PSM_STARTPOSITION",pos,"%6.2f","re");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
  		theBlock->AddMetaVariable("PSM_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
	
	Position = theConvertor->Transform(PLAPosition,"PLA","PSO");	
	Direction = theConvertor->Transform(PLADirection,"PLA","PSO");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("PSO_STARTPOSITION",pos,"%6.2f","km");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
		theBlock->AddMetaVariable("PSO_STARTDIRECTION",dir,"%6.2f","");
		dir.clear();
	}	
  	
	
	Position = theConvertor->Transform(PLAPosition,"PLA","PMAG");	
	Direction = theConvertor->Transform(PLADirection,"PLA","PMAG");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("PMAG_STARTPOSITION",pos,"%6.2f","km");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
		theBlock->AddMetaVariable("PMAG_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
  	
	
	
			
  
  	AddMetaVariableToCSVFile(G4String("PLAG_ALT"),
			         G4String("km"),
			         PLAGalt/km,
			         G4String("%6.2f"),
			         theBlock);
        AddMetaVariableToCSVFile(G4String("PLAG_LAT"),
			         G4String("degree"),
			         PLAGlat/degree,
			         G4String("%6.2f"),
			         theBlock);			   
  	AddMetaVariableToCSVFile(G4String("PLAG_LON"),
			         G4String("degree"),
			         PLAGlong/degree,
			         G4String("%6.2f"),
			         theBlock);
    /*	AddMetaVariableToCSVFile(G4String("PLAG_ZEN"),
			         G4String("degree"),
			         PLAGzenith/degree,
			         G4String("%6.2f"),
			         theBlock);			   
  	AddMetaVariableToCSVFile(G4String("PLAG_AZIM"),
			         G4String("degree"),
			         PLAGazimuth/degree,
			         G4String("%6.2f"),
			         theBlock);
      */	
	// particle type, energy, rigidity
  
  	if (!BlineCase) {
		G4ParticleDefinition* theParticleDefinition = 
					thePrimaryGenerator
   						->GetParticleSource()
							->GetParticleDefinition();
		G4String particle_name= theParticleDefinition->GetParticleName();
  		theBlock->AddMetaVariableStr("PARTICLE",particle_name);
		
		G4double energy = thePrimaryGenerator
   						->GetParticleSource()
							->GetParticleEnergy();
		AddMetaVariableToCSVFile(G4String("ENERGY"),
			         G4String("GeV"),
			         energy/GeV,
			         G4String("%.3e"),
			         theBlock);
		
		G4double m0 = theParticleDefinition->GetPDGMass();
		G4double E=energy+ m0;
		G4double rigidity =std::abs( std::sqrt( E*E - m0*m0)
				 		/theParticleDefinition->GetPDGCharge());		 
		
		AddMetaVariableToCSVFile(G4String("RIGIDITY"),
			         G4String("GV"),
			         rigidity/GV,
			         G4String("%.3e"),
			         theBlock);		 
				 
	}
	
	//definition of variable
  
  	theBlock->AddVariable("PLA_POS","Re",3,"position in PLA cartesian coordinate","%6.3f");
  	theBlock->AddVariable("PMAG_POS","Re",3,"position in PMAG cartesian coordinate","%6.3f");
  	theBlock->AddVariable("PSM_POS","Re",3,"position in PSM cartesian coordinate","%6.3f");
  	theBlock->AddVariable("PSO_POS","Re",3,"position in PSO cartesian coordinate","%6.3f");
       
  
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::RegisterTrajectoryPoint(G4ThreeVector PLAPosition)
{ if (RegisterParticleTrajectory){
	SpenvisCSV* theBlock = &(theParticleTrajectoryCSVBlocks.back());
	std::vector<double> pos;
	pos.clear();
	SpaceCoordinatePlanet* theConvertor = 
			SpaceCoordinatePlanet::GetInstance();
	G4ThreeVector Position =PLAPosition/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	Position = theConvertor->Transform(PLAPosition,"PLA","PMAG")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	Position = theConvertor->Transform(PLAPosition,"PLA","PSM")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
		
	Position = theConvertor->Transform(PLAPosition,"PLA","PSO")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	WriteData(pos,theBlock);
	
		
  	
  }	
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSpenvisManager::SaveTrajectoryCSVBlocks(G4String name)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  for (unsigned int i=0; i<theParticleTrajectoryCSVBlocks.size(); i++){
  		G4String block_name;
		std::stringstream astream;
		astream<<i+1;
		astream >>block_name;
		if (i+1  <10) block_name="0000"+block_name;
		else if  (i+1  <100) block_name="000"+block_name;
		else if  (i+1  <1000) block_name="00"+block_name;
		else if  (i+1  <10000) block_name="0"+block_name;
		theSpenvisCSVCollection->AddCSVBlock(block_name,
						  theParticleTrajectoryCSVBlocks[i]);
	
  }
  theSpenvisCSVCollection->OutputCollection(name);
  delete theSpenvisCSVCollection;
}
