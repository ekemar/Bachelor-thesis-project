#include "PlanetManager.hh"
#include "PlanetMessenger.hh"

//PlanetUnits
#include "PlanetUnits.hh"

//Atmospheric model
#include "PlanetAtmosphereModel.hh"
#include "EarthAtmosphereModel.hh"
#include "MarsAtmosphereModel.hh"
#include "MercuryAtmosphereModel.hh"
#ifdef USE_JUPITER
#include "JupiterAtmosphereModel.hh"
#endif



//Magnetic field
#include "PlanetMagneticField.hh"
#include "MarsMagneticField.hh"
#include "EarthMagneticField.hh"
#include "MercuryMagneticField.hh"
#ifdef USE_JUPITER
#include "JupiterMagneticField.hh"
#endif

//Soil
#include "PlanetSoil.hh"

#include "SpaceCoordinatePlanet.hh"


//Unit
#include "G4UnitsTable.hh"
PlanetManager *PlanetManager::instance  = 0;

////////////////////////////////////////////////////////////////////////////////
//
PlanetManager::PlanetManager(G4String aName){
	theAtmosphereModel=0;
	theMagneticField=0;
	theSoil = new PlanetSoil();
	instance = this;
	theMessenger = new PlanetMessenger(this);
	if (aName !="") SelectPlanet(aName);
}
////////////////////////////////////////////////////////////////////////////////
//
PlanetManager::~PlanetManager()
{ //if (theAtmosphereModel) delete  theAtmosphereModel;
  //if (theMagneticField) delete theMagneticField;
  delete theSoil;
}
///////////////////////////////////////////////////////////////////////////////
//
PlanetManager *PlanetManager::GetInstance()
{
  if (instance == 0) instance = new PlanetManager();
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetManager::SelectPlanet(G4String aPlanetName)
{       PlanetName="UNKNOWN";
	if (theAtmosphereModel) delete theAtmosphereModel;
	if (theMagneticField) delete theMagneticField;
	if (aPlanetName =="Mars") {
		gplanet= gmars;
		rplanet=rmars;
		Rplanet=rplanet;
		BuildUnitTable();
		theAtmosphereModel = new MarsAtmosphereModel();
		theMagneticField = new MarsMagneticField();
		PlanetName="Mars";
		rplanet_eq = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtEquator();
		rplanet_pole = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtPole();
		SpaceCoordinatePlanet::GetInstance()->SetPlanetName(PlanetName);
		
		return;
	}
	else if (aPlanetName =="Earth"){
		gplanet= gearth;
		rplanet=re;
		Rplanet=rplanet;
		BuildUnitTable();
		theAtmosphereModel = new EarthAtmosphereModel();
		theMagneticField = new EarthMagneticField();
		PlanetName="Earth";
		SpaceCoordinatePlanet::GetInstance()->SetPlanetName(PlanetName);
		rplanet_eq = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtEquator();
		rplanet_pole = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtPole();
		return;
	}
	else if (aPlanetName =="Mercury"){
		gplanet= gmerc;
		rplanet=rmerc;
		Rplanet=rplanet;
		BuildUnitTable();
		theAtmosphereModel = new MercuryAtmosphereModel();
		theMagneticField =  new MercuryMagneticField();
		PlanetName="Mercury";
		SpaceCoordinatePlanet::GetInstance()->SetPlanetName(PlanetName);
		rplanet_eq = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtEquator();
		rplanet_pole = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtPole();
		return;
	}
#ifdef USE_JUPITER
	else if (aPlanetName =="Jupiter"){
		gplanet= gjupiter;
		rplanet=rjupiter;
		Rplanet=rplanet;
		BuildUnitTable();
		theAtmosphereModel = new JupiterAtmosphereModel();
		theMagneticField =  new JupiterMagneticField();
		PlanetName="Jupiter";
		SpaceCoordinatePlanet::GetInstance()->SetPlanetName(PlanetName);
		rplanet_eq = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtEquator();
		rplanet_pole = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtPole();
		return;
	}
#endif	
	G4cout<<"The planet that you have specified is not considered ";
	G4cout<<" in planetocosmics \n"<<std::endl;
	G4cout<<"Our beautiful Earth will be selected ";
	theAtmosphereModel = new EarthAtmosphereModel();
	theMagneticField = 0;
	PlanetName="Earth";
	gplanet= gearth;
	rplanet=re;
	BuildUnitTable();
	Rplanet=rplanet;
	
	//Space coordinate
	SpaceCoordinatePlanet::GetInstance()->SetPlanetName(PlanetName);
	rplanet_eq = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtEquator();
	rplanet_pole = SpaceCoordinatePlanet::GetInstance()->GetRadiusAtPole();
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetManager::BuildUnitTable()
{     
 // Definition of new units in the unit table should be defined
 // at the beginning before the
 // instantiation of any maessenger and should be followed by
 // G4UnitDefinition::BuildUnitsTable();

  	new G4UnitDefinition("nanotesla","nT","Magnetic flux density",nT);
 	new G4UnitDefinition("hour","hour","Time",3600.*s);
  	new G4UnitDefinition("minute","minute","Time",60.*s);
  	new G4UnitDefinition("day","day","Time",24.*3600.*s);
  	new G4UnitDefinition("gigavolt","GV","Electric potential",1000.*megavolt);
  	new G4UnitDefinition("1/cm3","1/cm3","Number density",1./cm3);
  	new G4UnitDefinition("1/m3","1/m3","Number density",1./m3);
  	new G4UnitDefinition("/cm3","/cm3","Number density",1./cm3);
  	new G4UnitDefinition("/m3","/m3","Number density",1./m3);
  	new G4UnitDefinition("#/cm3","#/cm3","Number density",1./cm3);
  	new G4UnitDefinition("#/m3","#/m3","Number density",1./m3);
  	new G4UnitDefinition("mole/cm3","mole/cm3","Number density",mole/cm3);
  	new G4UnitDefinition("mole/m3","mole/m3","Number density",mole/m3);
  	new G4UnitDefinition("hectopascal","hPa","Pressure",100.*pascal);
  	new G4UnitDefinition("g/cm2","g/cm2","Depth",g/cm2);
  	new G4UnitDefinition("g/m2","g/m2","Depth",g/m2);
  	new G4UnitDefinition("kg/cm2","kg/cm2","Depth",kg/cm2);
  	new G4UnitDefinition("kg/m2","kg/m2","Depth",kg/m2);
  
  
  	new G4UnitDefinition("PeV/nuc","PeV/n","Energy per nucleon",PeV);
  	new G4UnitDefinition("TeV/nuc","TeV/n","Energy per nucleon",TeV);
  	new G4UnitDefinition("MeV/nuc","MeV/n","Energy per nucleon",MeV);
  	new G4UnitDefinition("GeV/nuc","GeV/n","Energy per nucleon",GeV);
  	new G4UnitDefinition("keV/nuc","keV/n","Energy per nucleon",keV);
  
  	new G4UnitDefinition("#/m2/s/sr","#/m2/s/sr","Integral dir flux",1./m2/s/sr);
  	new G4UnitDefinition("1/m2/s/sr","1/m2/s/sr","Integral dir flux",1./m2/s/sr);
  	new G4UnitDefinition("#/cm2/s/sr","#/cm2/s/sr","Integral dir flux",1./cm2/s/sr); 
  	new G4UnitDefinition("1/cm2/s/sr","1/cm2/s/sr","Integral dir flux",1./cm2/s/sr); 
  
  	new G4UnitDefinition("#/m2/s","#/m2/s","Integral omni flux",1./m2/s);
  	new G4UnitDefinition("1/m2/s","1/m2/s","Integral omni flux",1./m2/s);
  	new G4UnitDefinition("#/cm2/s","#/cm2/s","Integral omni flux",1./cm2/s); 
  	new G4UnitDefinition("1/cm2/s","1/cm2/s","Integral omni flux",1./cm2/s); 
  
  	G4String type ="Differential dir flux";
  
  	G4String unit ="1/m2/s/sr/MeV";
  	G4String symbol ="#/m2/s/sr/MeV";
  	G4double val =1./m2/s/sr/MeV;
  	new G4UnitDefinition(unit,symbol,type,val);
 
  	unit ="1/cm2/s/sr/MeV";
  	symbol ="#/cm2/s/sr/MeV";
  	val =1./cm2/s/sr/MeV;
  	new G4UnitDefinition(unit,symbol,type,val);
  
  	unit ="1/m2/s/sr/GeV";
  	symbol ="#/m2/s/sr/GeV";
  	val =1./m2/s/sr/GeV;
  	new G4UnitDefinition(unit,symbol,type,val);
 
  	unit ="1/cm2/s/sr/GeV";
  	symbol ="#/cm2/s/sr/GeV";
  	val =1./cm2/s/sr/GeV;
  	new G4UnitDefinition(unit,symbol,type,val);
  
  	type ="Differential omni flux";
  
  	unit ="1/m2/s/MeV";
  	symbol ="#/m2/s/MeV";
  	val =1./m2/s/MeV;
  	new G4UnitDefinition(unit,symbol,type,val);
 
  	unit ="1/cm2/s/MeV";
  	symbol ="#/cm2/s/MeV";
  	val =1./cm2/s/MeV;
  	new G4UnitDefinition(unit,symbol,type,val);
  
  	unit ="1/m2/s/GeV";
  	symbol ="#/m2/s/GeV";
  	val =1./m2/s/GeV;
  	new G4UnitDefinition(unit,symbol,type,val);
 
  	unit ="1/cm2/s/GeV";
  	symbol ="#/cm2/s/GeV";
  	val =1./cm2/s/GeV;
  	new G4UnitDefinition(unit,symbol,type,val);
 
  
        type ="Length";
	unit ="Rplanet";
	symbol ="Rplanet";
	new G4UnitDefinition(unit,symbol,type,rplanet);
	
	unit ="rplanet";
	symbol ="rplanet";
	new G4UnitDefinition(unit,symbol,type,rplanet);
	
	G4UnitDefinition::BuildUnitsTable();
	
}
