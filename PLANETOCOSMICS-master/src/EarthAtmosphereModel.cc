#include "EarthAtmosphereModel.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"

#include "G4PropagatorInField.hh"
//#include "IEarthAtmosphere.hh"
#include "EarthAtmosphereMessenger.hh"
#include "nrlmsise-00.hh"

EarthAtmosphereModel::EarthAtmosphereModel():
PlanetAtmosphereModel("Earth")
{ 
  
  ListOfAtmosphericModels.push_back("NRLMSISE00"); 
  ListOfAtmosphericModels.push_back("MSISE90");
  
  //Caution the messenger should be defined after the SetDefault()
  myMessenger = new EarthAtmosphereMessenger(this);
  
}
////////////////////////////////////////////////////////////////////////////////
//
EarthAtmosphereModel::~EarthAtmosphereModel()
{;
}


/////////////////////////////////////////
///////////////////////////////////////////
void EarthAtmosphereModel::
           ComputeAtmosphereTableFromModel(G4double alt_min,G4double alt_max)
{
  if (AtmosphericModel == "MSISE90" || AtmosphericModel == "NRLMSISE00"){
  	
	//clear vectors
   	//------------------
   	elements_in_atmosphere.clear();
   	atmosphere_mass_density_elements.clear();
   	atmosphere_altitude.clear();
   	atmosphere_depth.clear();
   	atmosphere_mass_density.clear();
   	atmosphere_number_density.clear();
   	atmosphere_temperature.clear();
   	atmosphere_pressure.clear();

       
   	O2_given=true;
   	N2_given=true;
	
	G4Element* elO;
        G4Element* elH;
        G4Element* elC;
        G4Element* elHe;
        G4Element* elAr;
        G4Element* elN; 
	
	elH = G4Element::GetElement("Hydrogen");
	elHe = G4Element::GetElement("Helium");
	elC = G4Element::GetElement("Carbon");
	elO = G4Element::GetElement("Oxygen");
	elAr = G4Element::GetElement("Argon");
	elN = G4Element::GetElement("Nitrogen");
   
   	elements_in_atmosphere.push_back(elH);
   	elements_in_atmosphere.push_back(elHe);
   	elements_in_atmosphere.push_back(elO);
   	elements_in_atmosphere.push_back(elN);
   	elements_in_atmosphere.push_back(elAr);
   
   	G4double avogadro =6.022e+23;
   	G4double mass_H  = elH->GetA()/avogadro;
   	G4double mass_He = elHe->GetA()/avogadro;
   	G4double mass_O  = elO->GetA()/avogadro;
   	G4double mass_N  = elN->GetA()/avogadro;
   	G4double mass_Ar = elAr->GetA()/avogadro;
  
  	// G4cout<<"test"<<std::endl;
   	for (unsigned int i=0; i<elements_in_atmosphere.size() ; i++)
      		atmosphere_mass_density_elements.push_back(std::vector<double> (atmosphere_altitude.size(),0.));
   
  	// G4cout<<"test1"<<std::endl;
   	// altitude vector 
   	G4double alt1,alt2;
   	alt1=alt_max;
   	alt2=alt_min;
   	if (alt_max > 500. *km) alt1=500.*km;
   	G4int nalt = int ((alt1-alt2)/(0.1*km));
   	atmosphere_altitude.push_back(alt1);
	
   	for (int i=0; i<nalt; i++)
        	atmosphere_altitude.push_back(alt1 - 0.1*km * (double) (i+1));
   
   	if (atmosphere_altitude[nalt]-alt2 >0.0001*km)
                           	atmosphere_altitude.push_back(alt2);
   			
   	//compute other vectors
   	//float d[8],dd[9];
   	//float t[2];
   	int yday=ReferenceDate.DayOfYear()+1;
   	float sec = (float) ReferenceDate.SecondFromDayStart();
   	float stl=sec/3600.+( (float) GEODETICLong)/15.;
   
   	G4double H_mdensity;
  	G4double He_mdensity;
   	G4double O_mdensity;
   	G4double N_mdensity;
   	G4double Ar_mdensity;
   	G4double ndensity;
  	G4double mdensity;
  	G4double temperature;
  	G4double pressure; 
   	G4double k_boltzmann =  1.381e-23*m*newton/kelvin;
	
	//C++ atmosphere - PVD
	nrlmsiseoutput atmo_output;
	nrlmsiseinput atmo_input;
  	nrlmsiseflags atmo_flags;
	aparray atmo_aph;
  	for (int i = 0; i < 7; i++) atmo_aph.a[i] = 100;
	atmo_flags.switches[0] = 0;
  	for (int i = 1; i < 24; i++) atmo_flags.switches[i] = 1;	
		
	atmo_input.doy = yday;
	atmo_input.year = 0; /* without effect */
  	atmo_input.sec = sec;
	atmo_input.g_lat = GEODETICLat;
	atmo_input.g_long = GEODETICLong;
	atmo_input.lst = stl;
	atmo_input.f107A = f107_a;
	atmo_input.f107 = f107;
	atmo_input.ap = Ap;
	
	nrlmsisecalc nrlmsise;
	
   	for (unsigned int i=0; i<atmosphere_altitude.size(); i++)
		{
         	if (AtmosphericModel == "NRLMSISE00" || AtmosphericModel == "MSISE90")
			{
			if (AtmosphericModel == "MSISE90") G4cout<<" USING NRLMSISE00 instead of MSISE90!"<<std::endl;
			
			float alt = (float) (atmosphere_altitude[i]/km); 
			if (alt == 32.5) alt = 32.499;
			atmo_input.alt=alt;
			
	  		//gtd7_(&yday, &sec, &alt, &GEODETICLat, &GEODETICLong, &stl, &f107_a, &f107, &Ap, &mm, dd, t);
			nrlmsise.gtd7(&atmo_input, &atmo_flags, &atmo_output);
			      
	  		H_mdensity=atmo_output.d[6]*mass_H/cm3;
	  		He_mdensity=atmo_output.d[0]*mass_He/cm3;
			//no anomalous oxygen taken into account as we are limited to <500km
	  		O_mdensity=(atmo_output.d[1]+2.*atmo_output.d[3])*mass_O/cm3;
	  		N_mdensity= (2.*atmo_output.d[2]+atmo_output.d[7])*mass_N/cm3;
	  		Ar_mdensity= atmo_output.d[4]*mass_Ar/cm3;
	  		ndensity = atmo_output.d[6] + atmo_output.d[0] + atmo_output.d[1] + atmo_output.d[2] + atmo_output.d[3] + atmo_output.d[4] + atmo_output.d[7];
	  		ndensity /= cm3;
	  		mdensity =  H_mdensity + He_mdensity + O_mdensity +  N_mdensity + Ar_mdensity;		
			}

		temperature=atmo_output.t[1]*kelvin;
	 	pressure = ndensity*k_boltzmann*temperature; //perfect gaz
	 	//G4cout<<" pressure "<<pressure/100./pascal<<std::endl; 
	 	//G4cout<<k_boltzmann<<std::endl;
	 	atmosphere_mass_density.push_back(mdensity);
	 	atmosphere_pressure.push_back(pressure);
	 	atmosphere_temperature.push_back(temperature);
	 	atmosphere_number_density.push_back(ndensity);
	 	atmosphere_mass_density_elements[0].push_back(H_mdensity);
	 	atmosphere_mass_density_elements[1].push_back(He_mdensity);
	 	atmosphere_mass_density_elements[2].push_back(O_mdensity);
	 	atmosphere_mass_density_elements[3].push_back(N_mdensity);
	 	atmosphere_mass_density_elements[4].push_back(Ar_mdensity);
     	}
		
  }
  
  //for all models the depth is computed in the same way
  
  atmosphere_depth.clear();
  atmosphere_depth.push_back(0.);


  for (unsigned int i=1; i<atmosphere_altitude.size() ; i++){
  	G4double dens1=atmosphere_mass_density[i];
        G4double dens0=atmosphere_mass_density[i-1];
	G4double h=atmosphere_altitude[i-1]-atmosphere_altitude[i];
	G4double depth_layer =0;
	if (std::abs(2.*(dens1-dens0)/(dens0+dens1))  <1e-12 ) depth_layer=(dens1+dens0)*h/2.;  
	else depth_layer=(dens1-dens0)*h/std::log(dens1/dens0);
	atmosphere_depth.push_back(atmosphere_depth[i-1]+depth_layer);
  } 
  WriteAtmosphereComposition("AtmoModel.txt"); 
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthAtmosphereModel::SetDefault()
{ 
  AtmosphericModel="NRLMSISE00";
  
  //Parameters for MSISE90 NRMLMSISE2001 models
  //-------------------------------------------- 
  ReferenceDate = DateAndTime(2000,1,1);
  GEODETICLong =0.;
  GEODETICLat =0.;
  Ap=4.;
  f107=150.;
  f107_a=150.;
}
