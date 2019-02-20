#include "SpaceCoordinatePlanet.hh"
#include "geomdefs.hh"
#include "time.h"
#include "PlanetUnits.hh"
#include "MarsData.hh"
#include "MercuryNewData.hh"
#include "EarthData.hh"
#include "SunData.hh"
#include "orbit_functions.hh"
#include "PlanetManager.hh"
#include "PlanetMagneticField.hh"

#ifdef USE_JUPITER
#include "JupiterData.hh"
#endif



#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include"vector"
#include"globals.hh"
#include"fstream"
#include <strstream>

extern "C" {
//#ifdef USE_SPICE
#include<SpiceCK.h>
#include<SpiceZpr.h>
//#endif
    double jyear_(void);
    double rpd_(void);
} 
SpaceCoordinatePlanet *SpaceCoordinatePlanet::instance = 0;

SpaceCoordinatePlanet::SpaceCoordinatePlanet()
{ // Set the start date
  ReferenceDate = DateAndTime(2000,1,1,12,0,0);
     
  verbose=0;
  //angle eccliptic and Earth equator at J2000 
  eps0_2000=23.439291111*degree;
  HAE2000toGEI2000=G4RotationMatrix(G4ThreeVector(1.,0.,0.),eps0_2000);
  GEI2000toHAE2000=HAE2000toGEI2000.inverse();
  SetPlanetName("Mars");
//#ifdef USE_SPICE	 
     UseSpice = true; 
//#endif
  Initialise();
}

////////////////////////////////////////////////////////////////////////////////
//
SpaceCoordinatePlanet::~SpaceCoordinatePlanet()
{
    ;
}
////////////////////////////////////////////////////////////////////////////////
//
SpaceCoordinatePlanet *SpaceCoordinatePlanet::GetInstance()
{
  if (instance == 0)
	instance = new SpaceCoordinatePlanet;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::SetSystemInAndOut(G4String sys_in, G4String sys_out) 
{ G4bool sys_in_ok =false;
  G4bool sys_out_ok =false;
  for (unsigned int i=0; i<ListOfCoordinateSystems.size();i++){
  	if (sys_in  == ListOfCoordinateSystems[i]) sys_in_ok=true;
  	if (sys_out  == ListOfCoordinateSystems[i]) sys_out_ok=true;
  }
  

  if (!sys_in_ok || !sys_out_ok ) {
	G4cout << sys_in<<" and/or "<< sys_out<< " are not good systems of coordinates" 
			<<std::endl;
	return;		
  }
  

  if (PlanetName =="Earth"){
        system_in = PlanetCoordinateSystemNameFromEarthCoordinateSystem(sys_in);
	system_out = PlanetCoordinateSystemNameFromEarthCoordinateSystem(sys_out);
  }
  else {
	system_in = sys_in;
	system_out = sys_out;
  }	
  
  ComputeSelectedMatrix();
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::SetPlanetName(G4String aName) 
{
 if (aName == "Mars" ){
 	 PlanetSpiceCode=499;
         PlanetIsMagnetic =false;
 	 PlanetName=aName;
	 a_planet = const_cast<double*>(a_mars);
         e_planet = const_cast<double*>(e_mars);
         lm_planet = const_cast<double*>(lm_mars);
         mu_planet = mu_mars;
	 Omega_planet = const_cast<double*>(Omega_mars);
         p_planet = const_cast<double*>(p_mars); //omega+Omega
         i_planet = const_cast<double*>(i_mars);
         pu_planet = const_cast<double*>( pu_mars);
 	 Ca_planet = const_cast<double*>(Ca_mars);
 	 Sa_planet = const_cast<double*>(Sa_mars);
 	 qu_planet = const_cast<double*>(qu_mars);
 	 Cl_planet = const_cast<double*>(Cl_mars);
 	 Sl_planet = const_cast<double*>(Sl_mars);
 	 alpha0_planet = const_cast<double*>(alpha0_mars);
 	 delta0_planet = const_cast<double*>(delta0_mars);
 	 W0_planet = const_cast<double*>(W0_mars);  
 	 radius_eq_planet = radius_eq_mars; 
 	 flattening_planet = flattening_mars;
	
 }
 
 else if (aName == "Mercury" ){
 	 PlanetSpiceCode=199;
 	 PlanetIsMagnetic=true;
 	 PlanetName=aName;
	 a_planet = const_cast<double*>(a_mercury);
         e_planet = const_cast<double*>(e_mercury);
         lm_planet = const_cast<double*>(lm_mercury);
         mu_planet = mu_mercury;
	 Omega_planet = const_cast<double*>(Omega_mercury);
         p_planet = const_cast<double*>(p_mercury); //omega+Omega
         i_planet = const_cast<double*>(i_mercury);
         pu_planet = const_cast<double*>( pu_mercury);
 	 Ca_planet = const_cast<double*>(Ca_mercury);
 	 Sa_planet = const_cast<double*>(Sa_mercury);
 	 qu_planet = const_cast<double*>(qu_mercury);
 	 Cl_planet = const_cast<double*>(Cl_mercury);
 	 Sl_planet = const_cast<double*>(Sl_mercury);
 	 alpha0_planet = const_cast<double*>(alpha0_mercury);
 	 delta0_planet = const_cast<double*>(delta0_mercury);
 	 W0_planet = const_cast<double*>(W0_mercury);  
 	 radius_eq_planet = radius_eq_mercury; 
 	 flattening_planet = flattening_mercury;
  }
  else if (aName == "Earth" ){
 	 PlanetSpiceCode=399;
         PlanetIsMagnetic =true;
 	 PlanetName=aName;
	 a_planet = const_cast<double*>(a_earth);
         e_planet = const_cast<double*>(e_earth);
         lm_planet = const_cast<double*>(lm_earth);
         mu_planet = mu_earth;
	 Omega_planet = const_cast<double*>(Omega_earth);
         p_planet = const_cast<double*>(p_earth); //omega+Omega
         i_planet = const_cast<double*>(i_earth);
         pu_planet = const_cast<double*>( pu_earth);
 	 Ca_planet = const_cast<double*>(Ca_earth);
 	 Sa_planet = const_cast<double*>(Sa_earth);
 	 qu_planet = const_cast<double*>(qu_earth);
 	 Cl_planet = const_cast<double*>(Cl_earth);
 	 Sl_planet = const_cast<double*>(Sl_earth);
 	 alpha0_planet = const_cast<double*>(alpha0_earth);
 	 delta0_planet = const_cast<double*>(delta0_earth);
 	 W0_planet = const_cast<double*>(W0_earth);  
 	 radius_eq_planet = radius_eq_earth; 
 	 flattening_planet = 1. - radius_pole_earth/radius_eq_earth;
	
  }

#ifdef USE_JUPITER	 
  else if (aName == "Jupiter" ){
 	 PlanetSpiceCode=699;
         PlanetIsMagnetic =true;
 	 PlanetName=aName;
	 a_planet = const_cast<double*>(a_jupiter);
         e_planet = const_cast<double*>(e_jupiter);
         lm_planet = const_cast<double*>(lm_jupiter);
         mu_planet = mu_jupiter;
	 Omega_planet = const_cast<double*>(Omega_jupiter);
         p_planet = const_cast<double*>(p_jupiter); //omega+Omega
         i_planet = const_cast<double*>(i_jupiter);
         pu_planet = const_cast<double*>( pu_jupiter);
 	 Ca_planet = const_cast<double*>(Ca_jupiter);
 	 Sa_planet = const_cast<double*>(Sa_jupiter);
 	 qu_planet = const_cast<double*>(qu_jupiter);
 	 Cl_planet = const_cast<double*>(Cl_jupiter);
 	 Sl_planet = const_cast<double*>(Sl_jupiter);
 	 alpha0_planet = const_cast<double*>(alpha0_jupiter);
 	 delta0_planet = const_cast<double*>(delta0_jupiter);
 	 W0_planet = const_cast<double*>(W0_jupiter);  
 	 radius_eq_planet = radius_eq_jupiter; 
 	 flattening_planet = 1. - radius_pole_jupiter/radius_eq_jupiter;
	
  }  
#endif
  else {
 	G4cout<<"No space coordinate convertor exist for the planet that "
	      <<"you have selected!"<<std::endl; 
	return;        
 }
 ListOfCoordinateSystems.clear();
 if (PlanetName =="Earth"){
 	ListOfCoordinateSystems.clear();
 	ListOfCoordinateSystems.push_back("GEO");
 	ListOfCoordinateSystems.push_back("GSM");
 	ListOfCoordinateSystems.push_back("GSE");
 	ListOfCoordinateSystems.push_back("MAG");
 	ListOfCoordinateSystems.push_back("GEODETIC");
 	ListOfCoordinateSystems.push_back("SM");
 }

 ListOfCoordinateSystems.push_back("PLA");
 ListOfCoordinateSystems.push_back("PLAG");
 ListOfCoordinateSystems.push_back("POI");
 ListOfCoordinateSystems.push_back("PEQS");
 ListOfCoordinateSystems.push_back("PSEQ");
 ListOfCoordinateSystems.push_back("PSO");
 ListOfCoordinateSystems.push_back("GEI");
 ListOfCoordinateSystems.push_back("HAE2000");
 if (PlanetIsMagnetic) {
 	ListOfCoordinateSystems.push_back("PMAG");
	ListOfCoordinateSystems.push_back("PSMAG");
	ListOfCoordinateSystems.push_back("PSM");
 }
 // G4cout<<"radius equator"<<radius_eq_planet;
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::SetDipoleAxisInPLA(G4ThreeVector axis)
{DipoleAxisInPLA=axis;
 if (PlanetIsMagnetic) ComputeMagneticTranformationMatrices(); 
 
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinatePlanet::Transform(G4ThreeVector vec_in,
						  G4String sys_in,
						  G4String sys_out)
{ G4bool sys_in_ok =false;
  G4bool sys_out_ok =false;
  for (unsigned int i=0; i<ListOfCoordinateSystems.size();i++){
  	if (sys_in  == ListOfCoordinateSystems[i]) sys_in_ok=true;
  	if (sys_out  == ListOfCoordinateSystems[i]) sys_out_ok=true;
  }
  

  if (!sys_in_ok || !sys_out_ok ) {
	G4cout << sys_in<<" and/or "<< sys_out<< " are not good systems of coordinates" <<G4endl;
	return G4ThreeVector(vec_in);
  }

  if (PlanetName =="Earth"){
         	system_in = PlanetCoordinateSystemNameFromEarthCoordinateSystem(sys_in);
		system_out = PlanetCoordinateSystemNameFromEarthCoordinateSystem(sys_out);
		
  }
  else {
		system_in = sys_in;
		system_out = sys_out;
  }
   
  ComputeSelectedMatrix();
  return Transform(vec_in);

}

////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinatePlanet::Transform(G4ThreeVector vec_in)
{
    if (system_in == system_out)
	return G4ThreeVector(vec_in);
    else
	return SelectedMatrix * vec_in;

}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinatePlanet::TransformPLAinPSM(G4ThreeVector pos_pla)
{  return PLAtoPSM * pos_pla;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector SpaceCoordinatePlanet::TransformPSMinPLA(G4ThreeVector pos_psm)
{  return PSMtoPLA * pos_psm;
}
/////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputePLAPositionFromPLAG
    (G4double altitude, G4double latitude, G4double longitude,
     G4ThreeVector & position) 
{   
    G4double EllipsoidLargeSemiaxis2=radius_eq_planet*radius_eq_planet;
    G4double EllipsoidSmallSemiaxis2=EllipsoidLargeSemiaxis2
    						*(1.-flattening_planet)
						*(1.-flattening_planet);
    G4double delta = 90. * degree - latitude;
    G4double sind = sin(delta);
    G4double cosd = cos(delta);
    G4double p = longitude;

    G4double eta = sqrt(EllipsoidLargeSemiaxis2 * sind * sind
			+ EllipsoidSmallSemiaxis2 * cosd * cosd);

    G4double xi =
	(EllipsoidLargeSemiaxis2 -
	 EllipsoidSmallSemiaxis2) * sind * cosd / eta;

    G4double r = sqrt(xi * xi + (eta + altitude) * (eta + altitude));

    G4double sbeta = xi / r;
    G4double cbeta = (eta + altitude) / r;
    G4double sth = cbeta * sind + sbeta * cosd;
    G4double cth = -sbeta * sind + cbeta * cosd;
    G4double sp = sin(p);
    G4double cp = cos(p);
    position = G4ThreeVector(sth * cp, sth * sp, cth);
    position.setMag(r);
}

//////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputePLADirectionAndPositionFromPLAG
    (G4double altitude, G4double latitude, G4double longitude,
     G4double zenith, G4double azimuth,
     G4ThreeVector & position, G4ThreeVector & direction) {
     
    
    G4double EllipsoidLargeSemiaxis2=radius_eq_planet*radius_eq_planet;
    G4double EllipsoidSmallSemiaxis2=EllipsoidLargeSemiaxis2
    						*(1.-flattening_planet)
						*(1.-flattening_planet);
    G4double delta = 90. * degree - latitude;
    G4double sind = sin(delta);
    G4double cosd = cos(delta);
    G4double p = longitude;
    G4double eta = sqrt(EllipsoidLargeSemiaxis2 * sind * sind
			+ EllipsoidSmallSemiaxis2 * cosd * cosd);
    G4double xi =
	(EllipsoidLargeSemiaxis2 -
	 EllipsoidSmallSemiaxis2) * sind * cosd / eta;
    G4double r = sqrt(xi * xi + (eta + altitude) * (eta + altitude));
    G4double sbeta = xi / r;
    G4double cbeta = (eta + altitude) / r;
    G4double sth = cbeta * sind + sbeta * cosd;
    G4double cth = -sbeta * sind + cbeta * cosd;
    G4double sp = sin(p);
    G4double cp = cos(p);
    position = G4ThreeVector(sth * cp, sth * sp, cth);
    position.setMag(r);
    G4double t = position.theta();
    G4ThreeVector VEllipsoid = G4ThreeVector(0., 0., 1.);
    VEllipsoid.setTheta(zenith);
    VEllipsoid.setPhi(180. * degree - azimuth);
    G4double VplaR, VplaT, VplaP;
    VplaR = cbeta * VEllipsoid.z() + sbeta * VEllipsoid.x();
    VplaT = -sbeta * VEllipsoid.z() + cbeta * VEllipsoid.x();
    VplaP = VEllipsoid.y();
    G4ThreeVector VPlaSphere = G4ThreeVector(VplaT, VplaP, VplaR);
    direction = -VPlaSphere.rotateY(t).rotateZ(p);
}

///////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputePLAGCoordinatesFromPLAPosition
    (const G4ThreeVector PLAposition,
     G4double & altitude, G4double & longitude, G4double & latitude)
{   G4double EllipsoidLargeSemiaxis2=radius_eq_planet*radius_eq_planet;
    G4double EllipsoidSmallSemiaxis2=EllipsoidLargeSemiaxis2
    						*(1.-flattening_planet)
						*(1.-flattening_planet);
    G4double rho = PLAposition.rho();
    G4double rho2 = rho * rho;
    G4double z = std::abs(PLAposition.z());
    G4double z2 = z * z;
    G4double r1 = std::sqrt(EllipsoidLargeSemiaxis2);
    G4double r2 = std::sqrt(EllipsoidSmallSemiaxis2);
    G4double precision = .0001;
    G4double cosa = rho / PLAposition.r();
    G4double sina = z / PLAposition.r();
    G4double alpha = std::acos(cosa);
    longitude = PLAposition.phi();

    if (r1 == r2) {
	altitude = PLAposition.r() - r1;
	if (sina < 0.)
	    latitude = -alpha;
	else
	    latitude = alpha;

    } else {
	//compute first b and t 

	G4double a = cosa * cosa / EllipsoidLargeSemiaxis2;
	a += sina * sina / EllipsoidSmallSemiaxis2;
	G4double b = cosa * rho / EllipsoidLargeSemiaxis2;
	b += sina * z / EllipsoidSmallSemiaxis2;
	b *= 2.;

	G4double c = rho2 / EllipsoidLargeSemiaxis2;
	c += z2 / EllipsoidSmallSemiaxis2;
	c -= 1.;

	G4double sqrt_rhot = std::sqrt(b * b - (4. * a * c));
	G4double t = (-b + sqrt_rhot) / (2. * a);
	G4double cosae = (rho + t * cosa) / r1;
	G4double sinae = (z + t * sina) / r2;
	G4double cosb = r2 * cosae;
	G4double sinb = r1 * sinae;
	G4double rb = std::sqrt(cosb * cosb + sinb * sinb);
	cosb /= rb;
	sinb /= rb;
	G4double beta = (std::acos(cosb));



	if (t >= 0) {		// point is below the planetary surface
	    cosa = cosb;
	    sina = sinb;
	    while (std::abs(alpha - beta) > precision) {
	    	/*std::cout<<"below"<<std::endl;
		std::cout<<alpha<<std::endl;
		std::cout<<beta<<std::endl;*/
		alpha = std::acos(cosa);
		a = cosa * cosa / EllipsoidLargeSemiaxis2;
		a += sina * sina / EllipsoidSmallSemiaxis2;
		b = cosa * rho / EllipsoidLargeSemiaxis2;
		b += sina * z / EllipsoidSmallSemiaxis2;
		b *= 2.;
		sqrt_rhot = std::sqrt(b * b - (4. * a * c));
		t = (-b + sqrt_rhot) / (2. * a);
		cosae = (rho + t * cosa) / r1;
		sinae = (z + t * sina) / r2;
		cosb = r2 * cosae;
		sinb = r1 * sinae;
		rb = std::sqrt(cosb * cosb + sinb * sinb);
		cosb /= rb;
		sinb /= rb;
		beta = (std::acos(cosb));
		cosa = cosb;
		sina = sinb;
	    }
	} else {		//point is above the planetary Surface 
	    G4ThreeVector direction1 =
		G4ThreeVector(rho, 0., z) - G4ThreeVector(r1, 0., 0.);
	    G4ThreeVector direction2 =
		G4ThreeVector(rho, 0., z) - G4ThreeVector(0., 0., r2);

	    G4double max_alpha =
		std::max(90. * degree - direction1.theta(),
			 90. * degree - direction2.theta());

	    G4double alpha1, alpha2;
	    alpha1 = alpha;


	    if (beta > max_alpha) {
		cosa = cos(max_alpha);
		sina = sin(max_alpha);
		alpha2 = max_alpha;
	    } else {
		cosa = cosb;
		sina = sinb;
		alpha2 = beta;
	    }
	    while (std::abs(alpha - beta) > precision) {
/*		std::cout<<alpha<<std::endl;
		std::cout<<beta<<std::endl;*/
		alpha = alpha2;
		a = cosa * cosa / EllipsoidLargeSemiaxis2;
		a += sina * sina / EllipsoidSmallSemiaxis2;
		b = cosa * rho / EllipsoidLargeSemiaxis2;
		b += sina * z / EllipsoidSmallSemiaxis2;
		b *= 2.;
		sqrt_rhot = std::sqrt(b * b - (4. * a * c));
		t = (-b + sqrt_rhot) / (2. * a);
		cosae = (rho + t * cosa) / r1;
		sinae = (z + t * sina) / r2;
		cosb = r2 * cosae;
		sinb = r1 * sinae;
		rb = std::sqrt(cosb * cosb + sinb * sinb);
		cosb /= rb;
		sinb /= rb;
		beta = (std::acos(cosb));
		if (alpha < beta) {
		    if (beta > max_alpha) {
			cosa = cos(max_alpha);
			sina = sin(max_alpha);
			alpha1 = alpha;
			alpha2 = alpha + max_alpha;

		    } else {
			cosa = cosb;
			sina = sinb;
			alpha1 = alpha;
			alpha2 = beta;
		    }
		} else {
		    alpha2 = (alpha1 + alpha2) / 2.;
		    cosa = cos(alpha2);
		    sina = sin(alpha2);

		}
	    }
	}

	altitude = -t;
	latitude = beta;
	if (PLAposition.z() < 0.)
	    latitude = -latitude;
	longitude = PLAposition.phi();
    }
}

///////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::SetReferenceDate(G4int year, G4int month,
						G4int day, G4int hour,
						G4int minute, G4int second)
{ReferenceDate =DateAndTime(year,month,day,hour,minute,second);
 Initialise();  
}

////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::SetReferenceDate(DateAndTime ref_date)
{
    ReferenceDate = ref_date;
    Initialise();
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::PlanetAndSun()
{ 
  //compute the time in Barycentric dynamical time used for computing
  // motion of planets around the barycenter of the solar system  
  DateAndTime TDB_date =ReferenceDate.TDB_from_UTC();

  //time form J2000 TDB in day, second and Julian century (365250 days) 
  double day= TDB_date.JulianDate()- 2451545.0;
  //double tsec= day*86400.;
  double t= day/365250.;
  double t2 =t*t;
  
  //compute the orbital parameters for planet orbit
  double a = a_planet[0]+t*a_planet[1]; // semi major axis in AU
  double e=  e_planet[0]+t*e_planet[1]+t2*e_planet[2]; //eccentricity 
  double lm=  (lm_planet[0]+(t*lm_planet[1]+t2*lm_planet[2])/3600.)*degree; //true longitude
  double omega_bar= (p_planet[0]+(t*p_planet[1]+t2*p_planet[2])/3600.)*degree; //argument of the pericenter
  double Omega= (Omega_planet[0]+(t*Omega_planet[1]+t2*Omega_planet[2])/3600.)*degree;//longitude of the node 
  omega_bar-=Omega;
  //inclination of the orbital plane with the eccliptic plane at J2000 TDB
  double inc=  (i_planet[0]+(t*i_planet[1]+t2*i_planet[2])/3600.)*degree;
 
  
  //trigonometric term for a and lm
  // see Numerical expressions for precession formulae and mean  elements 
  // for the Moon and the planets, Astronomy and Astrophysics, 282,663-683(1994)
  
  double mu=t*0.35953620;
  for (unsigned i =0;i<8;i++) {
  	G4double pu = mu*double(pu_planet[i]);
	a += 1.e-7*(double(Ca_planet[i])*std::cos(pu)+double(Sa_planet[i])*std::sin(pu));
	G4double qu= mu*double(qu_planet[i]);
	lm +=1.e-7*(double(Cl_planet[i])*std::cos(qu)+double(Sl_planet[i])*std::sin(qu));  
  }
  
  G4double pu = mu*double(pu_planet[8]);
  a += t*1.e-7*(double(Ca_planet[8])*std::cos(pu)+double(Sa_planet[8])*std::sin(pu));
  for (unsigned i =8;i<10;i++) {
	G4double qu= mu*double(qu_planet[i]);
	lm +=t*1.e-7*(double(Cl_planet[i])*std::cos(qu)+double(Sl_planet[i])*std::sin(qu));  
  }
  a*=au;
  
  
  //position of Planet in POI Planetary Orbital Inertial
  // In this system the xy plane contains the  orbit of the planet 
  // around the Sun at the period considered, with the  x axis corresponding to 
  // the periapsis of the orbit. The z axis close the orthogonal system.
  
  PlanetPosPOI = orbit_func::pos_in_orbit_system(mu_planet, a, e, lm, omega_bar, Omega);
  
  //The matrix from HAE200 to POI is the Euler matrix E(Omega,inc,omega_bar)
  //The Euler matrix convention used in CLHEP G4RotationMatrix 
  // is the same than in the apper from Fränz and Harper, Heliospheric
  // coordinate system 
  
  HAE2000toPOI =G4RotationMatrix(Omega,inc,omega_bar);
  POItoHAE2000 = HAE2000toPOI.inverse();
  
  PlanetPosHAE2000=HAE2000toPOI.inverse()*PlanetPosPOI;
  if (verbose > 0){
        G4cout<<"Sun Planet vector in HAE2000"<<std::endl;	
  	G4cout<<PlanetPosHAE2000/au<<std::endl;
  
  }
  
  //
  //Matrix from POI to PSO
  //
  G4ThreeVector XPSOinPOI =  -PlanetPosPOI/PlanetPosPOI.mag();
  XPSOinPOI.unit();
  G4ThreeVector ZPSOinPOI = G4ThreeVector(0.,0.,1.);
  G4ThreeVector YPSOinPOI =ZPSOinPOI.cross(XPSOinPOI); 
  G4RotationMatrix PSOtoPOI = G4RotationMatrix(XPSOinPOI,YPSOinPOI,ZPSOinPOI);
  
  
  //
  //Matrix from PSO to HAE2000
  //
  PSOtoHAE2000 = POItoHAE2000* PSOtoPOI;
  HAE2000toPSO = PSOtoHAE2000.inverse(); 
  if (verbose > 0){
	G4cout<<"HAE2000toPSO transformation matrix"<<std::endl;
	G4cout<<HAE2000toPSO.rowX()<<std::endl;
  	G4cout<<HAE2000toPSO.rowY()<<std::endl;
  	G4cout<<HAE2000toPSO.rowZ()<<std::endl;
  }
  
  
  
  //Orientation and rotation of Planet in GEI2000
  //---------------------------------------------
  //Data taken from Table 15.7 of Explanatory supplement to the  Astronomical
  //Almanac,1992
  //See also section 4.3 of Fränz, M and Harper, D., Heliospheric Coordinate
  //system, Corrected version of Plan. Space Science, 50, 217, 2002   
  
  G4double alpha0=alpha0_planet[0]+t*alpha0_planet[1];
  alpha0*=degree;
  G4double delta0=delta0_planet[0]+t*delta0_planet[1];
  delta0*=degree;
  G4double W0=W0_planet[0]+day*W0_planet[1];
  W0*=degree;
  G4RotationMatrix M=G4RotationMatrix(alpha0+90.*degree,90.*degree-delta0,W0);
  PLAtoGEI2000 = M.inverse();
  GEI2000toPLA = M;
  PlanetRotationAxisInGEI2000 = PLAtoGEI2000.colZ();
  if (verbose > 0){
	G4cout<<"GEI2000toPLA transformation matrix"<<std::endl;
	G4cout<<GEI2000toPLA.rowX()<<std::endl;
  	G4cout<<GEI2000toPLA.rowY()<<std::endl;
  	G4cout<<GEI2000toPLA.rowZ()<<std::endl;
  }
  
  //Orientation and rotation of Planet in GEI2000
  //---------------------------------------------
  //Data taken from Table 15.7 of Explanatory supplement to the  Astronomical
  //Almanac,1992
  //See also section 4.3 of Fränz, M and Harper, D., Heliospheric Coordinate
  //system, Corrected version of Plan. Space Science, 50, 217, 2002   
  
  alpha0=alpha0_sun[0]+t*alpha0_sun[1];
  alpha0*=degree;
  delta0=delta0_sun[0]+t*delta0_sun[1];
  delta0*=degree;
  W0=W0_sun[0]+day*W0_sun[1];
  W0*=degree;
  M=G4RotationMatrix(alpha0+90.*degree,90.*degree-delta0,W0);
  SunRotationAxisInGEI2000 = PLAtoGEI2000.colZ();
  	
 
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::TestMarinerOrbit(G4String input_file,G4String output_file)
{
 std::ifstream InputFile;
 InputFile.open(input_file, (std::ios::binary | std::ios::in));
 
 std::ofstream OutputFile;
 OutputFile.open(output_file, (std::ios::binary | std::ios::out)); 
 G4String time;
 G4ThreeVector Bpso,Bpseq,Ppso,Ppseq;
 G4double Br,Bphi,Bth,r,lon,lat;
 DateAndTime first_date;
 G4int i=0;
 
 OutputFile<<"Time [hour]"<<'\t'<<"X_PLA"<<'\t'<<'\t'<<"Y_PLA"<<'\t'<<'\t'<<"Z_PLA"
	           <<'\t'<<'\t'<<"X_PSO"<<'\t'<<'\t'<<"Y_PSO"<<'\t'<<'\t'<<"Z_PSO"<<'\t'<<'\t'
	           <<std::endl;
 while (!InputFile.eof()){
 G4double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,at;
 	
 	InputFile>>time
	         >>a1>>a2>>a3>>at
		 >>a4>>a5>>a6
		 >>Br>>Bth>>Bphi
		 >>a7>>a8>>a9
		 >>a10>>a11>>a12
		 >>r>>lat>>lon;
	lat=lat*degree;
	lon=lon*degree;	 
        Bpso=G4ThreeVector(a1,a2,a3);
	Bpseq=G4ThreeVector(a4,a5,a6);
	Ppso=G4ThreeVector(a7,a8,a9);
	Ppseq=G4ThreeVector(a10,a11,a12);
	std::strstream astream;
	astream<<G4String(time(0,4))<<'\t'
	       <<G4String(time(5,2))<<'\t'
	       <<G4String(time(8,2))<<'\t'
	       <<G4String(time(11,2))<<'\t'
	       <<G4String(time(14,2))<<'\t'
	       <<G4String(time(17,6));
	G4int year,day,mo,hour,min, isec,msec;
	G4double sec;
	
	
	astream>>year>>mo>>day>>hour>>min>>sec;
	isec=int(sec);
	msec=int(1000.*(sec-isec));
	
	DateAndTime theDate =DateAndTime(year,mo,day,hour,min,isec,msec); 
	
	SetReferenceDate(theDate);
	if (i==0) {
		first_date =theDate;
		i=1;
	}
	OutputFile.setf(std::ios::scientific);
        OutputFile.precision(4);            
	G4double thours = theDate.DifferenceInHours(first_date);
	G4ThreeVector PLApos = G4ThreeVector(0.,0.,1);
	PLApos.setRThetaPhi(r,90.*degree -lat,lon);
	G4ThreeVector PSOpos = Transform(PLApos,"PLA" ,"PSO");
	OutputFile<<thours<<'\t'<<PLApos.x()<<'\t'<<PLApos.y()<<'\t'<<PLApos.z()
	           <<'\t'<<PSOpos.x()<<'\t'<<PSOpos.y()<<'\t'<<PSOpos.z()
		   <<std::endl;
	       
	       
	
	
	 
		 
 }
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::TestImageOrbit(G4String input_file,G4String output_file)
{
 std::ifstream InputFile;
 InputFile.open(input_file, (std::ios::binary | std::ios::in));
 
 std::ofstream OutputFile;
 OutputFile.open(output_file, (std::ios::binary | std::ios::out)); 
 G4String time, day;
 G4ThreeVector Bpso,Bpseq,Ppso,Ppseq;
 DateAndTime first_date;
 G4int i=0;
 
 OutputFile<<"Time [hour]"<<'\t'<<"X_GSM"<<'\t'<<'\t'<<"Y_GSM"<<'\t'<<'\t'<<"Z_GSM"
	           <<std::endl;
 while (!InputFile.eof()){
 G4double a1,a2,a3,a4,a5,a6,a7,a8,a9;

 	InputFile>>time>>day
	         >>a1>>a2>>a3
		 >>a4>>a5>>a6
		 >>a7>>a8>>a9
		 ;

	std::strstream astream;
	astream<<G4String(time(0,2))<<'\t'
	       <<G4String(time(3,2))<<'\t'
	       <<G4String(time(6,4))<<'\t'
	       <<G4String(day(0,2))<<'\t'
	       <<G4String(day(3,2))<<'\t'
	       <<G4String(day(6,6));
	G4int year,day,mo,hour,min, isec,msec;
	G4double sec;
	
	
	astream>>day>>mo>>year>>hour>>min>>sec;
	
	
	
	
	isec=int(sec);
	msec=int(1000.*(sec-isec));
	
	DateAndTime theDate =DateAndTime(year,mo,day,hour,min,isec,msec); 
	
	//SetReferenceDate(theDate);
	PlanetManager::GetInstance()->GetMagneticField()->SetStartDate(theDate);
	if (i==0) {
		first_date =theDate;
		i=1;
	}
	OutputFile.setf(std::ios::scientific);
        OutputFile.precision(4);            
	G4double thours = theDate.DifferenceInHours(first_date);
	G4ThreeVector GCIpos = G4ThreeVector(a1,a2,a3);
	G4ThreeVector GSEpos = G4ThreeVector(a7,a8,a9);
	G4ThreeVector GSMpos = Transform(GSEpos,"GSE" ,"GSM");
	OutputFile<<thours
	           <<'\t'<<GSMpos.x()<<'\t'<<GSMpos.y()<<'\t'<<GSMpos.z()
		   <<std::endl;
	       
	       
	
	 
		 
 }
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::TestWindOrbit(G4String input_file,G4String output_file)
{
 std::ifstream InputFile;
 InputFile.open(input_file, (std::ios::binary | std::ios::in));
 
 std::ofstream OutputFile;
 OutputFile.open(output_file, (std::ios::binary | std::ios::out)); 
 G4String time, day;
 G4ThreeVector Bpso,Bpseq,Ppso,Ppseq;
 DateAndTime first_date;
 G4int i=0;
 OutputFile<<"Time [hour]"<<'\t'<<"X_GSM"<<'\t'<<'\t'<<"Y_GSM"<<'\t'<<'\t'<<"Z_GSM"
	           <<std::endl;
 
 while (!InputFile.eof()){
 G4double a1,a2,a3,a4,a5,a6,a7,a8,a9;
 	
 	InputFile>>time>>day
	         >>a1>>a2>>a3
		 >>a4>>a5>>a6
		 >>a7>>a8>>a9
		 ;
	std::strstream astream;
	astream<<G4String(time(0,2))<<'\t'
	       <<G4String(time(3,2))<<'\t'
	       <<G4String(time(6,4))<<'\t'
	       <<G4String(day(0,2))<<'\t'
	       <<G4String(day(3,2))<<'\t'
	       <<G4String(day(6,6));
	G4int year,day,mo,hour,min, isec,msec;
	G4double sec;
	
	
	astream>>day>>mo>>year>>hour>>min>>sec;
	
	
	
	isec=int(sec);
	msec=int(1000.*(sec-isec));
	
	DateAndTime theDate =DateAndTime(year,mo,day,hour,min,isec,msec); 
	
	//SetReferenceDate(theDate);
	PlanetManager::GetInstance()->GetMagneticField()->SetStartDate(theDate);
	if (i==0) {
		first_date =theDate;
		i=1;
	}
	OutputFile.setf(std::ios::scientific);
        OutputFile.precision(4);            
	G4double thours = theDate.DifferenceInHours(first_date);
	G4ThreeVector GCIpos = G4ThreeVector(a1,a2,a3);
	G4ThreeVector GSEpos = G4ThreeVector(a4,a5,a6);
	G4ThreeVector GSMpos = Transform(GSEpos,"GSE" ,"GSM");
	OutputFile<<thours
	           <<'\t'<<GSMpos.x()<<'\t'<<GSMpos.y()<<'\t'<<GSMpos.z()
		   <<std::endl;
	       
	       
	
	
	 
		 
 }
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::TestMarsOdysseyOrbit(G4String input_file,G4String output_file)
{std::ifstream InputFile;
 InputFile.open(input_file, (std::ios::binary | std::ios::in));
 
 std::ofstream OutputFile;
 OutputFile.open(output_file, (std::ios::binary | std::ios::out)); 
 G4String time;
 G4ThreeVector Bpso,Bpseq;
 G4double r,lon,lat;
 DateAndTime first_date;
 G4int i=0;
 
 OutputFile<<"Time [hour]"<<'\t'<<"X_PLA"<<'\t'<<'\t'<<"Y_PLA"<<'\t'<<'\t'<<"Z_PLA"
	           <<'\t'<<'\t'<<"X_PSO"<<'\t'<<'\t'<<"Y_PSO"<<'\t'<<'\t'<<"Z_PSO"<<'\t'<<'\t'
	           <<std::endl;
 while (!InputFile.eof()){
 G4double a4,a5,a6,a7,a8,a9,at;
 	
 	InputFile>>time
	         >>r>>lat>>lon>>at
		 >>a4>>a5>>a6
		 >>a7>>a8>>a9;
	lat=lat*degree;
	lon=lon*degree;	 
	std::strstream astream;
	astream<<G4String(time(0,4))<<'\t'
	       <<G4String(time(5,2))<<'\t'
	       <<G4String(time(8,2))<<'\t'
	       <<G4String(time(11,2))<<'\t'
	       <<G4String(time(14,2))<<'\t'
	       <<G4String(time(17,6));
	G4int year,day,mo,hour,min, isec,msec;
	G4double sec;
	
	
	astream>>year>>mo>>day>>hour>>min>>sec;
	isec=int(sec);
	msec=int(1000.*(sec-isec));
	
	DateAndTime theDate =DateAndTime(year,mo,day,hour,min,isec,msec); 
	
	SetReferenceDate(theDate);
	if (i==0) {
		first_date =theDate;
		i=1;
	}            
	G4double thours = theDate.DifferenceInHours(first_date);
	G4ThreeVector PLApos = G4ThreeVector(0.,0.,1);
	PLApos.setRThetaPhi(r,90.*degree -lat,lon);
	G4ThreeVector PSOpos = Transform(PLApos,"PLA" ,"PSO");
	OutputFile.setf(std::ios::scientific);
        OutputFile.precision(4); 
	OutputFile<<thours<<'\t'<<PLApos.x()<<'\t'<<PLApos.y()<<'\t'<<PLApos.z()
	           <<'\t'<<PSOpos.x()<<'\t'<<PSOpos.y()<<'\t'<<PSOpos.z()
	           <<std::endl;
	
	
        


	      
	
	 
		 
 }
}
////////////////////////////////////////////////////////////////////////////////
//
G4String SpaceCoordinatePlanet::PlanetCoordinateSystemNameFromEarthCoordinateSystem(G4String aName)
{ if (aName == "GEO") return "PLA";
  if (aName == "GEODETIC") return "PLAG";
  if (aName == "MAG") return "PMAG";
  if (aName == "GSM") return "PSM";
  if (aName == "GSE") return "PSO";
  if (aName == "SM") return "PSMAG";
  for (unsigned int i=0; i<ListOfCoordinateSystems.size();i++){
  	if ( aName == ListOfCoordinateSystems[i]) return aName;
  }
  
  return "";
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
//#ifdef USE_SPICE
void SpaceCoordinatePlanet::SpicePlanetAndSun()
{
 
  //Planet-Position from JPL-ephemeris
  //-----------------------------------
  SpiceInt h[13];
 
  //load frames and time file
  G4String leap_second_file = getenv("SPICE_LEAPSECOND_DATA");
  G4String frame_file = getenv("SPICE_FRAME_DATA");
  ldpool_c(leap_second_file);
  ldpool_c(frame_file);
  
  //load ephemerids
  G4String ephemerid_file =getenv("SPICE_EPHEMERID_DATA");
  spklef_c(ephemerid_file, &h[0]);
  
  
  //definition of time 
  //time definition in UTC 
  // transfert from reference date to string for cspice
    
    std::strstream astream;
    astream <<ReferenceDate.year<<'\t'
    	    <<ReferenceDate.month<<'\t'
	    <<ReferenceDate.day<<'\t'
	    <<ReferenceDate.hour<<'\t'
	    <<ReferenceDate.min<<'\t'
	    <<ReferenceDate.sec;
    G4String syear, smonth, sday, shour, smin, ssec;
    astream >> syear >> smonth >> sday >> shour >> smin >> ssec;
    G4String str1 =
	smonth + "/" + sday + "/" + syear + " " + shour + ":" + smin +
	":" + ssec;	
    const char *str2 = str1.data();
    SpiceChar *string = (SpiceChar *) str2;

    SpiceDouble et;
   
    
    // the time is converted in et = second in TDB system from J200
    //this is the time used in the ephemerides file
    str2et_c(string, &et);

    //Planet position from Sun in J2000
    SpiceDouble state[6];
    SpiceDouble lt;

    if (PlanetName =="Mars")
    	spkezr_c("mars",et,"ECLIPJ2000","none","sun",state,&lt);
    else if (PlanetName =="Mercury")
    	spkezr_c("mercury",et,"ECLIPJ2000","none","sun",state,&lt);	
    else if (PlanetName =="Earth")
    	spkezr_c("earth",et,"ECLIPJ2000","none","sun",state,&lt);
	
    
    PlanetPosHAE2000=G4ThreeVector(state[0]*km,state[1]*km,state[2]*km);
    if (verbose > 0){
        G4cout<<"Sun Planet vector in HAE2000 computed with Spice"<<std::endl;	
  	G4cout<<PlanetPosHAE2000/au<<std::endl;
  
    }
   
   //PSO system
   //-------------
    G4ThreeVector XPSOinHAE2000 = -PlanetPosHAE2000;
    XPSOinHAE2000.unit();
    G4ThreeVector ZPSOinHAE2000 =PlanetPosHAE2000.cross(G4ThreeVector(state[3],
    								    state[4],
								    state[5]));
   
    G4ThreeVector YPSOinHAE2000 = ZPSOinHAE2000.cross(XPSOinHAE2000);
    PSOtoHAE2000 = G4RotationMatrix(XPSOinHAE2000, YPSOinHAE2000, ZPSOinHAE2000);
    HAE2000toPSO = PSOtoHAE2000.inverse();
    if (verbose > 0){
	G4cout<<"HAE2000toPSO transformation matrix"<<std::endl;
	G4cout<<HAE2000toPSO.rowX()<<std::endl;
  	G4cout<<HAE2000toPSO.rowY()<<std::endl;
  	G4cout<<HAE2000toPSO.rowZ()<<std::endl;
    }
   //Orientation and rotation of Planet in GEI2000
   //---------------------------------------------
    
    SpiceDouble rot[3][3];
    tipbod_c("J2000", PlanetSpiceCode, et, rot);
    
    GEI2000toPLA =
        G4RotationMatrix(G4ThreeVector(rot[0][0], rot[1][0],rot[2][0]),
			 G4ThreeVector(rot[0][1], rot[1][1],rot[2][1]),
			 G4ThreeVector(rot[0][2], rot[1][2],rot[2][2]));
    PLAtoGEI2000 = GEI2000toPLA.inverse();
    PlanetRotationAxisInGEI2000 = PLAtoGEI2000.colZ();
    if (verbose > 0){
	G4cout<<"GEI2000toPLA transformation matrix from Spice"<<std::endl;
	G4cout<<GEI2000toPLA.rowX()<<std::endl;
  	G4cout<<GEI2000toPLA.rowY()<<std::endl;
  	G4cout<<GEI2000toPLA.rowZ()<<std::endl;
    }
    
    //Orientation and rotation of Sun in GEI2000
    //-------------------------------------------
    tipbod_c("J2000", 10, et, rot);
    SunRotationAxisInGEI2000 = G4ThreeVector(rot[2][0], rot[2][1],rot[2][2]);
    
    //PSO system
    
}
//#endif


////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputeTranformationMatrices()
{

  // Caution!!                                                                 
  // If a rotation matrix M12 rotates  the cartesian axis of system 1 such that 
  // that after the rotation these axes coincide with the axis of system 2 then
  // for any vector V we have  
  // V1 =M12  * V2 and  V2= M12.inverse() * V1 
  // where V1 and V2 are 
  // G4ThreeVector objects  defining the components of V in the system 1 and 2,
  // respectively. 

  //HAE2000toPLA PLAtoHAE2000                                                 
   
  HAE2000toPLA = GEI2000toPLA*HAE2000toGEI2000;
  PLAtoHAE2000 = HAE2000toPLA.inverse();
     

  // HAE2000toPES and HAE2000PES
  // We define PEQS as the Planetary Equatorial Solar system
  // In this system the z axis represents the rotation axis of the planet
  // The x-z plane contains the planet-sun axis, the x axis oriented in direction
  // of the sun. 
  // The y axis close the system Y=Z x X

  G4ThreeVector ZPEQSinHAE2000 =  GEI2000toHAE2000* PlanetRotationAxisInGEI2000;
  G4ThreeVector YPEQSinHAE2000 =  ZPEQSinHAE2000.cross(-PlanetPosHAE2000);
  YPEQSinHAE2000.unit();
  G4ThreeVector XPEQSinHAE2000 =  YPEQSinHAE2000.cross(ZPEQSinHAE2000);
  PEQStoHAE2000 = G4RotationMatrix(XPEQSinHAE2000, YPEQSinHAE2000, ZPEQSinHAE2000);
  HAE2000toPEQS = PEQStoHAE2000.inverse();

  // HAE2000toPSEQ and HAE2000PSEQ
  // We define PSEQ as the Planetary Solar Equatorial system
  // In this system the x axis represents the planet-sun line
  // The y axis is parallel too the sun equatorial plane
  // The Z axis close the system Z=X x Y  
  
  G4ThreeVector SunRotationAxisinHAE2000 = GEI2000toHAE2000* SunRotationAxisInGEI2000;
  G4ThreeVector XPSEQinHAE2000  = -PlanetPosHAE2000/PlanetPosHAE2000.mag();
  G4ThreeVector YPSEQinHAE2000  =  SunRotationAxisinHAE2000.cross(XPSEQinHAE2000);
  YPSEQinHAE2000 = YPSEQinHAE2000/YPSEQinHAE2000.mag();
  G4ThreeVector ZPSEQinHAE2000  =  XPSEQinHAE2000.cross(YPSEQinHAE2000);   
  PSEQtoHAE2000 = G4RotationMatrix(XPSEQinHAE2000, YPSEQinHAE2000, ZPSEQinHAE2000);
  HAE2000toPSEQ = PSEQtoHAE2000.inverse();
  
  //For magnetic planet only
  
  if (PlanetIsMagnetic) ComputeMagneticTranformationMatrices();
  	
	
  
     
    
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputeMagneticTranformationMatrices()
{// HAE2000toPMAG, PMAGtoHAE2000
  G4ThreeVector ZPMAGInHAE2000 = PLAtoHAE2000*DipoleAxisInPLA;
  if (verbose > 0){
  	G4cout<<"Dipole Axis "<<DipoleAxisInPLA<<std::endl;
  }
  G4ThreeVector  ZPEQSinHAE2000 = PEQStoHAE2000.colZ();
  
  if (PLAtoHAE2000.colZ().dot(ZPMAGInHAE2000)>(1.-1e-9)){
	HAE2000toPMAG=HAE2000toPLA;
	PMAGtoHAE2000=PLAtoHAE2000;
	
  }
  else if  (ZPEQSinHAE2000.dot(ZPMAGInHAE2000)<(-1.+1e-9)){
	PMAGtoHAE2000=G4RotationMatrix(-PLAtoHAE2000.colX(),
					        PLAtoHAE2000.colY(),
						ZPMAGInHAE2000);
	HAE2000toPMAG = PMAGtoHAE2000.inverse(); 				
	
  }
  else { 
	G4ThreeVector YPMAGInPLA = G4ThreeVector(-DipoleAxisInPLA.y(),
	                                               DipoleAxisInPLA.x(),
						       0.);
	YPMAGInPLA =YPMAGInPLA/YPMAGInPLA.mag();				       
	G4ThreeVector YPMAGInHAE2000 = PLAtoHAE2000*YPMAGInPLA;
	G4ThreeVector XPMAGInHAE2000 = YPMAGInHAE2000.cross(ZPMAGInHAE2000);
	PMAGtoHAE2000=G4RotationMatrix(XPMAGInHAE2000,
				       YPMAGInHAE2000,
				       ZPMAGInHAE2000);
	HAE2000toPMAG = PMAGtoHAE2000.inverse();
  }
	
// HAE2000toPSM, PSMtoHAE2000
  
  G4ThreeVector XPSMInHAE2000  = -PlanetPosHAE2000;
  XPSMInHAE2000  = XPSMInHAE2000/XPSMInHAE2000.mag();
  G4ThreeVector YPSMInHAE2000 = ZPMAGInHAE2000.cross(XPSMInHAE2000);
  YPSMInHAE2000=YPSMInHAE2000/YPSMInHAE2000.mag();
  G4ThreeVector ZPSMInHAE2000 = XPSMInHAE2000.cross(YPSMInHAE2000);
  PSMtoHAE2000 = G4RotationMatrix(XPSMInHAE2000, 
				    YPSMInHAE2000, 
				    ZPSMInHAE2000);
  HAE2000toPSM=PSMtoHAE2000.inverse();
  PLAtoPSM=HAE2000toPSM*PLAtoHAE2000;
  PSMtoPLA=PLAtoPSM.inverse();	
  G4ThreeVector ZPSM_in_PLA = PSMtoPLA.colZ();
  G4ThreeVector XPSM_in_PLA = PSMtoPLA.colX();
  TiltAngle= DipoleAxisInPLA.angle(ZPSM_in_PLA);
  if (verbose > 0){
   G4cout<<"TiltAngle "<<TiltAngle/degree<<std::endl;
  }
  G4double angle1 = DipoleAxisInPLA.angle(XPSM_in_PLA);
  if ((angle1/degree) > 90.)  TiltAngle= -TiltAngle; 	
  
  
  
  // HAE2000toPSMAG, PSMAGtoHAE2000
  G4ThreeVector XPSMAGInHAE2000 = YPSMInHAE2000.cross(ZPMAGInHAE2000);
  PSMAGtoHAE2000 = G4RotationMatrix(XPSMAGInHAE2000, 
				      YPSMInHAE2000, 
				      ZPMAGInHAE2000);
  HAE2000toPSMAG=PSMAGtoHAE2000.inverse();
}
////////////////////////////////////////////////////////////////////////////////
//
void SpaceCoordinatePlanet::ComputeSelectedMatrix()
{

    G4RotationMatrix HAE2000_to_out = G4RotationMatrix();
    G4RotationMatrix in_to_HAE2000 = G4RotationMatrix();
    SelectedMatrix = G4RotationMatrix();
    

    if (system_in == system_out)
	return;

    if (system_in == "PLA")
	in_to_HAE2000= PLAtoHAE2000;
    else if (system_in == "PEQS")
	in_to_HAE2000= PEQStoHAE2000;
    else if (system_in == "PSEQ")
	in_to_HAE2000= PSEQtoHAE2000;
    else if (system_in == "PSO")
	in_to_HAE2000= PSOtoHAE2000;
    else if (system_in == "PMAG")
	in_to_HAE2000= PMAGtoHAE2000;
    else if (system_in == "PSMAG")
	in_to_HAE2000= PSMAGtoHAE2000;
    else if (system_in == "PSM")
	in_to_HAE2000= PSMtoHAE2000;	
		

    if (system_out == "PLA")
	 HAE2000_to_out =  HAE2000toPLA;
    else if (system_out == "PEQS")
	 HAE2000_to_out =  HAE2000toPEQS;
    else if (system_out == "PSEQ")
	 HAE2000_to_out =  HAE2000toPSEQ;
    else if (system_out == "PSO")	 
	 HAE2000_to_out =  HAE2000toPSO;
    else if (system_out == "PMAG")	 
	 HAE2000_to_out =  HAE2000toPMAG;
    else if (system_out == "PSMAG")	 
	 HAE2000_to_out =  HAE2000toPSMAG;
    else if (system_out == "PSM")	 
	 HAE2000_to_out =  HAE2000toPSM;  	 

    SelectedMatrix =  HAE2000_to_out* in_to_HAE2000;

}

///////////////////////////////////////////////////////////////////////////////

void SpaceCoordinatePlanet::Initialise()
{
//#ifdef USE_SPICE	 
    if (UseSpice) {
    	SpicePlanetAndSun();
    }
    else {
    	PlanetAndSun();
    }	  	 	 
//#endif
//#ifndef USE_SPICE
//    PlanetAndSun();
//#endif
    ComputeTranformationMatrices();
    ComputeSelectedMatrix();
}

///////////////////////////////////////////////////////////////////////////////
