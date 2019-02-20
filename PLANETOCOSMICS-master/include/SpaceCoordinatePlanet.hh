#ifndef SPACECOORDINATEPLANET1_HH
#define SPACECOORDINATEPLANET1_HH 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4RotationMatrix.hh"
#include "DateAndTime.hh"


class SpaceCoordinatePlanet 
{
private: 
        static SpaceCoordinatePlanet* instance; 

private:
        SpaceCoordinatePlanet(); 

public:
	~SpaceCoordinatePlanet() ;
	
	static  SpaceCoordinatePlanet* GetInstance(); 
	
	//Set methods
	void SetReferenceDate(G4int year,G4int month ,G4int day,
	                       G4int hour,G4int minute,G4int second);
	void SetReferenceDate(DateAndTime ref_date);  		       
	void SetSystemInAndOut(G4String sys_in, G4String sys_out);
	void SetPlanetName(G4String aName);
	void SetDipoleAxisInPLA(G4ThreeVector axis);
	inline void SetVerbosity(G4int val){verbose=val;};
	
	
//#ifdef USE_SPICE	 
	inline void SetUseSpice(G4bool aBool){UseSpice= aBool;}; 	 
//#endif	         
	 
	
	//Get Methods
	inline DateAndTime GetReferenceDate() {return ReferenceDate;} 
	inline G4ThreeVector GetPlanetPositionInHAE2000(){return PlanetPosHAE2000;}
	inline G4ThreeVector GetPlanetRotationAxisInGEI2000(){return PlanetRotationAxisInGEI2000;}
	inline G4double GetTiltAngle(){return TiltAngle;}
	inline G4bool GetPlanetIsMagnetic(){return PlanetIsMagnetic;}
	
	
	inline std::vector< G4String > GetListOfCoordinateSystems()
					{return ListOfCoordinateSystems;}
	
	G4ThreeVector Transform(G4ThreeVector vec_in, G4String sys_in,
			        G4String sys_out);
	G4ThreeVector Transform(G4ThreeVector vec_in);
	G4ThreeVector TransformPLAinPSM(G4ThreeVector pos_pla);
	G4ThreeVector TransformPSMinPLA(G4ThreeVector pos_psm);				   
	
        void ComputePLAPositionFromPLAG
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4ThreeVector& position);
	 
	void ComputePLADirectionAndPositionFromPLAG
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4double zenith, G4double azimuth,
		     G4ThreeVector& position, G4ThreeVector& direction);				    	 
         
	void ComputePLAGCoordinatesFromPLAPosition
	             (const G4ThreeVector PLAposition, 
		      G4double& altitude, G4double& longitude, G4double& latitude);  	 	
	
	void TestMarinerOrbit(G4String input_file,G4String output_file);
	void TestMarsOdysseyOrbit(G4String input_file,G4String output_file);
	void TestImageOrbit(G4String input_file,G4String output_file);
	void TestWindOrbit(G4String input_file,G4String output_file);
	double GetRadiusAtEquator(){return radius_eq_planet;}
	double GetRadiusAtPole(){return (1.-flattening_planet)*radius_eq_planet;}
	 
private:
        
	
	G4String PlanetCoordinateSystemNameFromEarthCoordinateSystem(G4String aName);
	
	
	//Compute Planet position and orientation and Sun orientation 
	void SpicePlanetAndSun();
	void PlanetAndSun();
	
	
	
	//Compute the transformation coordinate matrices 
        void ComputeTranformationMatrices();
	void ComputeMagneticTranformationMatrices();
	
	//Compute the selected matrix of transformation according to 
	// system_in and system_out 
        void ComputeSelectedMatrix();
	
   	//Recompute everything after change of the reference date
	void Initialise();
	
	


//Attributes

private:
	//ReferenceDate
	DateAndTime ReferenceDate; // in UTC
	
	//PlanetName
	G4String PlanetName;

	G4int PlanetSpiceCode;
	
	//List of coordinate system
	std::vector<G4String> ListOfCoordinateSystems;
	
	
	//Magnetic 
	
	G4bool PlanetIsMagnetic;
	G4ThreeVector DipoleAxisInPLA; 
	
	
	//Coordinate system list
	//----------------------      
 
	//Rotation matrix
	//----------
	
	G4RotationMatrix HAEtoPLA,HAEtoPSO;
	
	//HAE2000 Heliocentric Aries Ecliptic  Sytem at Julian Day 2000. TDB 
	//GEI2000 Geocentric Equatorial Inertial at Julian Day 2000. TDB
	G4RotationMatrix HAE2000toGEI2000,GEI2000toHAE2000, SelectedMatrix;
	
	//PLA is the planetocentric system that is fixed to the planet 
	G4RotationMatrix GEI2000toPLA,PLAtoGEI2000,HAE2000toPLA, PLAtoHAE2000;
	
	//PEQS is Planetocentric Equatorial Solar System
	//XY plane parallel to the equator
	//XZ plane contains the PlanetSun line 
	//Z axis rotation axis of the planet
	G4RotationMatrix HAE2000toPEQS, PEQStoHAE2000;
	
	//POI is Planetary Orbital Inertial system
	//XY plane contains the  orbit of the planet 
        //X axis sun-perapsis line. 
	//Z axis close the orthogonal system.
	G4RotationMatrix HAE2000toPOI, POItoHAE2000;
	
	//PSEQ is Planetocentric Equatorial Solar System
	// X axis Sun-Planet line
	// XZ plane contains rotation axis of the sun z toward north pole
	// Y axis in the sun equator  
	G4RotationMatrix HAE2000toPSEQ, PSEQtoHAE2000;
	
	//PSO is Planetocentric Solar orbital system
	// X axis Planet-Sun direction
	// Y axis  as X lies in the orbit plane
	// Z close the system
	G4RotationMatrix HAE2000toPSO, PSOtoHAE2000;
	
	//PMAG is Planetocentric Magnetic coordinate
	// Z axis is the DipoleAxis
	// XZ plane contains the planet rotation axis
	// Y close the system
	G4RotationMatrix HAE2000toPMAG, PMAGtoHAE2000;
	
	//PSM is the Planetocentric Solar Magnetospheric system
	// X axis is the planet sun direction
	//  XZ plane conatins the magnetic dipole axis of the planet
	// Z point toward the magnetic dipole axis
	//Y closes the system
	G4RotationMatrix HAE2000toPSM, PSMtoHAE2000;
	G4double TiltAngle;
	
	
	//PSMAG is the Planetocentric Solar Magnetic system
	// Z axis is the magnetic dipole axis
	// XZ plane same than for PSM
	//Y axis same than for PSM
	G4RotationMatrix HAE2000toPSMAG, PSMAGtoHAE2000;
	
	//PSM is the Planetocentric Solar Magnetospheric system
	// X axis is the planet sun direction
	//  XZ plane conatins the magnetic dipole axis of the planet
	// Z point toward the magnetic dipole axis
	//Y closes the system
	G4RotationMatrix PLAtoPSM, PSMtoPLA;
	
	
	      	
	G4double eps0_2000;
	
	
	// Sun-Planet vector in HAE2000
	 G4ThreeVector PlanetPosHAE2000;
	 //Sun Planet Vector in POI 
	 G4ThreeVector PlanetPosPOI; 
	 
	 //Sun rotation axis in GEI
	 G4ThreeVector SunRotationAxisInGEI2000;
	  
//#ifdef USE_SPICE	 
	 G4bool UseSpice;	 
//#endif	 
	//Planet rotattion axis
	G4ThreeVector PlanetRotationAxisInGEI2000;
	
	//system in and out 
	G4String system_in;
	G4String system_out;
	
	G4int verbose;
	
	//Orbital element
	double *a_planet, *e_planet, *lm_planet, mu_planet, *Omega_planet;
        double *p_planet, *i_planet, *pu_planet, *Ca_planet, *Sa_planet;
	double *qu_planet, *Cl_planet, *Sl_planet, *alpha0_planet;
	double *delta0_planet, *W0_planet, radius_eq_planet, flattening_planet;

	
/*	double e_planet[3];
        double lm_planet[3];
        double mu_planet;
	double Omega_planet[3];
        double p_planet[3]; //omega+Omega
        double i_planet[3];
        double pu_planet[9];
 	double Ca_planet[9];
 	double Sa_planet[9];
 	double qu_planet[9];
 	double Cl_planet[9];
 	double Sl_planet[9];
//orientation data
 	double alpha0_planet[2];
 	double delta0_planet[2];
 	double W0_planet[2]; //flattening and radius of planet 
 	double radius_eq_planet; 
 	double flattening_planet ; */
	         
} ;

#endif
