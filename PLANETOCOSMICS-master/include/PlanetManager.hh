#ifndef PlanetManager_h 
#define PlanetManager_h 1

#include "globals.hh"

class PlanetAtmosphereModel;
class PlanetMagneticField;
class PlanetSoil;
class PlanetMessenger;
class PlanetManager 
{
public:	
	~PlanetManager();
	static PlanetManager* GetInstance();
	void SelectPlanet(G4String aPlanetName);
	inline G4String GetPlanetName(){return PlanetName;}
	inline PlanetMagneticField* GetMagneticField(){return theMagneticField;}
	inline PlanetAtmosphereModel* GetAtmosphereModel(){return theAtmosphereModel;}
	inline PlanetSoil* GetSoil(){return theSoil;}
	inline G4double GetRplanet(){return Rplanet;}
	inline G4double GetGplanet(){return gplanet;}
	inline G4double GetRplanetAtEquator(){return rplanet_eq;}
	inline G4double GetRplanetAtPole(){return rplanet_pole;}
	
	
private:
 	PlanetMessenger* theMessenger;
	PlanetMagneticField* theMagneticField;
	PlanetAtmosphereModel* theAtmosphereModel;
	PlanetSoil* theSoil; 
	G4String PlanetName;
	G4double rplanet,Rplanet,gplanet;
	G4double rplanet_eq,rplanet_pole;
private: 
        static PlanetManager* instance;
	

private:
	PlanetManager(G4String aName="");
	void BuildUnitTable();	
		 
};
#endif
