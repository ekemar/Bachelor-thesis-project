#ifndef PLANETOCOSPrimaryGeneratorAction_h
#define PLANETOCOSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "MYGeneralParticleSource.hh"


#ifdef USE_ROOT_SOURCE
class ROOTPartSourceGenerator;	
#endif

class G4ParticleGun;
class PLANETOCOSPrimaryMessenger;
class G4Event;
class G4PrimaryVertex;
class G4ParticleDefinition;


class PLANETOCOSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction 
{
  public:
    	PLANETOCOSPrimaryGeneratorAction();    
    	~PLANETOCOSPrimaryGeneratorAction();

  public:
    	void GeneratePrimaries(G4Event* anEvent);
    	void (PLANETOCOSPrimaryGeneratorAction::*pGeneratePrimaries)(G4Event* anEvent);
    	void GenerateStandardPrimaries(G4Event* anEvent);
	
    	void GeneratePrimariesForComputingCutOffRigidity(G4Event* anEvent);
    	void SetRigidity(G4double aRigidity);
    
// only for spherical geometry    
    	G4bool SetPositionAndDirection(const G4String CoordSys, 
                                       const G4double anAlt,
				       const G4double aLong,
				       const G4double aLat,
				       const G4double aZenith,
				       const G4double  anAzimuth); 
    
    	void SetPositionAndDirection(	const G4String CoordSys, 
                                 	const G4ThreeVector position,
				 	const G4ThreeVector direction); 				  
    
    	G4bool SetPosition(const G4String CoordSys, 
                     	   const G4double anAlt,
			   const G4double aLong,
                     	   const G4double aLat);
    
    	G4bool SetPosition(const G4String CoordSys, 
                     	   const G4ThreeVector aPosition);
			   
			   
			   
    
    	G4bool SetPositionOnDipoleMagneticShell(const G4String Reference, 
						const G4double L, 
		      				const G4double latitude, 
						const G4double longitude);  		      		     				 
    
    	G4bool SetDirection(const G4String CoordSys, 
                      	    const G4double aZenith,
			    const G4double anAzimuth);
    
    	G4bool SetDirection(const G4String CoordSys, 
                      	    const G4ThreeVector aDirection); 	
				 
// only for flat geometry
	
	void SetPosition( const G4ThreeVector aPosition);
	void SetDirection( const G4ThreeVector aDirection);
	void SetDirection( const G4double aZenith, const G4double anAzimuth);
	
// source altitude
	G4double GetSourceAltitude();	
	      
// For magnetic analysis  
    
        void AddValuesToRigidityVector(G4int nvalues,G4double val1, G4double step);
        void SetDefaultRigidityVector();
    
        void SelectTypeOfPrimaries(G4String aString);
        void PrintPrimaryVertexInformation(const G4PrimaryVertex* aPrimary);
        void PrintBfieldAtPrimary();
    
        inline void SetDetectPrimary(G4bool aBool) {DetectPrimary =aBool;}	
        inline G4double GetNbOfPrimariesForScaling(){return G4double(nb_of_primaries_for_scaling_flux);}
        inline void SetNbOfPrimariesForScaling(G4int n){nb_of_primaries_for_scaling_flux=n;}
        inline void ResetNbOfPrimariesForScaling()
                                  {nb_of_primaries_for_scaling_flux=0;}	      
 
//User spectrum
    
        void ReadUserSpectrum(G4String file_name);
	
	void RandomIsotropicDistribution(
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
	const G4int eventNumber);	
    
//modulated galactic cosmic ray flux
        void SelectModulatedGalacticFlux(G4ParticleDefinition*, 
                                     G4double phi_mod,
				     G4double Enuc_min,
				     G4double Enuc_max);
    	void SelectModulatedGalacticFluxAtSolMin(G4ParticleDefinition*,
                                             G4double Enuc_min,
				             G4double Enuc_max);
    	void SelectModulatedGalacticFluxAtSolMax(G4ParticleDefinition*,
                                             G4double Enuc_min,
				             G4double Enuc_max);
    	void SelectMeanModulatedGalacticFlux(G4ParticleDefinition*,
                                         G4double Enuc_min,
				         G4double Enuc_max);
					 
	inline void ResetRigidityIndex() {rigidity_index=0;}
     	inline G4int GetNumberOfRigidity() const {return rigidity_values.size();}
     	inline void ResetRigidityVector(){rigidity_values.clear();}
     	inline G4ThreeVector GetPLAPosition(){return PLAPosition;}
      	inline G4ThreeVector GetPLADirection(){return PLADirection;}
       	inline void GetPLAGPosition(G4double& PLAGalt,G4double& PLAGlat,G4double& PLAGlong)
                       				{PLAGalt=PLAGaltitude;
		        			 PLAGlat=PLAGlatitude;
						 PLAGlong=PLAGlongitude;}
						 
    	inline MYGeneralParticleSource* GetParticleSource()
                                  {return myParticleSource;}
    	inline void SetVerbosity(G4int n){verbosity=n;} 
	
  // Cutoff definition
  	void ReadCutOffVsDirection(G4String file_name);					 
        void SetCutoffRigidity(G4double aVal); 
        void SetConsiderCutoff(G4bool aBool); 
        inline void SetCutoffType(G4String aType){cutoff_type= aType;}


#ifdef USE_ROOT_SOURCE
	inline void SetUseRootSource(bool aBool){UseRootSource =aBool;}
	void SetVar1Type(G4String aString);
	void SetVar2Type(G4String aString);
	void ReadROOTFirstSourceHisto(G4String file_name,G4String path);
	void ReadROOTSecondSourceHisto(G4String file_name,G4String path);
#endif 
				    
				    
  private:
    	MYGeneralParticleSource* myParticleSource;
    	PLANETOCOSPrimaryMessenger* myMessenger;
    
    
    	G4ThreeVector PLAPosition;
    	G4ThreeVector PLADirection;
    	G4double PLAGaltitude, PLAGlongitude, PLAGlatitude;
    	G4double Rigidity;
    
    	//member data for defining the binning shema in rigidity for computing
    	// rigidity filter
    	//---------------------------------------------------------------
    	std::vector<G4double> rigidity_values;
    	G4int rigidity_index;
    	G4int verbosity;
    
    	G4bool  DetectPrimary;
    	G4int nb_of_primaries_for_scaling_flux;
    	G4double scaling_flux;
    
 	// type of geometry
    	//----------------------
    	G4bool TypeOfGeometry;
	
	
	//Cutoff rigidity
	//-----------------
	
	G4bool check_cutoff;
    	G4String  cutoff_type; //DIRINTERPOLATE, FIXED
    	G4double  cutoff_rigidity;                //fixed_cut_off
    	G4double max_cutoff,min_cutoff;
    	std::vector< double > cutoff_phi;
    	std::vector< double > cutoff_theta;
    	std::vector< std::vector< double > > cutoff_rigidities;
    	G4bool consider_cutoff;
    	G4bool cutoff_hasbeen_already_defined;
    
private:    
    	void ComputePrimaries(G4Event* anEvent, G4double& energy,
                     G4double& rigidity, G4ThreeVector& position,
		     G4double& zenith, G4double& azimuth );
    	
	G4double ComputeCutoff(G4double zenith_angle, 
                           G4double azimuth_angle, 
			   G4ThreeVector position);
	
	G4double  ModulatedDifferentialGalacticFlux(G4ParticleDefinition*,
                                                G4double energy,
						G4double phi_mod);
#ifdef USE_ROOT_SOURCE
        void ComputePrimariesFromROOTHisto(G4Event* anEvent, G4double& energy,
                                           G4double& rigidity, G4ThreeVector& position,
			                   G4double& zenith, G4double& azimuth );
	bool UseRootSource;
	ROOTPartSourceGenerator* theROOTSourceGenerator;
	G4String var1_type;
	G4String var2_type;
	
	
	
#endif 
};

#endif


