#ifndef PlanetAtmosphereModel_h
#define PlanetAtmosphereModel_h 1

#include "globals.hh"
#include"G4strstreambuf.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "PlanetUnits.hh"


class PLANETOCOSGeometryConstruction;
class PlanetAtmosphereModel 
{
  public:
    PlanetAtmosphereModel(G4String aName);
    ~PlanetAtmosphereModel();

  public:      
     
    void ReadAtmosphereComposition(G4String);
    void WriteAtmosphereComposition(G4String);
    void ComputeAtmosphericLayers(G4double alt_min,G4double alt_max,
                                  G4double perc_depth, G4double min_tickness,
				  G4double max_tickness,
				  const std::vector<double> flux_detector_altitudes, 
                                  const std::vector<double> flux_detector_depth,
				  std::vector< double >& Altitudes,
				  std::vector< double >& Depth,
				  std::vector< int >& IndexDetectionBoundary,
				  std::vector< G4Material* >& Material);
   
   void ComputeAtmosphericLayers1(G4double alt_min,G4double alt_max,
                                           G4double perc_depth, G4double min_tickness,
				           G4double max_tickness,
				           const std::vector<double> flux_detector_altitudes, 
                                           const std::vector<double> flux_detector_depth,
				           std::vector< double >& Altitudes,
				           std::vector< double >& Depth,
				           std::vector< int >& IndexDetectionLayers,
					   bool detect_at_top,
					   bool detect_at_ground,
				           std::vector< G4Material* >& Materials) ; 
    
    //Set methods	
    
    void SetAtmosphericModel(G4String aName);
    inline void SetGeometry(PLANETOCOSGeometryConstruction* aGeometry)
    						{theGeometry = aGeometry;}
    virtual void SetDefault()=0;
    //Get method
    inline std::vector< G4String > GetListOfAtmosphericModels(){return ListOfAtmosphericModels;} 
  protected:
     
    //  Variable for atmospheric definition by table   
    //----------------------------------------------
    
    std::vector< std::vector<double> >    atmosphere_mass_density_elements;
    std::vector< G4Element* >  elements_in_atmosphere;
    std::vector< double >  atmosphere_altitude;
    std::vector< double >  atmosphere_mass_density;
    std::vector< double >  atmosphere_number_density; 
    std::vector< double >  atmosphere_temperature;
    std::vector< double >  atmosphere_pressure;
    std::vector< double>   atmosphere_depth;
     
    G4double altitude_unit;
    G4double mass_density_unit;
    G4double number_density_unit;
    G4double temperature_unit;
    G4double pressure_unit;
    G4bool number_of_particle_composition;
    G4bool O2_given;
    G4bool N2_given;
     
    //Type of atmospheric model
    //------------------------
    G4String AtmosphericModel; //TABLE
    std::vector<G4String> ListOfAtmosphericModels;
    
    //Planet name
    //-----------
    G4String PlanetName;
    
    //Pointer on the geometry
    //---------------------
    PLANETOCOSGeometryConstruction* theGeometry;
    
   
     
  protected:
    
    void ComputeAtmosphericLayersFromTable(G4double alt_min,G4double alt_max,
                                           G4double perc_depth, G4double min_tickness,
				           G4double max_tickness,
				           const std::vector<double> flux_detector_altitudes, 
                                           const std::vector<double> flux_detector_depth,
				           std::vector< double >& Altitudes,
				           std::vector< double >& Depth,
				           std::vector< int >& IndexDetectionBoundary,
				           std::vector< G4Material* >& Materials) ;
    //In this method detection lyear are defined for detectioin os 2nd particle in atmosphere
    //-------------------------------------------------------------------------------------------
    
    void ComputeAtmosphericLayersFromTable1(G4double alt_min,G4double alt_max,
                                           G4double perc_depth, G4double min_tickness,
				           G4double max_tickness,
				           const std::vector<double> flux_detector_altitudes, 
                                           const std::vector<double> flux_detector_depth,
				           std::vector< double >& Altitudes,
				           std::vector< double >& Depth,
				           std::vector< int >& IndexDetectionLayers,
					   bool detect_at_top,
					   bool detect_at_ground,
				           std::vector< G4Material* >& Materials) ;
    G4bool WithMaterial;					   
    virtual void ComputeAtmosphereTableFromModel(G4double alt_min,G4double alt_max)=0;				   
    
    
    std::vector <G4Element* > 
                 ComputeCompositionFromChemicalFormula(G4String  aFormula);				    
      
     
    G4double InterpolateAltitudeAtDepth(const std::vector< double> altitudes,
                                        const std::vector< double> depths,
					const std::vector< double> densities,
					G4double aDepth);
    
    G4double InterpolateDepthAtAltitude(G4double anAltitude);
    
    
    	
      				    
};

#endif

