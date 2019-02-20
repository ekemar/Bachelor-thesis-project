#ifndef PLANETOCOSGeometryConstruction_h
#define PLANETOCOSGeometryConstruction_h 1

#include "globals.hh"
#include"G4strstreambuf.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "EarthAtmosphereModel.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class PLANETOCOSGeometryMessenger;
class PLANETOCOSSD;
class PLANETOCOSSoilSD;
class G4UserLimits;
class PlanetMagneticField;
class EarthAtmosphereModel;
class PlanetAtmosphereModel;
class G4Material;

class PLANETOCOSGeometryConstruction : public G4VUserDetectorConstruction
{
  public:
    PLANETOCOSGeometryConstruction();
    ~PLANETOCOSGeometryConstruction();

  public:
  
    G4VPhysicalVolume* Construct();
    
    inline void ReadAtmosphereCompositionTable(G4String file_name)
                     {theAtmosphericModel->ReadAtmosphereComposition(file_name);}

    inline G4LogicalVolume* GetMagnetosphereLogicalVolume(){return logicMagnetosphere;}
    inline G4double GetZpositionForAtmosphereBottom() const
                                        {return ZpositionForAtmosphereBottom;}
    inline G4double GetZpositionAtZeroAltitude()const
                                        {return ZpositionAtZeroAltitude;} 
    
    void SetAtmosphereMaxStepLength(G4double);
    void SetMagnetosphereMaxStepLength(G4double);
    void SetMagnetosphereMaxTrackLength(G4double);
    void SetMagnetosphereMaxTrackDuration(G4double);
    
    void DefineDefaultGeometry();
    void UpdateGeometry();
    
    inline void RemoveAllDetectors(){detector_altitudes.clear();
			             detector_depths.clear();}
    
    inline void AddDetectorAtAltitude(G4double alt)
                                    {detector_altitudes.push_back(alt);}
    inline void AddDetectorAtDepth(G4double depth)
                                    {detector_depths.push_back(depth);}
     
    ///Set methods
    inline void SetAtmosphereHmax(G4double hmax){atmosphere_hmax = hmax;}
    inline void SetAtmosphereHmin(G4double hmin){atmosphere_hmin = hmin;}
    inline void SetDepthPercent(G4double aVal){depth_percent=aVal;}
    inline void SetMaxThickness(G4double thick_max){max_thickness = thick_max;}
    inline void SetMinThickness(G4double thick_min){min_thickness = thick_min;}
    inline void SetGeometryType(G4String astring){GeometryType = astring;}
    inline void SetHalfLength(G4double aVal){half_length = aVal;}
    inline void SetMagnetosphereH(G4double aVal){MagnetosphereH= aVal;}
    inline void SetPlanetH(G4double aVal){PlanetH = aVal;}
    inline void SetGeometryVerbosity(G4int val){geometry_verbosity = val;}
    inline void SetConsiderAtmosphere(G4bool aVal){ConsiderAtmosphere =aVal;}
    inline void SetDetectionInAtmosphereInMiddleOfALayer(G4bool aVal){ DetectionInAtmosphereInMiddleOfALayer=aVal;}
    inline void SetDetectionBelowSoilLayers(G4bool aVal){detect_below_soil_layers =aVal;}
     
     //Get methods
    
    inline G4double GetHalfLength(){return half_length;} 
    inline G4double GetAtmosphereHmax()const{return atmosphere_hmax;} 
    inline G4String GetGeometryType()const{return GeometryType;}
    inline G4int GetUpBoundaryIsADetector(G4int i){return UpBoundaryIsADetector[i];}
    inline G4int GetLowBoundaryIsADetector(G4int i){return LowBoundaryIsADetector[i];}
    
    //IsAnAtmoDetectionLayer was neede for a test porbably not more needed
    inline G4int GetIsAnAtmosphericDetectionLayer(G4int i){return IsAnAtmoDetectionLayer[i];}
    inline PlanetAtmosphereModel* GetAtmosphericModel(){return theAtmosphericModel;}
    inline std::vector<double> GetAtmoAltitudes(){return AtmoAltitudes;}
    inline std::vector<double> GetAtmoDepths(){return AtmoDepths;}
    inline std::vector<double> GetSoilDepthsInLengthUnit(){return SoilDepthsInLengthUnit;}
    inline std::vector<double> GetSoilDepthsInDepthUnit(){return SoilDepthsInDepthUnit;}
    
    
    
    
    
    inline std::vector<G4LogicalVolume*> GetAtmosphereLogicalVolumeVector()
                                                {return logicAtmosphere;}
    
    double GetAtmosphereMaxStepLength();
    double GetMagnetosphereMaxStepLength();
    double GetMagnetosphereMaxTrackLength();
    double GetMagnetosphereMaxTrackDuration();
    						
						
    G4int GetNbOfLayers();
    G4LogicalVolume* GetLogicalLayer(G4int n);
    
    //Elements
    
    void CreateElement(G4String el_name,
    		       G4String el_symbol,
		       G4double z,
		       G4double a);
    void ListElements();
     
    
    
     
  private:
    
    PLANETOCOSGeometryMessenger*  detectorMessenger; 
   
//Element table O, H, He,Ar,N 
//--------------
     
    G4Element* elO;
    G4Element* elH;
    G4Element* elC;
    G4Element* elHe;
    G4Element* elAr;
    G4Element* elN; 
    

    
//  Atmospheric model  
//-------------------------------------------

    PlanetAtmosphereModel*  theAtmosphericModel;

//geometry parameters
//---------------
     
    G4String GeometryType; //Spherical,Euclidean
    G4double atmosphere_hmax; // top of the atmosphere 
    G4double atmosphere_hmin; // bottom of the atmosphere 
    G4double MagnetosphereH; // Magnetosphere extension beyong atmosphere
    G4double PlanetH;  // Planet surface thickness for flat geometry
    G4double depth_percent;
    G4double half_length;  //half length of the different volume boxes 
                          // in x and y. Only for Euclidean geometry 
    G4double max_thickness; //maximum tickness of am atmospheric layer
    G4double min_thickness; //minimum tickness of an atmospheric layer
    G4double percent_depth; //percent of total atmospheric depth for a layer
    std::vector<double> detector_altitudes; 
    			//altitudes where flux should be detected in atmosphere 			       
    std::vector<double> detector_depths; //depth corresponding to the flux_detector_altitudes
    
    std::vector<double> detector_altitudes_for_histo; 			       
    std::vector<double> detector_depths_for_histo; 
    
    std::vector<G4Material*> AtmoMaterials;
    std::vector<double> AtmoAltitudes;
    std::vector<double> AtmoDepths;
    
    std::vector<double> SoilDepthsInLengthUnit;
    std::vector<double> SoilDepthsInDepthUnit;
    
    std::vector<int > UpBoundaryIsADetector;
    std::vector<int > LowBoundaryIsADetector;
    std::vector<int> IsAnAtmoDetectionLayer;
    
    G4double ZpositionForAtmosphereBottom; // coordinate of atmosphere bottom in world 
    G4double ZpositionAtZeroAltitude;
    std::vector< G4LogicalVolume* >   logicAtmosphere;  
    G4double* LimitsOfAtmosphereLayers;
    G4UserLimits* theMagnetosphereUserLimits;
    G4UserLimits* theAtmosphereUserLimits;
    G4Material* Vacuum; 
    G4LogicalVolume* logicMagnetosphere; 
    G4int geometry_verbosity;
    G4bool detect_below_soil_layers;
    
    
    //Atmosphere 
     
     
    //magnetic field
    PlanetMagneticField* fMagneticField;
     
    //sensitive detector for atmosphere
    PLANETOCOSSD* atmoSD;
    PLANETOCOSSoilSD* soilSD;
    
    
    //Atmosphere
    
    G4bool ConsiderAtmosphere;
    G4bool DetectionInAtmosphereInMiddleOfALayer;
    
    //MaxStepLength
    G4double MagnetosphereMaxStepLength;
    G4double AtmosphereMaxStepLength;

    
    
     						     
private:
    G4VPhysicalVolume* ReConstructGeometry(); 
    
         			    
};

#endif

