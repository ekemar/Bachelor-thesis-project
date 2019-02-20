#ifndef PLANETOCOSSoilSD_h
#define PLANETOCOSSoilSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "PLANETOCOSEdepHit.hh"

class PLANETOCOSGeometryConstruction;
class G4Step;
class G4HCofThisEvent;
class G4StepPoint;


#include"G4ios.hh"

class PLANETOCOSSoilSD : public G4VSensitiveDetector
{

  public:
      PLANETOCOSSoilSD(G4String name, PLANETOCOSGeometryConstruction*);
      ~PLANETOCOSSoilSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      
      void clear();
      void DrawAll();
      void PrintAll();
      
      //methods
      void RegisterEkinOfKilledParticle(const G4Track* aTrack);
      
      //SetMethods
      inline void SetZGroundLevel(G4double aVal){ZGroundLevel = aVal;}
                                                   
      inline void SetpDepthsInDepthUnit(std::vector<double >* apointer)
                                      {pDepthsInDepthUnit = apointer;}
				      
      inline void SetpDepthsInLengthUnit(std::vector<double >* apointer)
                         	      {pDepthsInLengthUnit = apointer;}		          
      inline void SetDetectEdepVsThickness(G4bool abool)
                                            {DetectEdepVsThickness = abool;}
      inline void SetDetectEdepVsDepth(G4bool abool)
                                            {DetectEdepVsDepth = abool;}
      
      inline void SetGeometryType(G4String geo_type)
                                            {GeometryType = geo_type;}
					    
					    
      inline void SetNbLayersAboveSoil(G4int n){nb_layers_above_soil=n;}				    
      //GetMethods
      inline PLANETOCOSEdepHitsCollection* GetSoilEdepHitCollection()
      							{return SoilEdepHitCollection;}				  
  private:
      PLANETOCOSEdepHitsCollection* SoilEdepHitCollection; 
      G4double ZGroundLevel;
      PLANETOCOSGeometryConstruction* fDetector;
      G4String GeometryType;
      std::vector< double > * pDepthsInDepthUnit; //pointer on altitude vector of soil layers
      std::vector< double > * pDepthsInLengthUnit; // pointer on depth vector of soil layers
      G4bool DetectEdepVsDepth;
      G4bool DetectEdepVsThickness;
      G4int nb_layers_above_soil;
 
       
 
};




#endif

