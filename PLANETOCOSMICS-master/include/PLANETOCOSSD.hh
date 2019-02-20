#ifndef PLANETOCOSSD_h
#define PLANETOCOSSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "PLANETOCOSFluxHit.hh"
#include "PLANETOCOSPrimaryHit.hh"
#include "PLANETOCOSEdepHit.hh"
#include "PLANETOCOSPostTrackHit.hh"

class PLANETOCOSGeometryConstruction;
class G4Step;
class G4HCofThisEvent;
class G4StepPoint;


#include"G4ios.hh"

class PLANETOCOSSD : public G4VSensitiveDetector
{

  public:
      PLANETOCOSSD(G4String name, PLANETOCOSGeometryConstruction*);
      ~PLANETOCOSSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      
      void clear();
      void DrawAll();
      void PrintAll();
      
      //methods
      void RegisterEkinOfKilledParticle(const G4Track* aTrack);
      
      //SetMethods
      inline void SetZSeaLevel(G4double aVal){ZSeaLevel = aVal;}
                                                   
      inline void SetpAltitudes(std::vector<double >* apointer)
                                      {pAltitudes = apointer;}
				      
      inline void SetpDepths(std::vector<double >* apointer)
                         	      {pDepths = apointer;}		          
      inline void SetDetectEdepvsAltitude(G4bool abool)
                                            {DetectEdepvsAltitude = abool;}
      inline void SetDetectEdepvsDepth(G4bool abool)
                                            {DetectEdepvsDepth = abool;}
      
      inline void SetGeometryType(G4String geo_type)
                                            {GeometryType = geo_type;}
					    
      inline void SetPrimaryHit(PLANETOCOSPrimaryHit* aHit)
                                            {primaryHit = aHit;}
					    
      inline void SetNbMagnetoLayers(G4int n){nb_magneto_layers=n;}				    
      //GetMethods
      inline PLANETOCOSFluxHitsCollection* GetDetectorFluxCollection()
                                          		{return DetectorFluxHitCollection;}
      inline PLANETOCOSPrimaryHitsCollection* GetPrimaryFluxCollection()
                                          		{return PrimaryFluxHitCollection;}				  
      inline PLANETOCOSPostTrackHitsCollection* GetPostTrackHitCollection()
                                          		{return PostTrackHitCollection;}
      
      inline PLANETOCOSEdepHitsCollection* GetEdepHitCollection()
      							{return EdepHitCollection;}				  
  private:
      PLANETOCOSFluxHitsCollection* DetectorFluxHitCollection;
      PLANETOCOSPrimaryHitsCollection* PrimaryFluxHitCollection;
      PLANETOCOSEdepHitsCollection* EdepHitCollection; 
      PLANETOCOSPostTrackHitsCollection* PostTrackHitCollection;
      G4double ZSeaLevel;
      PLANETOCOSGeometryConstruction* fDetector;
      G4String GeometryType;
      std::vector< double > * pAltitudes; //pointer on atmospheric Altitudes in geometry
      std::vector< double > * pDepths; // pointer on atmospheric  depth in geometry
      G4bool DetectEdepvsAltitude;
      G4bool DetectEdepvsDepth;
      PLANETOCOSPrimaryHit* primaryHit;
      G4int nb_magneto_layers;
 
       
 
};




#endif

