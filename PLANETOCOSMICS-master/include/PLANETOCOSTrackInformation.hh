#ifndef PLANETOCOSTrackInformation_h
#define PLANETOCOSTrackInformation_h 1
#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "globals.hh"

////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSTrackInformation : public G4VUserTrackInformation
{
  public:
     PLANETOCOSTrackInformation(G4bool aBool);
     ~PLANETOCOSTrackInformation();
 
  public:   
     
     virtual inline void Print() const{};
   
     inline void AddOneOutwardCrossingOfLastDetector()
     				{NbOfInwardCrossingOfLastDetector++;}
     inline void AddOneInwardCrossingOfLastDetector()
     				{NbOfOutwardCrossingOfLastDetector++;}
     inline void AddOneCrossingOfEquator()
     				{NbOfCrossingOfEquator++;}
     inline void SetWasAlreadyInMagnetosphere(G4bool aBool)
     				{WasAlreadyInMagnetosphere=aBool;}
     inline void SetLowestAltNeededForThisTrack(G4double val)
     				{LowestAltNeededForThisTrack =val;}			
     inline void SetHasBeenAlreadyUpwardDetected(G4bool aBool)
     				{HasBeenAlreadyUpwardDetected=aBool;} 
     inline void SetWasProducedInTheCusp(G4bool aBool)
     				{WasProducedInTheCusp=aBool;}
     
     inline void SetDetectionEnergy(G4double anEnergy)
     				{DetectionEnergy= anEnergy;}
				
#ifdef TEST_ELEC_AT_BOUNDARY 				
     inline void SetNbLastDetector(G4int n){nb_last_detector=n;}
     inline void SetNStepAtLastDetector(G4int n){nstep_at_last_detector=n;}
     inline void SetCosthAtLastDetector(G4double  costh){costh_at_last_detector=costh;}	
#endif     			
     inline G4int GetNbOfInwardCrossingOfLastDetector()
     				{return NbOfInwardCrossingOfLastDetector;}
     inline G4int GetNbOfOutwardCrossingOfLastDetector()
     				{return NbOfOutwardCrossingOfLastDetector;}
     inline G4int GetNbOfEquatorCrossing()
     				{return NbOfCrossingOfEquator;}
     inline void AddNbOfPlanetTurn(G4double aVal){NbOfPlanetTurn+=aVal;}
     inline G4double GetNbOfPlanetTurn() {return NbOfPlanetTurn;}
     inline G4bool GetWasProducedInTheCusp(){return WasProducedInTheCusp;}	
     inline G4bool GetWasAlreadyInMagnetosphere() {return WasAlreadyInMagnetosphere;}
     inline G4double GetLowestAltNeededForThisTrack()
     				{return LowestAltNeededForThisTrack;}
				 
     inline G4bool GetHasBeenAlreadyUpwardDetected()
     				{return HasBeenAlreadyUpwardDetected;}
     inline double GetDetectionEnergy(){return DetectionEnergy;}
#ifdef TEST_ELEC_AT_BOUNDARY     
     inline G4int GetNStepAtLastDetector(){return nstep_at_last_detector;}
     inline G4int GetNbLastDetector(){return nb_last_detector;}
     inline G4double GetCosthAtLastDetector(){return costh_at_last_detector;}				
#endif      
     
      private:
     G4int NbOfInwardCrossingOfLastDetector;
     G4int NbOfOutwardCrossingOfLastDetector;
     G4int NbOfCrossingOfEquator;
     G4double NbOfPlanetTurn; // positive/negative for eastward/westward turn
     G4bool WasProducedInTheCusp;
     G4bool WasAlreadyInMagnetosphere;
#ifdef TEST_ELEC_AT_BOUNDARY 
     G4int nb_last_detector;
     G4double costh_at_last_detector;
     G4int  nstep_at_last_detector;
#endif     
     //lowest altitude needed to obtain the actual track
     G4double LowestAltNeededForThisTrack;
     G4bool HasBeenAlreadyUpwardDetected;
     G4double DetectionEnergy;
     
    
  
};
////////////////////////////////////////////////////////////////////////////////
#endif
