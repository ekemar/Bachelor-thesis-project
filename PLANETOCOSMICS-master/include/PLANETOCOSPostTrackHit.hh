#ifndef PLANETOCOSPostTrackHit_h
#define PLANETOCOSPostTrackHit_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include <vector>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "G4HCofThisEvent.hh"

typedef std::vector<G4ThreeVector> Fluxhit;
////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSPostTrackHit : public G4VHit
{
public:
  PLANETOCOSPostTrackHit (G4double aWeight, G4int aPDGCode,
  		      G4double aStartEnergy,
		      G4double aStopEnergy,
		      G4double aLifeTime,
                      G4double aNbOfPlanetTurn, 
		      G4int aNbOfEquatorCrossing, 
		      G4int aNbOfOutwardCrossing, 
		      G4int aNbOfInwardCrossing);
  ~PLANETOCOSPostTrackHit ();
  PLANETOCOSPostTrackHit (const PLANETOCOSPostTrackHit&);
  const PLANETOCOSPostTrackHit& operator= (const PLANETOCOSPostTrackHit&);
  int operator== (const PLANETOCOSPostTrackHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  inline G4double GetStopEnergy(){return StopEnergy;}
  inline G4double GetStartEnergy(){return StartEnergy;}
  inline G4double GetNbOfPlanetTurn(){return NbOfPlanetTurn;}
  inline G4int GetNbOfEquatorCrossing(){return NbOfEquatorCrossing;}
  inline G4int GetNbOfOutwardCrossing(){return NbOfOutwardCrossing;}
  inline G4int GetNbOfInwardCrossing(){return NbOfInwardCrossing;}
   inline G4double GetLifeTime(){return LifeTime;}
  
  G4double GetWeight() {return Weight;}
  G4int GetPDGCode() {return PDGCode;}
  
  void Draw () {};
  void Print () {};
  void clear () {};
  void DrawAll () {};
  void PrintAll () {};

private:
  G4double Weight;
  G4int PDGCode;
  G4double StartEnergy;
  G4double StopEnergy;
  G4double LifeTime;
  G4double NbOfPlanetTurn;
  G4int NbOfEquatorCrossing;
  G4int NbOfOutwardCrossing;
  G4int NbOfInwardCrossing;
 
  

};

typedef G4THitsCollection<PLANETOCOSPostTrackHit> PLANETOCOSPostTrackHitsCollection;

extern G4Allocator<PLANETOCOSPostTrackHit> PLANETOCOSPostTrackHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PLANETOCOSPostTrackHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) PLANETOCOSPostTrackHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PLANETOCOSPostTrackHit::operator delete(void *aHit)
{
  PLANETOCOSPostTrackHitAllocator.FreeSingle((PLANETOCOSPostTrackHit*) aHit);
}


////////////////////////////////////////////////////////////////////////////////
#endif
