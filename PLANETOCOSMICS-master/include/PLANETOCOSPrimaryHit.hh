#ifndef PLANETOCOSPrimaryHit_h
#define PLANETOCOSPrimaryHit_h 1
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
class PLANETOCOSPrimaryHit : public G4VHit
{
public:
  PLANETOCOSPrimaryHit ();
  PLANETOCOSPrimaryHit (G4double aWeight, G4int aPDGCode,G4double anEnergy,
                    G4double aPhi, G4double aTheta, G4double LatitudeOrY,
                    G4double LongitudeOrX);
  ~PLANETOCOSPrimaryHit ();
  PLANETOCOSPrimaryHit (const PLANETOCOSPrimaryHit&);
  const PLANETOCOSPrimaryHit& operator= (const PLANETOCOSPrimaryHit&);
  int operator== (const PLANETOCOSPrimaryHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  inline G4double GetPhi(){return Phi;}
  inline G4double GetTheta(){return Theta;}
  inline G4double GetEnergy(){return Energy;}
  inline G4double GetWeight() {return Weight;}
  inline G4int GetPDGCode() {return PDGCode;}
  
  inline G4double GetLatitudeOrY(){return LatitudeOrY;}
  inline G4double GetLongitudeOrX(){return LongitudeOrX;}  
  
  void Draw () {};
  void Print () {};
  void clear () {};
  void DrawAll () {};
  void PrintAll () {};

private:
  G4double Weight;
  G4int PDGCode;
  G4double Energy;
  G4double Phi;
  G4double Theta;
  
  G4double LatitudeOrY;
  G4double LongitudeOrX;
};

typedef G4THitsCollection<PLANETOCOSPrimaryHit> PLANETOCOSPrimaryHitsCollection;

extern G4Allocator<PLANETOCOSPrimaryHit> PLANETOCOSPrimaryHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PLANETOCOSPrimaryHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) PLANETOCOSPrimaryHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PLANETOCOSPrimaryHit::operator delete(void *aHit)
{
  PLANETOCOSPrimaryHitAllocator.FreeSingle((PLANETOCOSPrimaryHit*) aHit);
}


////////////////////////////////////////////////////////////////////////////////
#endif
