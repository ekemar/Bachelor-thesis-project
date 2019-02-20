#ifndef PLANETOCOSFluxHit_h
#define PLANETOCOSFluxHit_h 1
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
class PLANETOCOSFluxHit : public G4VHit
{
public:
  PLANETOCOSFluxHit ();
  PLANETOCOSFluxHit (G4double aWeight, G4int aPDGCode,G4double anEnergy,
                    G4double aPhi, G4double aTheta, G4double aLatOrY,
		    G4double aLongOrX, G4int nBoundDet);
  ~PLANETOCOSFluxHit ();
  PLANETOCOSFluxHit (const PLANETOCOSFluxHit&);
  const PLANETOCOSFluxHit& operator= (const PLANETOCOSFluxHit&);
  int operator== (const PLANETOCOSFluxHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  inline G4double GetPhi(){return Phi;}
  inline G4double GetTheta(){return Theta;}
  inline G4double GetLongitudeOrX(){return LongitudeOrX;}
  inline G4double GetLatitudeOrY(){return LatitudeOrY;}
  inline G4double GetEnergy(){return Energy;}
  G4double GetWeight() {return Weight;}
  G4int GetnBoundaryDetector() {return nBoundaryDetector;}
  G4int GetPDGCode() {return PDGCode;}
 
  
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
  G4int nBoundaryDetector; //detection layer number if -1 correspond to primary flux
  G4double LatitudeOrY;
  G4double LongitudeOrX;
    
};

typedef G4THitsCollection<PLANETOCOSFluxHit> PLANETOCOSFluxHitsCollection;

extern G4Allocator<PLANETOCOSFluxHit> PLANETOCOSFluxHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PLANETOCOSFluxHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) PLANETOCOSFluxHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PLANETOCOSFluxHit::operator delete(void *aHit)
{
  PLANETOCOSFluxHitAllocator.FreeSingle((PLANETOCOSFluxHit*) aHit);
}


////////////////////////////////////////////////////////////////////////////////
#endif
