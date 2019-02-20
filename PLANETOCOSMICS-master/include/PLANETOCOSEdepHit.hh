#ifndef PLANETOCOSEdepHit_h
#define PLANETOCOSEdepHit_h 1
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
class PLANETOCOSEdepHit : public G4VHit
{
public:
  PLANETOCOSEdepHit ();
  ~PLANETOCOSEdepHit ();
  PLANETOCOSEdepHit (const PLANETOCOSEdepHit&);
  const PLANETOCOSEdepHit& operator= (const PLANETOCOSEdepHit&);
  int operator== (const PLANETOCOSEdepHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  //SetMethods
  inline void SetWeight(G4double val) { Weight= val;};
  inline void SetEdep(G4double val)   { Edep= val;};
  inline void SetAltitude1(G4double val) { Altitude1= val;};
  inline void SetAltitude2(G4double val) { Altitude2= val;};
  inline void SetDepth1(G4double val) { Depth1= val;};
  inline void SetDepth2(G4double val) { Depth2= val;};
#ifdef  TEST_EDEP  
  inline void SetParticleTypeCode(G4int code) { ParticleTypeCode= code;}
#endif   

 
 //GetMethods
  inline G4double GetWeight() {return Weight;};
  inline G4double GetEdep() {return Edep;};
  inline G4double GetAltitude1() {return Altitude1;};
  inline G4double GetAltitude2() {return Altitude2;};
  inline G4double GetDepth1() {return Depth1;};
  inline G4double GetDepth2() {return Depth2;};
#ifdef  TEST_EDEP  
  inline int GetParticleTypeCode() { return ParticleTypeCode;}
#endif    
 

  void Draw () {};
  void Print () {};
  void clear () {};
  void DrawAll () {};
  void PrintAll () {};

private:
  G4double Edep;
  G4double Weight;
  G4double Altitude1,Altitude2;
  G4double Depth1,Depth2;
#ifdef  TEST_EDEP 
  G4int ParticleTypeCode ;
#endif   
 

};

typedef G4THitsCollection<PLANETOCOSEdepHit> PLANETOCOSEdepHitsCollection;

extern G4Allocator<PLANETOCOSEdepHit> PLANETOCOSEdepHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PLANETOCOSEdepHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) PLANETOCOSEdepHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PLANETOCOSEdepHit::operator delete(void *aHit)
{
  PLANETOCOSEdepHitAllocator.FreeSingle((PLANETOCOSEdepHit*) aHit);
}


////////////////////////////////////////////////////////////////////////////////
#endif
