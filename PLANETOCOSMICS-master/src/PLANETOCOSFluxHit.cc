////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSFluxHit.hh"

G4Allocator<PLANETOCOSFluxHit> PLANETOCOSFluxHitAllocator;
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSFluxHit::PLANETOCOSFluxHit ():Weight(1.),PDGCode(0),Energy(0.),Phi(0.),Theta(0.),nBoundaryDetector(-2)
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSFluxHit::PLANETOCOSFluxHit (G4double aWeight, G4int aPDGCode,G4double anEnergy,
                              G4double aPhi, G4double aTheta,  G4double aLatOrY,
		              G4double aLongOrX, G4int nBoundDet)
{ Weight  = aWeight;
  PDGCode = aPDGCode;
  Energy  = anEnergy;
  Phi     = aPhi;
  Theta   = aTheta;
  LatitudeOrY = aLatOrY;
  LongitudeOrX = aLongOrX;
  nBoundaryDetector = nBoundDet;
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSFluxHit::~PLANETOCOSFluxHit ()
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSFluxHit::PLANETOCOSFluxHit (const PLANETOCOSFluxHit& right) : G4VHit()
{
  Energy     = right.Energy;
  Weight   = right.Weight;
  Phi = right.Phi;
  Theta = right.Theta;
  LatitudeOrY = right.LatitudeOrY;
  LongitudeOrX = right.LongitudeOrX;
  nBoundaryDetector   = right.nBoundaryDetector;
  PDGCode = right.PDGCode;
}
////////////////////////////////////////////////////////////////////////////////
//
const PLANETOCOSFluxHit& PLANETOCOSFluxHit::operator= (const PLANETOCOSFluxHit& right)
{ Energy     = right.Energy;
  Weight   = right.Weight;
  Phi = right.Phi;
  Theta = right.Theta;
  LatitudeOrY = right.LatitudeOrY;
  LongitudeOrX = right.LongitudeOrX;
  nBoundaryDetector   = right.nBoundaryDetector;
  
  PDGCode = right.PDGCode;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//
int PLANETOCOSFluxHit::operator== (const PLANETOCOSFluxHit& ) const
{
  return 0;
}
