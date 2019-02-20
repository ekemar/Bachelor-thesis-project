////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSPrimaryHit.hh"

G4Allocator<PLANETOCOSPrimaryHit> PLANETOCOSPrimaryHitAllocator;
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryHit::PLANETOCOSPrimaryHit ():Weight(1.),PDGCode(0),Energy(0.),Phi(0.),Theta(0.),LatitudeOrY(0.),LongitudeOrX(0.)
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryHit::PLANETOCOSPrimaryHit (G4double aWeight, G4int aPDGCode,G4double anEnergy,
                                    G4double aPhi, G4double aTheta,
				    G4double aLatitudeOrY, G4double aLongitudeOrX)
{ Weight  = aWeight;
  PDGCode = aPDGCode;
  Energy  = anEnergy;
  Phi     = aPhi;
  Theta   = aTheta;
  
  LatitudeOrY = aLatitudeOrY; 
  LongitudeOrX = aLongitudeOrX;  
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryHit::~PLANETOCOSPrimaryHit ()
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPrimaryHit::PLANETOCOSPrimaryHit (const PLANETOCOSPrimaryHit& right) : G4VHit()
{
  Energy     = right.Energy;
  Weight   = right.Weight;
  Phi = right.Phi;
  Theta = right.Theta;
  PDGCode = right.PDGCode;
  
  LatitudeOrY = right.LatitudeOrY;
  LongitudeOrX = right.LongitudeOrX;  
}
////////////////////////////////////////////////////////////////////////////////
//
const PLANETOCOSPrimaryHit& PLANETOCOSPrimaryHit::operator= (const PLANETOCOSPrimaryHit& right)
{ Energy     = right.Energy;
  Weight   = right.Weight;
  Phi = right.Phi;
  Theta = right.Theta;
  PDGCode = right.PDGCode;
  
  LatitudeOrY = right.LatitudeOrY;
  LongitudeOrX = right.LongitudeOrX;   
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//
int PLANETOCOSPrimaryHit::operator== (const PLANETOCOSPrimaryHit& ) const
{
  return 0;
}
