////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSEdepHit.hh"

G4Allocator<PLANETOCOSEdepHit> PLANETOCOSEdepHitAllocator;
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEdepHit::PLANETOCOSEdepHit ():Edep(0.),Weight(1.),Altitude1(0.),
 Altitude2(0.),Depth1(0.),Depth2(0.)
#ifdef  TEST_EDEP  
,ParticleTypeCode(0)
#endif
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEdepHit::~PLANETOCOSEdepHit ()
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEdepHit::PLANETOCOSEdepHit (const PLANETOCOSEdepHit& right) : G4VHit()
{
  Edep     = right.Edep;
  Weight   = right.Weight;
  Altitude1   =right.Altitude1;
  Altitude2   =right.Altitude2;
  Depth1 = right.Depth1;
  Depth2 =  right.Depth2;
#ifdef  TEST_EDEP  
  ParticleTypeCode=right.ParticleTypeCode;
#endif

}
////////////////////////////////////////////////////////////////////////////////
//
const PLANETOCOSEdepHit& PLANETOCOSEdepHit::operator= (const PLANETOCOSEdepHit& right)
{
  Edep     = right.Edep;
  Weight   = right.Weight;
  Altitude1   =right.Altitude1;
  Altitude2   =right.Altitude2;
  Depth1 = right.Depth1;
  Depth2 =  right.Depth2;
#ifdef  TEST_EDEP  
  ParticleTypeCode=right.ParticleTypeCode;
#endif
  
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//
int PLANETOCOSEdepHit::operator== (const PLANETOCOSEdepHit& ) const
{
  return 0;
}
