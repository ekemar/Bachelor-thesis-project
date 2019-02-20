////////////////////////////////////////////////////////////////////////////////
//
#include "PLANETOCOSPostTrackHit.hh"

G4Allocator<PLANETOCOSPostTrackHit> PLANETOCOSPostTrackHitAllocator;

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPostTrackHit::PLANETOCOSPostTrackHit (G4double aWeight, G4int aPDGCode,
  		                        G4double aStartEnergy,
		                        G4double aStopEnergy,
		                        G4double aLifeTime,
                                        G4double aNbOfPlanetTurn,
					G4int aNbOfEquatorCrossing, 
		                        G4int aNbOfOutwardCrossing, 
		                        G4int aNbOfInwardCrossing)
{ Weight=aWeight;
  PDGCode = aPDGCode;
  StartEnergy = aStartEnergy;
  StopEnergy = aStopEnergy;
  LifeTime = aLifeTime;
  NbOfPlanetTurn = aNbOfPlanetTurn; 
  NbOfEquatorCrossing = aNbOfEquatorCrossing;
  NbOfOutwardCrossing = aNbOfOutwardCrossing;
  NbOfInwardCrossing = aNbOfInwardCrossing; 
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPostTrackHit::~PLANETOCOSPostTrackHit ()
{}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSPostTrackHit::PLANETOCOSPostTrackHit (const PLANETOCOSPostTrackHit& right) : G4VHit()
{ Weight= right.Weight;
  PDGCode = right.PDGCode;
  StartEnergy = right.StartEnergy;
  StopEnergy = right.StopEnergy;
  LifeTime = right.LifeTime;
  NbOfPlanetTurn = right.NbOfPlanetTurn; 
  NbOfEquatorCrossing = right.NbOfEquatorCrossing;
  NbOfOutwardCrossing = right.NbOfOutwardCrossing;
  NbOfInwardCrossing = right.NbOfInwardCrossing; 
}
////////////////////////////////////////////////////////////////////////////////
//
const PLANETOCOSPostTrackHit& PLANETOCOSPostTrackHit::operator= (const PLANETOCOSPostTrackHit& right)
{ Weight= right.Weight;
  PDGCode = right.PDGCode;
  StartEnergy = right.StartEnergy;
  StopEnergy = right.StopEnergy;
  LifeTime = right.LifeTime;
  NbOfPlanetTurn = right.NbOfPlanetTurn;
  NbOfEquatorCrossing = right.NbOfEquatorCrossing; 
  NbOfOutwardCrossing = right.NbOfOutwardCrossing;
  NbOfInwardCrossing = right.NbOfInwardCrossing; 
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//
int PLANETOCOSPostTrackHit::operator== (const PLANETOCOSPostTrackHit& ) const
{
  return 0;
}
