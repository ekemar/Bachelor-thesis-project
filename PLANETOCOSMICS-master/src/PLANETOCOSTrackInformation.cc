#include "PLANETOCOSTrackInformation.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSTrackInformation::PLANETOCOSTrackInformation (G4bool aBool)
{WasProducedInTheCusp = aBool;
 WasAlreadyInMagnetosphere = false; 
 NbOfInwardCrossingOfLastDetector = 0;
 NbOfOutwardCrossingOfLastDetector =0;
 NbOfPlanetTurn =0.;
 NbOfCrossingOfEquator =0;
 HasBeenAlreadyUpwardDetected = false;
#ifdef TEST_ELEC_AT_BOUNDARY
 nb_last_detector =-1;
 costh_at_last_detector = -2.;
 nstep_at_last_detector=-1;
#endif 
 
}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSTrackInformation::~PLANETOCOSTrackInformation ()
{;}
