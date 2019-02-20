#ifndef PLANETUNITS_HH
#define PLANETUNITS_HH 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              MAGCOSUnits.hh
//
// Version:		VERSION_NUMBER
// Date:		LAST_DATE
// Author:		L. Desorgher
// Organisation:	ORGANISATION_NAME
// Project:		PROJECT_NAME
// Customer:		CUSTOMER_NAME
// Contract:		CONTRACT_NAME
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This interface define new units of length and magtnetic field used in
// magnetocosmics. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "globals.hh"
//Earth
static const G4double re=6371.2*km;
static const G4double Re=6371.2*km;
static const G4double gearth=1.; //in earth gravity
static const G4double nT=1.e-9*tesla;
static const G4double nanotesla=1.e-9*tesla;
static const G4double GV=GeV;
static const G4double MV=MeV;
static const G4double erg =1e-7 *joule;
//Mercury
static const G4double rmerc=2439.7*km;
static const G4double Rmerc=2439.7*km;
static const G4double gmerc=0.38;//in earth gravity
static const G4double au =1.4959787e8*km;
//Mars
static const G4double Rmars =3397.*km;
static const G4double rmars =3397.*km;
static const G4double gmars =0.38;//in earth gravity
//Jupiter
static const G4double Rjupiter =71492.*km;
static const G4double rjupiter =71492.*km;
static const G4double gjupiter =2.36;//in earth gravity




#endif
