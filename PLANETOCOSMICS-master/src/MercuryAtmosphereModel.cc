#include "MercuryAtmosphereModel.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"

#include "G4PropagatorInField.hh"


MercuryAtmosphereModel::MercuryAtmosphereModel():
PlanetAtmosphereModel("Mercury")
{ 
}
////////////////////////////////////////////////////////////////////////////////
//
MercuryAtmosphereModel::~MercuryAtmosphereModel()
{;
}
////////////////////////////////////////////////////////////////////////////////
//
void MercuryAtmosphereModel::
           ComputeAtmosphereTableFromModel(G4double ,G4double )
{;  
}
////////////////////////////////////////////////////////////////////////////////
//
void MercuryAtmosphereModel::SetDefault()
{ AtmosphericModel="TABLE";
 
}
