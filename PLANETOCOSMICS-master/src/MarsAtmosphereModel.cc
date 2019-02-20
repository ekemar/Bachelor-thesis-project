#include "MarsAtmosphereModel.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"

#include "G4PropagatorInField.hh"


MarsAtmosphereModel::MarsAtmosphereModel():
PlanetAtmosphereModel("Mars")
{ ;
}
////////////////////////////////////////////////////////////////////////////////
//
MarsAtmosphereModel::~MarsAtmosphereModel()
{;
}
////////////////////////////////////////////////////////////////////////////////
//
void MarsAtmosphereModel::
           ComputeAtmosphereTableFromModel(G4double ,G4double )
{;  
}
////////////////////////////////////////////////////////////////////////////////
//
void MarsAtmosphereModel::SetDefault()
{ AtmosphericModel="TABLE";
  G4String nameFile = getenv("MARS_ATMO_DATA")+G4String("/marsgram_atmo_table_0N180E.txt");
  ReadAtmosphereComposition(nameFile);
 
}
