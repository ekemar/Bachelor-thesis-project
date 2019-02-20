#include"MarsMagneticFieldMessenger.hh"
#include"MarsMagneticField.hh"
#include"G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//units
#include "G4UnitsTable.hh"

MarsMagneticFieldMessenger::MarsMagneticFieldMessenger (MarsMagneticField* aField )
:PlanetMagneticFieldMessenger(aField) 
{ 
 theField = aField;
}
////////////////////////////////////////////////////////////////////////////////
//
MarsMagneticFieldMessenger::~MarsMagneticFieldMessenger()
{;
} 
////////////////////////////////////////////////////////////////////////////////
//
void MarsMagneticFieldMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ bool result;
  result  = SetMotherNewValue(command, newValues);  
}
