#ifndef PLANETOCOSVModularPhysicsList_h
#define PLANETOCOSVModularPhysicsList_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSVModularPhysicsList: public G4VModularPhysicsList
{

public:

  void CleanAllPhysics () {
    G4PhysConstVector::iterator itr;
    for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
      delete (*itr);
    }
    physicsVector->clear();
  };
  
};
////////////////////////////////////////////////////////////////////////////////
#endif
