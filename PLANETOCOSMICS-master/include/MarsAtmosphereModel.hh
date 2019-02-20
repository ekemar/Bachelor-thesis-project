#ifndef MarsAtmosphereModel_h
#define MarsAtmosphereModel_h 1

#include "globals.hh"
#include"G4strstreambuf.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "DateAndTime.hh"
#include "PlanetAtmosphereModel.hh"


class MarsAtmosphereModel : public PlanetAtmosphereModel
{
  public:
     MarsAtmosphereModel();
    ~MarsAtmosphereModel();  
     
     void SetDefault();

  private:
  
    virtual void ComputeAtmosphereTableFromModel(G4double alt_min,G4double alt_max);				   
    
    
    
      				    
};

#endif

