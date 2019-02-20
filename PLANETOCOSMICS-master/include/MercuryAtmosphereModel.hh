#ifndef MercuryAtmosphereModel_h
#define MercuryAtmosphereModel_h 1

#include "globals.hh"
#include"G4strstreambuf.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "DateAndTime.hh"
#include "PlanetAtmosphereModel.hh"


class MercuryAtmosphereModel : public PlanetAtmosphereModel
{
  public:
     MercuryAtmosphereModel();
    ~MercuryAtmosphereModel();  
     
     void SetDefault();

  private:
  
    virtual void ComputeAtmosphereTableFromModel(G4double alt_min,G4double alt_max);				   
    
    
    
      				    
};

#endif

