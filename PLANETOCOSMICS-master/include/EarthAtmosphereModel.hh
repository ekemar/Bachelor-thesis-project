
#ifndef EarthAtmosphereModel_h
#define EarthAtmosphereModel_h 1
#include "globals.hh"
#include"G4strstreambuf.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "DateAndTime.hh"
#include "PlanetAtmosphereModel.hh"

class EarthAtmosphereMessenger;
class EarthAtmosphereModel : public PlanetAtmosphereModel
{
  public:
    EarthAtmosphereModel();
    ~EarthAtmosphereModel();

  public:      
     
    
    
    //Set methods	
    
    inline void SetReferenceDate(DateAndTime aDate ){ReferenceDate =aDate;}
    inline void SetGEODETICLongitude(G4double aVal)
                                       {GEODETICLong = (float) aVal ;}
    inline void SetGEODETICLatitude(G4double aVal)
                                       {GEODETICLat = (float) aVal ;}
    
    inline void SetAp(G4double aVal){Ap = (float) aVal;}
    inline void SetF107(G4double aVal){f107 = (float) aVal;}
    inline void SetF107A(G4double aVal){f107_a = (float) aVal;}	
    void SetDefault();		   
      
     
  private:
     
    //Messenger
    
    EarthAtmosphereMessenger* myMessenger;
    //MSIS model parameters
    //---------------------
    DateAndTime ReferenceDate; //to be redefined
    float GEODETICLong, GEODETICLat;
    float Ap, f107, f107_a;
     
    
  private:
  
    virtual void ComputeAtmosphereTableFromModel(G4double alt_min,G4double alt_max);				   
    
    
    
      				    
};

#endif

