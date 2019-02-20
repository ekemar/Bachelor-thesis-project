#ifndef PLANETOCOSSteppingAction_h
#define PLANETOCOSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <map>
#include "G4ParticleTableIterator.hh"
#include "vector"



typedef G4ParticleTableIterator<G4String, G4double>::Map PartDictionary;
typedef G4ParticleTableIterator<G4String, G4double> PartDicIterator; 

class PLANETOCOSEventAction;
class PLANETOCOSSteppingActionMessenger;
class PLANETOCOSSteppingAction : public G4UserSteppingAction
{
  public:
    PLANETOCOSSteppingAction();
    virtual ~PLANETOCOSSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);
    
    void AddUntrackedParticle(G4String aParticleName, G4double aEkin);
    void RemoveUntrackedParticle(G4String);
    void ListUntrackedParticles();
    
    inline PartDictionary* GetUntrackedParticleDic() const
                                      {return UntrackedParticleDic;}
				      
    inline double GetLastLowestAltitude()const {return LastLowestAltitude;}
    inline double GetLowestAltitudeNeeded()const {return LowestAltitudeNeeded;}
    inline void SetLastLowestAltitude(G4double aVal){LastLowestAltitude = aVal;} 
    inline void SetDetectFlux(G4bool aBool){DetectFlux=aBool;}
    inline void SetStopUpFluxAtSelectedBoundary(G4bool aBool)
    			           {StopUpFluxAtSelectedBoundary=aBool;}
    inline void SetStopDownFluxAtSelectedBoundary(G4bool aBool)
    			           {StopDownFluxAtSelectedBoundary=aBool;}
    inline void SetNameOfStopBoundaryForUpFlux(G4String aName) 
    				{NameOfStopBoundaryForUpFlux=aName;}
    inline void SetNameOfStopBoundaryForDownFlux(G4String aName) 
    				{NameOfStopBoundaryForDownFlux=aName;}
    inline void SetPrimaryAlreadyInMagnetosphere(G4bool aBool)
                                   {PrimaryAlreadyInMagnetosphere =aBool;} 
    inline void SetStopAtMagnetopause(G4bool aBool){StopAtMagnetopause=aBool;}
    inline void SetMagnetopauseOutFactor(G4double aVal)
                                     {MagnetopauseOutFactor=aVal;}
    inline void SetMaxNumberOfTurnAroundThePlanet(G4double aVal)
    				     {MaxNumberOfTurnAroundThePlanet =aVal;}
    inline G4double GetMaxNumberOfTurnAroundThePlanet(){return MaxNumberOfTurnAroundThePlanet;}
    inline void SetStopAltitudeForUpward(G4double aVal)
    				     {StopAltitudeForUpward =aVal;}
    inline void SetStopAltitudeForDownward(G4double aVal)
    				     {StopAltitudeForDownward =aVal;}
    inline void SetStopDownwardPrimary(G4bool aBool)
    				     {StopDownwardPrimary =aBool;} 
    inline void SetStopUpwardPrimary(G4bool aBool)
    				     {StopUpwardPrimary =aBool;}
    
    inline void SetpAltitudes(std::vector<double >* apointer)
                                      {pAltitudes = apointer;}				     

#ifdef ANALYSIS_SOIL

    inline void  ResetLowestZorRReached(){LowestZorRReached = 1e14;}	
    inline double GetLowestZorRReached() const {return  LowestZorRReached;}	
#endif
  private:
    
    PLANETOCOSEventAction* eventAction;
    PartDictionary* UntrackedParticleDic;    
    PLANETOCOSSteppingActionMessenger* theMessenger; 
    G4bool DetectFlux ;
    G4int last_trackID;
    G4int last_stepNb;
    G4int nbad, nStep;
    
    
    G4bool StopUpFluxAtSelectedBoundary;
    G4bool StopDownFluxAtSelectedBoundary;
    G4bool LastStepWasOnBoundary;
    G4String NameOfStopBoundaryForUpFlux;
    G4String NameOfStopBoundaryForDownFlux;
   
    G4bool PrimaryAlreadyInMagnetosphere;
    G4bool StopAtMagnetopause;
    G4double MagnetopauseOutFactor;
    G4double LastLowestAltitude;
    G4double LowestAltitudeNeeded;
    G4double MaxNumberOfTurnAroundThePlanet;
    
    
    G4double StopAltitudeForUpward;
    G4bool StopUpwardPrimary;
    
    G4bool StopDownwardPrimary;
    G4double StopAltitudeForDownward;
    
    
    std::vector< double > * pAltitudes; //pointer on Altitudes in geometry
   
    

#ifdef ANALYSIS_SOIL

    G4double LowestZorRReached;	

#endif
   
    
};

#endif
