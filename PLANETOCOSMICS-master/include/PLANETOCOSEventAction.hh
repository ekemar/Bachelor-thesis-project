#ifndef PLANETOCOSEVENTACTION_HH
#define PLANETOCOSEVENTACTION_HH 

#include "G4UserEventAction.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
class G4Event;
class PLANETOCOSEventMessenger;
class G4Polyline;
class G4Polymarker;
class PLANETOCOSSteppingAction;





class PLANETOCOSEventAction : public G4UserEventAction
{
  public:

   	PLANETOCOSEventAction();
   	~PLANETOCOSEventAction();

  public:
   	void BeginOfEventAction(const G4Event*);
   	void EndOfEventAction(const G4Event*);
   
   	typedef G4ParticleTableIterator<G4String, G4Colour>::Map PartColourDictionary;
   	typedef G4ParticleTableIterator<G4String, G4Colour> PartColourIterator;
   	void AddParticleToBeDrawn(G4String aParticleName, G4Colour aCol);
  	void RemoveParticleToBeDrawn(G4String aParticleName);
   
  	void DrawTrajectoriesAndFieldLines(G4double zoom, G4double theta, G4double phi); 
  	void ResetVectorObjectToBeDrawn();
//   	void TraceMagnetopauseLine(G4double theta);

  private:
   	PLANETOCOSEventMessenger* theMessenger;
	PLANETOCOSSteppingAction* theSteppingAction;
   	G4Colour DrawColour;
	G4Colour DrawColourForBline;
   	G4bool DrawTrajectory;
   	G4bool DrawPoints;
  	G4double PointSize;
   	G4String DrawingCoordinateSystem;
   	std::vector<G4VisAttributes*> TrajectoryVisAttributes;
   	std::vector<G4Polyline>  TrajectoryPolyline;
   	std::vector<G4Polymarker>  TrajectoryPoints;
   	PartColourDictionary* ParticleToBeDrawn;
   	G4bool vis_secondary_mode;
	
    
 //inline methods
  public:
        inline void SetDrawColourForBline(G4Colour aColour){ DrawColourForBline = aColour;}
   	inline void SetDrawColour(G4Colour aColour){ DrawColour = aColour;}
   	inline void SetDrawTrajectory(G4bool aBool){DrawTrajectory=aBool;}
   	inline void SetDrawPoints(G4bool aBool){DrawPoints=aBool;}
   	inline void SetPointSize(G4double aVal){PointSize=aVal;}
	inline void SetSteppingAction(PLANETOCOSSteppingAction* aSteppingAction)
					       {theSteppingAction = aSteppingAction;}
   	inline G4bool GetDrawTrajectory(){return DrawTrajectory;}
   
        inline void SetDrawingCoordinateSystem(G4String aCoordinateSystem){DrawingCoordinateSystem = aCoordinateSystem;}
   	inline void SetVisSecondaryMode(G4bool aBool){vis_secondary_mode = aBool;}			
 	 
  			
  
};

#endif

    
