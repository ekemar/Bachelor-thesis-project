//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolEventAction.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines the EventAction used during tracing of magnetic field
//  lines by the BlineTool a tool 
// During the EndOfEventAction it stores magneticfield lines as polyline 
//  and Polymarker object.
// These polyline and Polymarker objects can be drawn  by using the drawn the         
// DrawFieldLines() method
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef BlineToolEVENTACTION_HH
#define BlineToolEVENTACTION_HH 
#include "G4UserEventAction.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

class G4Event;
class BlineToolEventMessenger;
class G4Polyline;
class G4Polymarker;
class BlineTool;

class BlineToolEventAction : public G4UserEventAction
{
  public:

   BlineToolEventAction(BlineTool* aBlineTool);
   ~BlineToolEventAction();

  public:
   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);
   
   void DrawFieldLines(G4double zoom, G4double theta, G4double phi); 
   void ResetVectorObjectToBeDrawn();
  private:
  
   BlineTool* fBlineTool;
   
   
   G4Colour DrawColour;
   G4bool DrawBline;
   G4bool DrawPoints;
   G4double PointSize;
   std::vector<G4VisAttributes*> TrajectoryVisAttributes;
   std::vector<G4Polyline>  TrajectoryPolyline;
   std::vector<G4Polymarker>  TrajectoryPoints;
   
    
    
 //inline methods
  public:
   inline void SetDrawColour(G4Colour aColour){ DrawColour = aColour;}
   
   // not yet valid could be in the future ????
  /*inline void SetDrawLineWidth(G4double aVal){
                 TrajectoryVisAttributes.SetLineWidth(aVal);}
   // not yet valid could be in the future ????
   inline void SetDrawLineStyle(G4VisAttributes::LineStyle aStyle){
                 TrajectoryVisAttributes.SetLineStyle(aStyle);}	*/	  		 
   
   inline void SetDrawBline(G4bool aBool){DrawBline=aBool;}
   inline void SetDrawPoints(G4bool aBool){DrawPoints=aBool;}
   inline void SetPointSize(G4double aVal){PointSize=aVal;}
   
   inline G4bool GetDrawBline(){return DrawBline;}
   
		   
};

#endif

    
