//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineTool.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines a tool to trace and visualise magnetic field lines
// To use this tool in your G4 application you should  
//                  create an instance of this class somewhere in your code
// It will only work if you have declared a G4MagneticField field object to the
//                                               FieldManager      
// 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// BlineEquation(G4MagneticField* MagField )
//    Constructor.
//
// ~BlineEquation()
//    Destructor.
//

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//

#ifndef BlineTool_h
#define BlineTool_h 1

#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4UserRunAction.hh"
#include"vector"
#include"G4VisAttributes.hh"
#include"G4Polyline.hh"
#include"G4Polymarker.hh"
#include"G4FieldManager.hh"
#include"G4IdentityTrajectoryFilter.hh"
#include"G4VUserPrimaryGeneratorAction.hh"
#include"G4MagneticField.hh" 
#include"G4ChordFinder.hh"
#include<vector>
#include"BlineToolSteppingAction.hh"

class BlineToolMessenger;
class BlineToolSteppingAction;
class BlineToolEventAction;
class BlineToolPrimaryGeneratorAction;
class BlineEquation;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4UserStackingAction;
class G4UserRunAction;
class G4UserEventAction;
class G4VUserPrimaryGeneratorAction;

class BlineTool : public G4UserRunAction 
{ public:
  
   BlineTool(G4VUserPrimaryGeneratorAction* aUserPrimaryGeneratorAction,
             G4UserStackingAction* aUserStackingAction = 0,
	     G4UserSteppingAction* aUserSteppingAction = 0,
	     G4UserTrackingAction* aUserTrackingAction = 0,
	     G4UserEventAction* aUserEventAction =0,
	     G4UserRunAction* aUserRunAction =0);
   virtual ~BlineTool();
   
   static BlineTool* GetInstance();
  
  public:
  
   virtual void BeginOfRunAction(const G4Run* aRun);
   virtual void EndOfRunAction(const G4Run* aRun);
  
   void ComputeBlines(G4int nlines);

//set methods
   inline void SetMaxTrackingStep(G4double max_step){MaxTrackingStep=max_step;}; 
 
 //get methods
   inline   BlineToolEventAction* GetEventAction(){return fEventAction;};
   
   //stop volume methods
   inline void 	AddStopVolume(G4String aVolumeName){
   	       		fSteppingAction->AddStopVolume(aVolumeName);	
	       	}
   inline void 	RemoveStopVolume(G4String aVolumeName){
   			fSteppingAction->RemoveStopVolume(aVolumeName);
		}
   inline void 	ClearStopVolumeVector(){
  			fSteppingAction->ClearStopVolumeVector();
		}
   		       
 private:
 //private methods
  // bool CheckMagneticFields();
   void ResetChordFinders(); 
 
 
 private:
 
   static BlineTool* instance;
 
 //atrributes
   
   BlineToolMessenger* fMessenger;
   BlineToolSteppingAction* fSteppingAction;
   BlineToolEventAction* fEventAction;
   BlineToolPrimaryGeneratorAction* fPrimaryGeneratorAction;
   G4double MaxTrackingStep;
   G4bool was_ResetChordFinders_already_called;
 
   //user defined primary generator action
   G4VUserPrimaryGeneratorAction* fUserPrimaryAction;
   G4UserRunAction* fUserRunAction;
   G4UserStackingAction* fUserStackingAction;
   G4UserSteppingAction* fUserSteppingAction;
   G4UserTrackingAction* fUserTrackingAction;
   G4UserEventAction* fUserEventAction;

   //       ChordFinders, detector fields, equation of motions, and field 
   //manager for the different locale and global magnetic fields      
   std::vector <G4ChordFinder* >   vecChordFinders;
   std::vector <G4FieldManager* > vecFieldManagers;
   std::vector <G4MagneticField* > vecMagneticFields;
   std::vector <BlineEquation*> vecEquationOfMotion;
   
   
   
   
   
   
   
   
   //Trajectory filter used to store intermediate points
   //---------------------------------------------------
  //no more used 
    //G4IdentityTrajectoryFilter* fTrajectoryFilter;
   
   // bline is a vector of G4ThreeVector position 
 //  std::vector< G4ThreeVector > aBline;
 
 
   
   
  
  
};












#endif
