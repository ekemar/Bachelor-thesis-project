//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolMessenger.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
#include"BlineToolMessenger.hh"
#include"BlineTool.hh"
#include"BlineToolEventAction.hh"

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

BlineToolMessenger::BlineToolMessenger
            (BlineTool* aBlineTool )
{ theBlineTool = aBlineTool;


 // G4UIparameter* param;
  
  
  
  BlineToolDir=new G4UIdirectory("/BlineTool/");
  BlineToolDir->SetGuidance("Interactive commands to trace and visualise magnetic field lines");
 
 
// commands

  BlineCmd = new G4UIcmdWithAnInteger("/BlineTool/ComputeBline",this);
  BlineCmd->SetGuidance("Compute magnetic field lines");
  BlineCmd->SetParameterName("nb_of_lines",false);
  BlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetMaxTrackingStepCmd = new G4UIcmdWithADoubleAndUnit("/BlineTool/SetMaxStepLength",this); 
  SetMaxTrackingStepCmd->SetGuidance("Set the maximum legth of tracking step when integrating magnetic field line");
  SetMaxTrackingStepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawColourCmd = new
             G4UIcmdWith3Vector("/BlineTool/SetColour",this);
  SetDrawColourCmd->SetGuidance("Set the colour drawing trajectories  and magnetic field lines");
  SetDrawColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  
  
  SetDrawBlineCmd = new
             G4UIcmdWithABool("/BlineTool/StockLines",this);
  SetDrawBlineCmd->SetGuidance("If true field lines are stocked in lines to be drawn ");
  SetDrawBlineCmd->SetParameterName("StockLines",false);
  SetDrawBlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawPointsCmd = new
             G4UIcmdWithABool("/BlineTool/StockPoints",this);
  SetDrawPointsCmd->SetGuidance("If true step field line points  are stocked in vector of points to be drawn");
  SetDrawPointsCmd->SetParameterName("StockPoints",false);
  SetDrawPointsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPointSizeCmd = new
             G4UIcmdWithADouble("/BlineTool/SetPointSize",this);
  SetPointSizeCmd->SetGuidance("Set the size of points for  drawing");
  SetPointSizeCmd->SetParameterName("StepSize",false);
  SetPointSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DrawCmd = new
             G4UIcmdWithoutParameter("/BlineTool/Show",this);
  DrawCmd->SetGuidance("Show the stored magnetic field  lines ");
  DrawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ResetCmd = new
             G4UIcmdWithoutParameter("/BlineTool/ResetMaterialToBeDrawn",this);
  ResetCmd->SetGuidance("Clear the vectors of lines and points to be drawn ");
  ResetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  
}

BlineToolMessenger::~BlineToolMessenger()
{
}		  

void BlineToolMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if (command == BlineCmd) theBlineTool->ComputeBlines(1);
  
  
  else if( command == SetMaxTrackingStepCmd ) 
            theBlineTool->SetMaxTrackingStep(SetMaxTrackingStepCmd
                                           ->GetNewDoubleValue(newValues));
  
  else if( command == SetDrawBlineCmd ) 
            theBlineTool->GetEventAction()->SetDrawBline(SetDrawBlineCmd
                              ->GetNewBoolValue(newValues));
  else if( command == SetDrawColourCmd ) 
        {G4ThreeVector vec=SetDrawColourCmd->GetNew3VectorValue(newValues);
	 theBlineTool->GetEventAction()->SetDrawColour(G4Colour(vec.x(),vec.y(),vec.z()));}			        	      
   
  else if( command == SetDrawPointsCmd ) 
            theBlineTool->GetEventAction()->SetDrawPoints(SetDrawPointsCmd
                              ->GetNewBoolValue(newValues));
   
  else if( command == SetPointSizeCmd ) 
            theBlineTool->GetEventAction()->SetPointSize(SetPointSizeCmd
                              ->GetNewDoubleValue(newValues));			      
  else if( command == DrawCmd ) 
                 theBlineTool->GetEventAction()->DrawFieldLines(.5,45.,45.);
  
  else if( command == ResetCmd ) theBlineTool->GetEventAction()->ResetVectorObjectToBeDrawn();
  
	    				      
}












