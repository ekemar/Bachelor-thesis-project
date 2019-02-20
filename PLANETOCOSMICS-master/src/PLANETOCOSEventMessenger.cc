#include "G4RunManager.hh"
#include "PLANETOCOSEventMessenger.hh"
#include "PLANETOCOSEventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
// #include "G4VViewer.hh"

PLANETOCOSEventMessenger::PLANETOCOSEventMessenger(PLANETOCOSEventAction * myAct)
:myAction(myAct)
{ 

  myDrawDir = new G4UIdirectory("/PLANETOCOS/DRAW/");
  myDrawDir->SetGuidance("Define parameters for drawing trajectories");
  
 //Unfortunately line width and line style are not yet consider 
 // in G4SceneHandler when  drawing G4Polyline 
 
 /* SetDrawLineWidthCmd = new
             G4UIcmdWithADouble("/PLANETOCOS/DRAW/SetLineWidth",this);
  SetDrawLineWidthCmd->SetGuidance("Set the  line width for drawing trajectories  and magnetic field lines");
  SetDrawLineWidthCmd->SetParameterName("LineWidth",false);
  SetDrawLineWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawLineStyleCmd = new
             G4UIcmdWithAnInteger("/PLANETOCOS/DRAW/SetLineStyle",this);
  SetDrawLineStyleCmd->SetGuidance("Set the  line Style for drawing trajectories  and magnetic field lines");
  SetDrawLineStyleCmd->SetParameterName("LineStyle",false);
  SetDrawLineStyleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
  
  SetDrawColourForBlineCmd = new
             G4UIcmdWith3Vector("/PLANETOCOS/DRAW/SetColourForBline",this);
  SetDrawColourForBlineCmd->SetGuidance("Set the colour used for drawing  magnetic field lines");
  SetDrawColourForBlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  SetDrawingCoordinateSystemCmd = new
             G4UIcmdWithAString("/PLANETOCOS/DRAW/SetCoordinateSystem",this);
  SetDrawingCoordinateSystemCmd->SetGuidance("Set the reference coordinate system for drawing trajectories and field line");
  SetDrawingCoordinateSystemCmd->SetParameterName("LineWidth",false);
  SetDrawingCoordinateSystemCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  
  SetDrawTrajectoryCmd = new
             G4UIcmdWithABool("/PLANETOCOS/DRAW/DrawTrajectory",this);
  SetDrawTrajectoryCmd->SetGuidance("If true trajectories and field line are traced");
  SetDrawTrajectoryCmd->SetParameterName("DrawTrajectory",false);
  SetDrawTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawPointsCmd = new
             G4UIcmdWithABool("/PLANETOCOS/DRAW/DrawPoints",this);
  SetDrawPointsCmd->SetGuidance("If true step points are drawn on trajectories and field line with are traced");
  SetDrawPointsCmd->SetParameterName("DrawPoints",false);
  SetDrawPointsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPointSizeCmd = new
             G4UIcmdWithADouble("/PLANETOCOS/DRAW/SetPointSize",this);
  SetPointSizeCmd->SetGuidance("Set the size of step points for  drawing");
  SetPointSizeCmd->SetParameterName("StepSize",false);
  SetPointSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DrawCmd = new
             G4UIcmdWithoutParameter("/PLANETOCOS/DRAW/Show",this);
  DrawCmd->SetGuidance("Show the stored particle trajectories and magnetic field  lines ");
  DrawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ResetCmd = new
             G4UIcmdWithoutParameter("/PLANETOCOS/DRAW/Reset",this);
  ResetCmd->SetGuidance("Show the stored particle trajectories and magnetic field  lines ");
  ResetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
 /* G4String cmd_name = "/PLANETOCOS/DRAW/MagnetopauseLine";
  TraceMagnetopauseLineCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  TraceMagnetopauseLineCmd->SetGuidance("Trace a line on the magnetopause");
  TraceMagnetopauseLineCmd->SetParameterName("theta",true);
  TraceMagnetopauseLineCmd->SetDefaultValue(0.);
  TraceMagnetopauseLineCmd->SetUnitCategory("Angle");
  TraceMagnetopauseLineCmd->SetDefaultUnit("degree"); */
  
  AddParticleToBeDrawnCmd = new
             G4UIcommand("/PLANETOCOS/DRAW/SetParticleTrajectoryColour",this);
  AddParticleToBeDrawnCmd
    ->SetGuidance("Define the colour that will be used to draw the next computed trajectories of a given type of particle. ");
 AddParticleToBeDrawnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G4UIparameter* param;
  param = new G4UIparameter("Pname",'s',false);
  param->SetDefaultValue("e-");
  AddParticleToBeDrawnCmd->SetParameter(param);
  
  param = new G4UIparameter("R",'d',false);
  param->SetDefaultValue("1.0");
  AddParticleToBeDrawnCmd->SetParameter(param);
  
  param = new G4UIparameter("G",'d',true);
  param->SetDefaultValue("0.0");
  AddParticleToBeDrawnCmd->SetParameter(param);
  
  param = new G4UIparameter("B",'d',true);
  param->SetDefaultValue("0.0");
  AddParticleToBeDrawnCmd->SetParameter(param);
  
  
  RemoveParticleToBeDrawnCmd = new
             G4UIcmdWithAString("/PLANETOCOS/DRAW/DoNotDrawParticleTrajectory",this);
  RemoveParticleToBeDrawnCmd
    ->SetGuidance("Remove the given particle from the list of particles for wich the trajectories should be registered for a later visualisation. ");
  RemoveParticleToBeDrawnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
 
  
  
}

PLANETOCOSEventMessenger::~PLANETOCOSEventMessenger()
{ delete   SetDrawTrajectoryCmd;
  delete   SetDrawColourForBlineCmd;
  delete   SetDrawingCoordinateSystemCmd;
  delete   SetDrawPointsCmd;
  delete   SetPointSizeCmd;
 
}

void PLANETOCOSEventMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  
 /* if( command == SetDrawLineWidthCmd ) 
            myAction->SetDrawLineWidth(SetDrawLineWidthCmd
                              ->GetNewDoubleValue(newValues));
   if( command == SetDrawLineStyleCmd ) 
          {G4int i=SetDrawLineStyleCmd
                              ->GetNewIntValue(newValues);
	    if (i<3)		      
            myAction->SetDrawLineStyle(G4VisAttributes::LineStyle(i)); 
	  }	*/
   if( command == SetDrawTrajectoryCmd ) 
            myAction->SetDrawTrajectory(SetDrawTrajectoryCmd
                              ->GetNewBoolValue(newValues));
   else if( command == SetDrawingCoordinateSystemCmd ) 
            myAction->SetDrawingCoordinateSystem(newValues);
   
   
   else if( command == SetDrawColourForBlineCmd ) 
        {G4ThreeVector vec=SetDrawColourForBlineCmd->GetNew3VectorValue(newValues);
	 myAction->SetDrawColourForBline(G4Colour(vec.x(),vec.y(),vec.z()));}			        	      
   
   else if( command == SetDrawPointsCmd ) 
            myAction->SetDrawPoints(SetDrawPointsCmd
                              ->GetNewBoolValue(newValues));
   
   else if( command == SetPointSizeCmd ) 
            myAction->SetPointSize(SetPointSizeCmd
                              ->GetNewDoubleValue(newValues));			      
   else if( command == DrawCmd ) 
                 myAction->DrawTrajectoriesAndFieldLines(.5,45.,45.);
   else if( command == ResetCmd ) myAction->ResetVectorObjectToBeDrawn();
   
   
   /*else if (command == TraceMagnetopauseLineCmd){
         myAction->TraceMagnetopauseLine(TraceMagnetopauseLineCmd
	                                       ->GetNewDoubleValue(newValues)); 
	}*/ 
   else if( command == AddParticleToBeDrawnCmd ){
  	const char* paramString=newValues;
      	G4double R,G,B;
    	char pname[30];
       	std::istringstream is((char*)paramString);
      	is >> pname >> R  >> G >> B;
      	G4String particleName=pname;
	G4Colour aColour = G4Colour(R,G,B);
        myAction->AddParticleToBeDrawn(particleName,aColour);
  }
  
  else if( command == RemoveParticleToBeDrawnCmd ) 
            myAction->RemoveParticleToBeDrawn(newValues);
  	
 		       
}

