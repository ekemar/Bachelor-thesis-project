//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolEventAction.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 

#include "BlineToolEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "BlineTool.hh"



BlineToolEventAction::BlineToolEventAction(BlineTool* aBlineTool)
{fBlineTool=aBlineTool;
}

BlineToolEventAction::~BlineToolEventAction()
{;
}

void BlineToolEventAction::BeginOfEventAction(const G4Event*)
{;}

void BlineToolEventAction::EndOfEventAction(const G4Event* evt)
{ 
  
  
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) 
   {n_trajectories = trajectoryContainer->entries(); 
   
   //visualisation
   //---------------------
  
    if (DrawBline || DrawPoints)
      {
       G4int n_point = 
             (*(evt->GetTrajectoryContainer()))[0]->GetPointEntries();
       
       G4Polyline pPolyline;
       G4Polymarker stepPoints;
       TrajectoryVisAttributes.push_back(new G4VisAttributes(DrawColour));
       stepPoints.SetMarkerType(G4Polymarker::circles);
       stepPoints.SetScreenSize(PointSize);
       stepPoints.SetFillStyle(G4VMarker::filled);
       stepPoints.SetVisAttributes(TrajectoryVisAttributes.back());
     
     
       for(G4int i=0; i<n_point; i++)
         { G4ThreeVector pos= 
          ((G4TrajectoryPoint*)
	       ((*(evt->GetTrajectoryContainer()))[0]->GetPoint(i)))
                                                      ->GetPosition();
	  if (DrawBline)  pPolyline.push_back( pos);
	  if (DrawPoints)  stepPoints.push_back(pos);
	 }
	
       
       pPolyline.SetVisAttributes(TrajectoryVisAttributes.back());
       
       TrajectoryPolyline.push_back(pPolyline);
       TrajectoryPoints.push_back(stepPoints);
       
        
       
     }
    
    
   }
} 

void BlineToolEventAction::DrawFieldLines
                                  (G4double , G4double , G4double
				  )
 {   unsigned int nline = TrajectoryPolyline.size();
     unsigned int npoints =TrajectoryPoints.size();
     G4cout<<"nline "<<nline<<std::endl;
     G4cout<<"npoints "<<npoints<<std::endl;
     
     
     G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
     if (!pVVisManager) 
      {G4cout<<"For visualising magnetic";
       G4cout<<" field lines you should select a visualisation driver"<<G4endl; 
       return;
      }
     
     if (nline ==0)
      {G4cout<<"There is no magnetic field line to visualise"<<G4endl;
       return;
      }
    
     
     for (unsigned int i=0;i<nline;i++)
               pVVisManager->Draw(TrajectoryPolyline[i]);
     for (unsigned int i=0;i<npoints;i++)
               pVVisManager->Draw(TrajectoryPoints[i]); 
      
    
    //((G4VisManager*)pVVisManager)->GetCurrentViewer()-> DrawView ();   
   //  ((G4VisManager*)pVVisManager)->GetCurrentViewer()->ShowView();
      
 }
void BlineToolEventAction::ResetVectorObjectToBeDrawn()
 {TrajectoryVisAttributes.clear();
  TrajectoryPolyline.clear();
  TrajectoryPoints.clear();
 }

