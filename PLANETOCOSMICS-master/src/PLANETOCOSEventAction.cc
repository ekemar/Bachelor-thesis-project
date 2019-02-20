#include "PLANETOCOSEventAction.hh"
#include "PLANETOCOSEventMessenger.hh"
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
#include "SpaceCoordinatePlanet.hh"
#include "PlanetUnits.hh"
#include "PlanetMagneticField.hh"

#include "PLANETOCOSApplicationScenario.hh"
#include "G4UserRunAction.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4UserSteppingAction.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "PlanetManager.hh"
#include "DurationManager.hh"
#include "BlineTool.hh"



///////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEventAction::PLANETOCOSEventAction()
{ SetDrawColour(G4Colour(1.,0.,0.));
  SetDrawColourForBline(G4Colour(1.,0.,0.)); 
  DrawingCoordinateSystem="PLA";
  DrawTrajectory=true;
  DrawPoints=false;
  PointSize=1;
  theMessenger = new PLANETOCOSEventMessenger(this);
  vis_secondary_mode =false;
  ParticleToBeDrawn=new PartColourDictionary();
}
///////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSEventAction::~PLANETOCOSEventAction()
{
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::BeginOfEventAction(const G4Event* evt)
{ 
   DurationManager* theDurationManager = DurationManager::GetInstance();
   if (!theDurationManager->CheckDurationAtBeginOfEvent()) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The run will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
	
  	
   }
  
   theSteppingAction->SetPrimaryAlreadyInMagnetosphere(false);	
  	
   PLANETOCOSAnalysisManager::GetInstance()->BeginOfEventAction(evt);

}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::EndOfEventAction(const G4Event* evt)
{  
   SetVisSecondaryMode(false);
   //if secondary mode then vis_secondary_mode will be set to true in the next command
   PLANETOCOSAnalysisManager::GetInstance()->EndOfEventAction(evt);
 

   G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if(trajectoryContainer)  n_trajectories = trajectoryContainer->entries(); 
   
   const PLANETOCOSGeometryConstruction* theGeometry =
          dynamic_cast< const PLANETOCOSGeometryConstruction* > 
	  	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   							
   G4double ZpositionForAtmosphereBottom = theGeometry->GetZpositionForAtmosphereBottom();
   G4String GeometryType =  theGeometry->GetGeometryType();
    
   //visualisation
   //---------------------
   if (n_trajectories && (DrawTrajectory || DrawPoints) ){
       SpaceCoordinatePlanet* theCoordinateConvertor
                        = SpaceCoordinatePlanet::GetInstance();
			
       theCoordinateConvertor->SetSystemInAndOut("PLA",DrawingCoordinateSystem);
       //G4cout<<n_trajectories<<std::endl;
       for(G4int i=0; i<n_trajectories; i++){
      		G4Trajectory* aTrajectory = ( G4Trajectory*)
                    			(*(evt->GetTrajectoryContainer()))[i]; 
       		G4String aParticleName = aTrajectory ->GetParticleName();		  
       		G4int n_point = aTrajectory ->GetPointEntries();
		if (!vis_secondary_mode || ParticleToBeDrawn->find(aParticleName) 
						!= ParticleToBeDrawn->end() || aParticleName == "chargedgeantino"){
			G4Polyline pPolyline;
       			G4Polymarker stepPoints;
			G4Colour colour = DrawColour;
			if (vis_secondary_mode){
				colour=(*ParticleToBeDrawn)[aParticleName];
				if (aTrajectory->GetTrackID() ==1 &&
				   ParticleToBeDrawn->find("Primary")!= ParticleToBeDrawn->end())
				    colour=(*ParticleToBeDrawn)["Primary"];
			}
			if (aParticleName == "chargedgeantino")
					colour=	DrawColourForBline;
			TrajectoryVisAttributes.push_back(new G4VisAttributes(colour));
       			stepPoints.SetMarkerType(G4Polymarker::circles);
       			stepPoints.SetScreenSize(PointSize);
       			stepPoints.SetFillStyle(G4VMarker::filled);
       			stepPoints.SetVisAttributes(TrajectoryVisAttributes.back());
			bool AtLeastOneValidPoint =false;
			for(G4int ii=0; ii<n_point; ii++){ 
				G4ThreeVector pos_pla= 
          			((G4TrajectoryPoint*)
	       				aTrajectory->GetPoint(ii))
                                                      ->GetPosition();
				G4ThreeVector pos =pos_pla;
				G4double r_or_z=pos.z();
				if (GeometryType =="SPHERICAL"){ 		      
	  			            pos= theCoordinateConvertor->Transform(pos_pla);					      
          				    r_or_z= pos.mag();	  
				}
				
				
				G4bool DrawInSoil = false;
				if (r_or_z >= ZpositionForAtmosphereBottom ||  DrawInSoil){
					/*G4cout<<(r_or_z-ZpositionForAtmosphereBottom)/km<<std::endl;
					G4cout<<ZpositionForAtmosphereBottom/km<<std::endl;*/
					AtLeastOneValidPoint =true;
					if (DrawTrajectory)  pPolyline.push_back( pos);
					if (DrawPoints)  stepPoints.push_back(pos);
	    			}
			}
			if (AtLeastOneValidPoint){
       				pPolyline.SetVisAttributes(TrajectoryVisAttributes.back());
       				TrajectoryPolyline.push_back(pPolyline);
       				TrajectoryPoints.push_back(stepPoints);
			}	
       		}
    
    
   	}
	
	

 }

 SetVisSecondaryMode(false);
 
}
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::DrawTrajectoriesAndFieldLines
                                  (G4double, G4double, G4double)
{ 
  	       
  unsigned int nline = TrajectoryPolyline.size();
  unsigned int npoints =TrajectoryPoints.size();
     
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager){ 
   	G4cout<<"For visualising trajectories and magnetic";
       	G4cout<<" field lines you should select a visualisation driver"<<G4endl; 
       	return;
  }
  if (nline ==0){
   	G4cout<<"There is no particle trajectory to visualise"<<G4endl;
       	return;
  }
  ((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()-> ClearStore ();
     
  G4Scene* aScene = ((G4VisManager*)pVVisManager)
                            ->GetCurrentSceneHandler()->GetScene();
  if (!aScene){
  	aScene = new G4Scene();
	((G4VisManager*)pVVisManager)
                            ->GetCurrentSceneHandler()->SetScene(aScene);
  } 			    
  
  //aScene->Clear();
        
  //((G4VisManager*)pVVisManager)->GetCurrentScene()-> Clear(); 
  /* ((G4VisManager*)pVVisManager)->GetCurrentViewer()-> SetView ();
  ((G4VisManager*)pVVisManager)->GetCurrentViewer()-> ClearView ();*/
   
  G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/drawVolume");
   G4cout<<"Nlines "<<nline<<std::endl; 
  for (unsigned int i=0;i<nline;i++)
  	/*((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()
	                             ->AddPrimitive(TrajectoryPolyline[i]);*/
          pVVisManager->Draw(TrajectoryPolyline[i]);
  
  for (unsigned i=0;i<npoints;i++)
        	((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()
	                             	     ->AddPrimitive(TrajectoryPoints[i]); 	       
	       
      
 }
 
///////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::ResetVectorObjectToBeDrawn()
 {TrajectoryVisAttributes.clear();
  TrajectoryPolyline.clear();
  TrajectoryPoints.clear();
 }
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::AddParticleToBeDrawn(G4String aParticleName, 
					      G4Colour aColour)
{ G4ParticleDefinition* aParticleDefinition =
              G4ParticleTable::GetParticleTable()->FindParticle(aParticleName);
  if (aParticleDefinition || aParticleName =="Primary") (*ParticleToBeDrawn)[aParticleName]=aColour;
  else 
     G4cout <<" the  particle does not exist in the particle table "  <<G4endl;	    
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSEventAction::RemoveParticleToBeDrawn(G4String aParticleName)
{ if (ParticleToBeDrawn->find(aParticleName) != ParticleToBeDrawn->end())
    		ParticleToBeDrawn->erase(ParticleToBeDrawn->find(aParticleName));
  else
     	G4cout <<" the  particle was not in the list of untracked particles "  <<G4endl;
}

