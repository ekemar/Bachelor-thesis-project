//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineTool.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
#include"BlineTool.hh"
#include"BlineToolMessenger.hh"
#include"BlineToolPrimaryGeneratorAction.hh"
#include"BlineToolEventAction.hh"
#include"BlineToolSteppingAction.hh"
#include"BlineEquation.hh"

#include"G4PropagatorInField.hh"
#include"G4TransportationManager.hh"
#include"G4FieldManager.hh"
#include"G4UImanager.hh"
#include"G4ParticleDefinition.hh"
#include"G4Circle.hh"
#include"G4Event.hh"
#include"G4UnitsTable.hh"
#include"G4ios.hh"
#include"G4RunManager.hh"
#include"G4CashKarpRKF45.hh"
#include"G4FieldManager.hh"
#include"G4LogicalVolumeStore.hh"

 
BlineTool* BlineTool::instance = 0;

//////////////////////////////////////////////////////////////////
BlineTool::BlineTool(G4VUserPrimaryGeneratorAction* aUserPrimaryGeneratorAction,
             G4UserStackingAction* aUserStackingAction ,
	     G4UserSteppingAction* aUserSteppingAction ,
	     G4UserTrackingAction* aUserTrackingAction ,
	     G4UserEventAction* aUserEventAction ,
	     G4UserRunAction* aUserRunAction)
{fUserPrimaryAction=aUserPrimaryGeneratorAction;
 fUserStackingAction=aUserStackingAction;
 fUserTrackingAction=aUserTrackingAction;
 fUserEventAction=aUserEventAction;
 fUserSteppingAction=aUserSteppingAction;
 fUserRunAction=aUserRunAction;

 fMessenger = new BlineToolMessenger(this);
 fSteppingAction = new BlineToolSteppingAction(this) ;
 fEventAction = new  BlineToolEventAction(this);
 fPrimaryGeneratorAction = new BlineToolPrimaryGeneratorAction();
 fPrimaryGeneratorAction->SetUserPrimaryAction(fUserPrimaryAction);
  //no more needed
  //fTrajectoryFilter = new G4IdentityTrajectoryFilter();
 MaxTrackingStep =10000000.*m;
 was_ResetChordFinders_already_called=false;
 
 
}
///////////////////////////////////////////////////////////////////////
BlineTool::~BlineTool()
{delete fMessenger;
 delete fSteppingAction;
 delete fEventAction; 
 delete fPrimaryGeneratorAction;
 for (unsigned int  i=0; i< vecEquationOfMotion.size();i++)
        { if (vecEquationOfMotion[i]) delete vecEquationOfMotion[i];
	  if (vecChordFinders[i]) delete vecChordFinders[i];
	}
 
}
////////////////////////////////////////////////////////////////////////////////
//
BlineTool* BlineTool::GetInstance()
{
  return instance;
}		  
////////////////////////////////////////////////////////////////////

void BlineTool::BeginOfRunAction(const G4Run*)
{ ;
}  
///////////////////////////////////////////////////////////////////////
void BlineTool::EndOfRunAction(const G4Run* )
{;
}
////////////////////////////////////////////////////////////////
void BlineTool::ComputeBlines(G4int n_of_lines)
{ //the first time ResetChordFinders should be called
   
   if (!was_ResetChordFinders_already_called)
       {ResetChordFinders();
        was_ResetChordFinders_already_called=true;
       } 


  //Replace the user action by the BlineTool actions
  
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  
  
  theRunManager->SetUserAction(this);
  theRunManager->SetUserAction(fSteppingAction);
  theRunManager->SetUserAction(fPrimaryGeneratorAction);
  theRunManager->SetUserAction(fEventAction);
  G4UserTrackingAction* aNullTrackingAction =0;
  theRunManager->SetUserAction(aNullTrackingAction);
  G4UserStackingAction* aNullStackingAction =0;
  theRunManager->SetUserAction(aNullStackingAction);
  
  
 
  
  
// replace the user defined chordfinder by  the element of 
//                                                  vecChordFinders  
  
  std::vector< G4ChordFinder*> user_chord_finders;
  std::vector<double>  user_largest_acceptable_step;
  for (unsigned int i=0;i<vecChordFinders.size();i++)
     {user_largest_acceptable_step.push_back(-1.);
      if (vecChordFinders[i])
       {user_chord_finders.push_back(vecFieldManagers[i]->GetChordFinder());
        vecChordFinders[i]->SetDeltaChord(user_chord_finders[i]
	                                              ->GetDeltaChord());
        vecFieldManagers[i]->SetChordFinder(vecChordFinders[i]);
		   
       }
      else user_chord_finders.push_back(0); 
     }
     
   // I have tried to use the smooth line filter ability  but I could not obtain 
   // a smooth trajectory in the G4TrajectoryContainer after an event
   // Another solution for obtaining a smooth trajectory is to limit
   // the LargestAcceptableStep in the G4PropagatorInField object. This is the
   // solution I used. 
	
	   //old solution
           // G4TransportationManager::GetTransportationManager()
           //          ->GetPropagatorInField()->SetTrajectoryFilter(fTrajectoryFilter);
   
          //new solution
          // set the largest_acceptable_step to max_step:length
       
        G4double previous_largest_acceptable_step =
        G4TransportationManager::GetTransportationManager()
                      ->GetPropagatorInField()
	                    ->GetLargestAcceptableStep();
         
	G4TransportationManager::GetTransportationManager()
                      ->GetPropagatorInField()
	                    ->SetLargestAcceptableStep(MaxTrackingStep); 
 	G4cout<<"MaxTrackingStep "<<MaxTrackingStep/km<<std::endl;
    
    //Start the integration of n_of_lines different magnetic field lines
    
     for (int i=0; i<n_of_lines;i++)
      {// for each magnetic field lines we integrate once backward and once
       // forward from the same starting point
      
       //backward integration
       
        for (unsigned int  i=0; i< vecEquationOfMotion.size();i++)
          { if (vecEquationOfMotion[i]) 
	           vecEquationOfMotion[i]->
		          SetBackwardDirectionOfIntegration(true);
	   
	  }
        theRunManager->BeamOn(1);
       
       //forward integration
        for (unsigned int  i=0; i< vecEquationOfMotion.size();i++)
          { if (vecEquationOfMotion[i]) 
	           vecEquationOfMotion[i]->
		          SetBackwardDirectionOfIntegration(false);
	   
	  }
        theRunManager->BeamOn(1); 
  
       }
  
   
  
  // Remove trajectory filter to PropagatorInField
    //It was for old solution when using smooth trajectory filter
       /*G4TransportationManager::GetTransportationManager()
            ->GetPropagatorInField()->SetTrajectoryFilter(0);*/

 
 	    
  
  
  //back to User defined actions and other parameters
  //--------------------------------
  
      G4TransportationManager::GetTransportationManager()
                      ->GetPropagatorInField()
	           ->SetLargestAcceptableStep(previous_largest_acceptable_step);
  
  
  
      //return to User actions
      theRunManager->SetUserAction(fUserRunAction);
      theRunManager->SetUserAction(fUserEventAction);
      theRunManager->SetUserAction(fUserPrimaryAction);
      theRunManager->SetUserAction(fUserSteppingAction);
      theRunManager->SetUserAction(fUserStackingAction);
      theRunManager->SetUserAction(fUserTrackingAction);
  
  
     //set user defined chord finders and largest accpetable step 
      for (unsigned int i=0;i<vecFieldManagers.size();i++)
        {if (user_chord_finders[i])
	   vecFieldManagers[i]->SetChordFinder(user_chord_finders[i]);
	} 
}  



////////////////////////////////////////////////////////////////
/*bool BlineTool::CheckMagneticFields()
{//Check FieldManagers
  
  if (vecFieldManagers[0] != 
        G4TransportationManager::GetTransportationManager()
	                                            ->GetFieldManager())
                                                                return false;
  if (vecMagneticFields[0] != 
     G4TransportationManager::GetTransportationManager()
                                   ->GetFieldManager()
				      ->GetDetectorField()) return false;
  G4LogicalVolumeStore* theVolumeStore 
                           = G4LogicalVolumeStore::GetInstance(); 						  
   
  std::vector<G4FieldManagers*> LogicalVolumeFields;
  unsigned int j=0;
  for (unsigned int i=0; i<theVolumeStore.size();i++)
       {if (theVolumeStore[i]->GetFieldManager())
              {j++;
	       if (j >= vecFieldManagers.size()) return false;
	       if (vecFieldManagers[j] != 
	               theVolumeStore[i]->GetFieldManager())
		                                            return false;
	       if (vecMagneticFields[j] != 
	               theVolumeStore[i]->GetFieldManager()->
				               GetDetectorField())
		                                            return false;
	       }
	}       						     						     
    
  if (j<vecFieldManagers.size()) return false;
  
  
 return true; 
   						  
 
        
}*/



////////////////////////////////////////////////////////////////
void BlineTool::ResetChordFinders()
{for (unsigned int i=0; i<vecEquationOfMotion.size();i++)
   {delete vecEquationOfMotion[i];
    delete vecChordFinders[i];
   } 
   
 vecChordFinders.clear();
 vecFieldManagers.clear();
 vecMagneticFields.clear();
 vecEquationOfMotion.clear();
 
 
 //global field
 vecChordFinders.push_back(0);
 vecMagneticFields.push_back(0);
 vecEquationOfMotion.push_back(0);
 
 vecFieldManagers.push_back
               (G4TransportationManager::GetTransportationManager()
	                                          ->GetFieldManager());
 if (vecFieldManagers[0])
     {vecMagneticFields[0]= (G4MagneticField*) 
                               vecFieldManagers[0]->GetDetectorField();
      if (vecMagneticFields[0])
         {vecEquationOfMotion[0]= new BlineEquation(vecMagneticFields[0]);
	  G4CashKarpRKF45*   pStepper  = 
	                     new G4CashKarpRKF45(vecEquationOfMotion[0]);
	  G4MagInt_Driver* pIntgrDriver 
	                       = new G4MagInt_Driver(0.1*km, pStepper, 
                                                  pStepper->GetNumberOfVariables());		     
	  vecChordFinders[0] = new G4ChordFinder(pIntgrDriver);
	 }      
      
     } 
     
  //local fields   
     
  G4LogicalVolumeStore* theVolumeStore 
                           = G4LogicalVolumeStore::GetInstance(); 						  
   
  unsigned int j=0;
  for (unsigned int i=0; i<theVolumeStore->size();i++)
       {if ((*theVolumeStore)[i]->GetFieldManager())
         {j++;
	  vecFieldManagers.push_back(((*theVolumeStore)[i])->GetFieldManager());
	  vecMagneticFields.push_back(
	                              (G4MagneticField*)
				          vecFieldManagers[j]->GetDetectorField());
	  vecEquationOfMotion.push_back(0);
	  vecChordFinders.push_back(0);
	  if (vecMagneticFields[j])
	     {vecEquationOfMotion[j]= new BlineEquation(vecMagneticFields[j]);
	      G4CashKarpRKF45*   pStepper  = 
	                     new G4CashKarpRKF45(vecEquationOfMotion[j]);
	      G4MagInt_Driver* pIntgrDriver 
	                       = new G4MagInt_Driver(0.1*km, pStepper, 
                                                  pStepper->GetNumberOfVariables());		     
	      vecChordFinders[j] = new G4ChordFinder(pIntgrDriver);
	     } 
	  }
	}             
}





