#include "PLANETOCOSSteppingActionMessenger.hh"
#include "PLANETOCOSSteppingAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4RunManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSteppingActionMessenger::PLANETOCOSSteppingActionMessenger(PLANETOCOSSteppingAction * msa)
:myAction(msa)
{ G4String cmd_name;
  G4String guidance;

  //directories
  //----------
 
  myStopBoundaryDir = new G4UIdirectory("/PLANETOCOS/STOPCONDITION/"); 
  myStopBoundaryDir->SetGuidance("Definition of condition for stopping particles");
  
  //commands
  //----------
  
  //Commands for stopping particle below user defined energy
  cmd_name = "/PLANETOCOS/STOPCONDITION/SetStoppingEnergy";
  AddUntrackedParticleCmd = new G4UIcommand(cmd_name,this);
  guidance ="DEfine the stopping energy of a selected type of particle.";
  AddUntrackedParticleCmd->SetGuidance(guidance); 
  AddUntrackedParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G4UIparameter* param;
  param = new G4UIparameter("Pname",'s',false);
  param->SetDefaultValue("e-");
  AddUntrackedParticleCmd->SetParameter(param);
  
  param = new G4UIparameter("E",'d',false);
  param->SetDefaultValue("0.0");
  AddUntrackedParticleCmd->SetParameter(param);
  
  param = new G4UIparameter("EUnit",'s',true);
  param->SetDefaultValue("MeV");
  AddUntrackedParticleCmd->SetParameter(param);
  
  
  
  cmd_name = "/PLANETOCOS/STOPCONDITION/DesactivateStopEnergyCondition";
  RemoveUntrackedParticleCmd = new G4UIcmdWithAString(cmd_name,this);
  guidance ="Desactivate the stopping energy condition for the slected particle. ";
  RemoveUntrackedParticleCmd->SetGuidance(guidance);
  RemoveUntrackedParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  //Commands for controlling stopping at user defined boundary
 /* cmd_name="/PLANETOCOS/STOPBOUNDARY/StopUpFluxAtSelectedBoundary";
  SetStopUpFluxAtSelectedBoundaryCmd = new G4UIcmdWithABool(cmd_name,this);
  guidance ="If true(false) particles moving upward  are (not) stopped at a selected boundary";  
  SetStopUpFluxAtSelectedBoundaryCmd->SetGuidance(guidance);
  SetStopUpFluxAtSelectedBoundaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name="/PLANETOCOS/STOPBOUNDARY/StopDownFluxAtSelectedBoundary";
  SetStopDownFluxAtSelectedBoundaryCmd = new G4UIcmdWithABool(cmd_name,this);
  guidance ="If true(false) particles moving downward  are (not) stopped at a selected boundary";  
  SetStopDownFluxAtSelectedBoundaryCmd->SetGuidance(guidance);
  SetStopDownFluxAtSelectedBoundaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name = "/PLANETOCOS/STOPBOUNDARY/SetNameOfStopBoundaryForUpFlux";
  SetNameOfStopBoundaryForUpFluxCmd = new G4UIcmdWithAString(cmd_name,this);
  guidance ="Defined the name of the volume at the bottom boundary of which";
  guidance +="upward moving particle should be stopped";
  SetNameOfStopBoundaryForUpFluxCmd->SetGuidance(guidance);
  SetNameOfStopBoundaryForUpFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/PLANETOCOS/STOPBOUNDARY/SetNameOfStopBoundaryForDownFlux";
  SetNameOfStopBoundaryForDownFluxCmd = new G4UIcmdWithAString(cmd_name,this);
  guidance ="Defined the name of the volume at the bottom boundary of which";
  guidance +="downward moving particle should be stopped";
  SetNameOfStopBoundaryForDownFluxCmd->SetGuidance(guidance);
  SetNameOfStopBoundaryForDownFluxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  */
  
  //Commands for controlling stopping at magnetopause
  cmd_name="/PLANETOCOS/STOPCONDITION/StopAtMagnetopause";
  SetStopAtMagnetopauseCmd = new G4UIcmdWithABool(cmd_name,this);
  guidance ="If true(false) the particle is (not) stopped outside the magnetopause";  
  SetStopAtMagnetopauseCmd->SetGuidance(guidance);
  SetStopAtMagnetopauseCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  				
  cmd_name="/PLANETOCOS/STOPCONDITION/SetMagnetopauseScaling";
  SetMagnetopauseOutFactorCmd = new  G4UIcmdWithADouble(cmd_name,this);
  guidance="Define the scalling factor of the magnetopause (by default 1)";
  SetMagnetopauseOutFactorCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  //Commands for stopping particle after some turns around the planet 
    
  cmd_name="/PLANETOCOS/STOPCONDITION/SetMaxNumberOfTurnAroundThePlanet";
  SetMaxNumberOfTurnAroundThePlanetCmd = new  G4UIcmdWithADouble(cmd_name,this);
  guidance="Set the maximum number of revolution that a particle can make around the planet roation axis. ";
  SetMaxNumberOfTurnAroundThePlanetCmd->SetGuidance(guidance);
  SetMaxNumberOfTurnAroundThePlanetCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  //Commands for controlling stopping at user defined altitude
  
  cmd_name="/PLANETOCOS/STOPCONDITION/SetStopAltitudeForUpwardFlux";
  SetStopAltitudeForUpwardCmd = new  G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetStopAltitudeForUpwardCmd->SetUnitCategory("Length");
  guidance="Define the altitude at which upward flux should be stopped";
  SetStopAltitudeForUpwardCmd->SetGuidance(guidance);
  SetStopAltitudeForUpwardCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name="/PLANETOCOS/STOPCONDITION/SetStopAltitudeForDownwardFlux";
  SetStopAltitudeForDownwardCmd = new  G4UIcmdWithADoubleAndUnit(cmd_name,this);
  SetStopAltitudeForDownwardCmd->SetUnitCategory("Length");
  guidance="Define the altitude for which Downward flux should be stopped";
  SetStopAltitudeForDownwardCmd->SetGuidance(guidance);
  SetStopAltitudeForDownwardCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  cmd_name="/PLANETOCOS/STOPCONDITION/StopAlsoUpwardPrimary";
  SetStopUpwardPrimaryCmd = new G4UIcmdWithABool(cmd_name,this);
  guidance="If true the primaries moving upward are also stopped when they reach a user defined altitude ";
  SetStopUpwardPrimaryCmd->SetGuidance(guidance);
  SetStopUpwardPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);   
  
  cmd_name="/PLANETOCOS/STOPCONDITION/StopAlsoDownwardPrimary";
  SetStopDownwardPrimaryCmd = new G4UIcmdWithABool(cmd_name,this);
  guidance="If true the primaries moving downward are also stopped when they reach a user defined altitude";
  SetStopDownwardPrimaryCmd->SetGuidance(guidance);
  SetStopDownwardPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSteppingActionMessenger::~PLANETOCOSSteppingActionMessenger()
{
    delete AddUntrackedParticleCmd;
    delete RemoveUntrackedParticleCmd;
     
    
  /*  delete SetStopUpFluxAtSelectedBoundaryCmd;
    delete SetStopDownFluxAtSelectedBoundaryCmd;
    delete SetNameOfStopBoundaryForUpFluxCmd;
    delete SetNameOfStopBoundaryForDownFluxCmd;*/
    
    
    delete SetStopAtMagnetopauseCmd;
    delete SetMagnetopauseOutFactorCmd; 
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSteppingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ //Commands for stopping particle below user defined energy
  if( command == AddUntrackedParticleCmd ) {
 	const char* paramString=newValue;
        G4double Ekin;
        char pname[30];
	char unts[30];
        std::istringstream is((char*)paramString);
        is >> pname >> Ekin  >> unts;
        G4String particleName=pname;
	G4String Eunit = unts;
	G4double uv = AddUntrackedParticleCmd->ValueOf(Eunit);
	Ekin=Ekin*uv;
        myAction->AddUntrackedParticle(particleName,Ekin);
  }
  else if( command == RemoveUntrackedParticleCmd ) 
        myAction->RemoveUntrackedParticle(newValue);
  
  //Commands for controlling stopping at user defined boundary
 /* else if( command ==  SetStopUpFluxAtSelectedBoundaryCmd ){
  	G4bool aBool = SetStopUpFluxAtSelectedBoundaryCmd->GetNewBoolValue(newValue);
   	myAction->SetStopUpFluxAtSelectedBoundary(aBool);
  }	
  else if( command ==  SetStopDownFluxAtSelectedBoundaryCmd ){
  	G4bool aBool = SetStopDownFluxAtSelectedBoundaryCmd->GetNewBoolValue(newValue);
   	myAction->SetStopDownFluxAtSelectedBoundary(aBool);
  }
  else if( command == SetNameOfStopBoundaryForUpFluxCmd ) 
        myAction->SetNameOfStopBoundaryForUpFlux(newValue);	
  else if( command == SetNameOfStopBoundaryForDownFluxCmd ) 
        myAction->SetNameOfStopBoundaryForDownFlux(newValue);	*/	
  
  //Commands for controlling stopping at magnetopause
  else if( command ==  SetStopAtMagnetopauseCmd ){
  	G4bool aBool = SetStopAtMagnetopauseCmd->GetNewBoolValue(newValue);
   	myAction->SetStopAtMagnetopause(aBool);
  }
  else if( command ==  SetMagnetopauseOutFactorCmd ){
  	G4double aVal = SetMagnetopauseOutFactorCmd->GetNewDoubleValue(newValue);
   	myAction->SetMagnetopauseOutFactor(aVal);
  }
   //Commands for stopping particle after some turns around the planet
  else if( command ==  SetMaxNumberOfTurnAroundThePlanetCmd ){
  	G4String geometry_type=
   	 	dynamic_cast< const PLANETOCOSGeometryConstruction* >
   			  (G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetGeometryType();
  	if (geometry_type !="SPHERICAL") {
		G4cout<<"This command can be used only in the case of spherical geometry"<<std::endl;
		return;
	}
	G4double aVal = SetMaxNumberOfTurnAroundThePlanetCmd->GetNewDoubleValue(newValue);
   	myAction->SetMaxNumberOfTurnAroundThePlanet(aVal);
  }
  else if( command ==  SetStopAltitudeForUpwardCmd ){
  	G4double aVal = SetStopAltitudeForUpwardCmd->GetNewDoubleValue(newValue);
   	myAction->SetStopAltitudeForUpward(aVal);
  }
  else if( command ==  SetStopAltitudeForDownwardCmd ){
  	G4double aVal = SetStopAltitudeForDownwardCmd->GetNewDoubleValue(newValue);
   	myAction->SetStopAltitudeForDownward(aVal);
  }
  else if( command ==  SetStopUpwardPrimaryCmd ){
  	G4bool aBool = SetStopUpwardPrimaryCmd->GetNewBoolValue(newValue);
   	myAction->SetStopUpwardPrimary(aBool);
  }
  else if( command ==  SetStopDownwardPrimaryCmd ){
  	G4bool aBool = SetStopDownwardPrimaryCmd->GetNewBoolValue(newValue);
   	myAction->SetStopDownwardPrimary(aBool);
  }
  
  		    
}



