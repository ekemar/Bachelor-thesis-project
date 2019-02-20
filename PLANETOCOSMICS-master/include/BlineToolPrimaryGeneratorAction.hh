//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolPrimaryGeneratorAction.hh			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 
// DESCRIPTION
// -----------
//
// This class defines the primary generator actions used by the BlineTool
//    when tracing magnetic field lines 
// It generates primary vertex by using the user defined primary generator
//  actions and replace only the  definition of the particle  to be tracked 
//  by a Charged Geantino.
// By this way the user control the way he wants to define start position
// for field line tracking in its own primary generator action
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifndef BlineToolPrimaryGeneratorAction_h
#define BlineToolPrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "vector"



class G4ParticleGun;
class G4GeneralParticleSource;
class BlineToolPrimaryMessenger;
class G4Event;
class BlineTool;



class BlineToolPrimaryGeneratorAction : 
          public G4VUserPrimaryGeneratorAction 
{
  public:
    BlineToolPrimaryGeneratorAction();    
    ~BlineToolPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void SetUserPrimaryAction(G4VUserPrimaryGeneratorAction* anAction)
                                             {fUserPrimaryAction=anAction;};
   
    
    		      
 

  private:
   
    G4VUserPrimaryGeneratorAction* fUserPrimaryAction;
    G4bool FirstPartOfBline;
    G4ThreeVector BlineStartPosition;
    G4double T0;
};

#endif


