//////////////////////////////////////////////////////////////////////////////////////// 
///		Module: 	BlineToolPrimaryGeneratorAction.cc			     ///
///		Author: 	Laurent Desorgher				     /// 
///		Version: 	1.0						     /// 
///		Last Date:	2003-10-06 08:23:30                                  ///
//////////////////////////////////////////////////////////////////////////////////////// 


#include "BlineToolPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "G4ChargedGeantino.hh"

BlineToolPrimaryGeneratorAction::BlineToolPrimaryGeneratorAction()
{ fUserPrimaryAction = 0;
  FirstPartOfBline =true;
}

BlineToolPrimaryGeneratorAction::~BlineToolPrimaryGeneratorAction()
{ 
}
/////////////////////////////////////////////////////////////////////////////
void BlineToolPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (!fUserPrimaryAction)
      {G4cout<<"The user primary action is not defined it will not work"<<std::endl;
       return;
      }



 
 //for the first part of a bline the start position and time are defined
 // by using the USER primary action while for the second part the precdent
 // values are taken 
  
  if (FirstPartOfBline)
     {//set the position and time defined by using the USER primary action 
      G4Event* tmpEvent =new G4Event();    
      fUserPrimaryAction->GeneratePrimaries(tmpEvent);
      BlineStartPosition = tmpEvent->GetPrimaryVertex()->GetPosition();
      T0 = tmpEvent->GetPrimaryVertex()->GetT0();
      delete tmpEvent;
     }
  FirstPartOfBline =!FirstPartOfBline;   
  
  
  G4PrimaryVertex* primary_vertex = 
                              new G4PrimaryVertex(BlineStartPosition, T0);
	                         
				 
 
   				 
  
 //define the particle to be tracked as Charged Geantino
    
   G4ChargedGeantino* particle_definition =G4ChargedGeantino::ChargedGeantino();
   
   
   G4double mass =  particle_definition->GetPDGMass();
   G4double energy = 10000.*MeV + mass;
   G4double pmom = sqrt(energy*energy-mass*mass);
   // the momentum direction and energy do not have an effect in tracing of 
   // bline but stilll need to be defined 
   
   G4double px = 0.;
   G4double py = 0.;
   G4double pz = pmom;

   G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
   particle->SetMass( mass );
   particle->SetCharge(particle_definition->GetPDGCharge());
   primary_vertex->SetPrimary( particle );
   
   anEvent->AddPrimaryVertex( primary_vertex );

}




 


