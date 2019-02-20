#ifndef PLANETOCOSHadronicPhysics_h
#define PLANETOCOSHadronicPhysics_h 1
/////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VPhysicsConstructor.hh"


//Builders

#include "G4PiKBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4NeutronBuilder.hh"
#include "G4MiscLHEPBuilder.hh"
#include "G4StoppingHadronBuilder.hh"

class MYBinaryProtonBuilder;
class MYBinaryNeutronBuilder;



////////////////////////////////////////////////////////////////////////////////
//
class PLANETOCOSHadronicPhysics : public G4VPhysicsConstructor
{
  public: 
    PLANETOCOSHadronicPhysics (const G4String& name ="QGSP_BIC_HP");
    virtual ~PLANETOCOSHadronicPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();
    
    inline void SetHadronicPhysicsType(G4String name)
                                 {HadronicPhysicsType = name;};
    void SetMaxAandZForFermiBreakUp(G4int anA,G4int aZ);
    void SetMinEForMultiFrag(G4double anE);
  private:
  
   G4String HadronicPhysicsType;
   void ConstructGeneral();
   
   
   void ConstructLHEPBasedPhysics(G4bool bic_flag,G4bool bert_flag,G4bool hp_flag);
   void ConstructLHEP();
   void ConstructLHEP_HP();
   void ConstructLHEP_BIC();
   void ConstructLHEP_BIC_HP();
   void ConstructLHEP_BERT();
   void ConstructLHEP_BERT_HP();
   
   void ConstructPartonStringBasedPhysics(G4bool ftf_flag, G4bool chips_flag,
   						G4bool bic_flag,G4bool bert_flag,G4bool hp_flag);
   void ConstructQGSP(); 
   void ConstructQGSP_HP();
   void ConstructQGSP_BIC();
   void ConstructQGSP_BIC_HP();
   void ConstructQGSP_BERT();
   void ConstructQGSP_BERT_HP();
   
   void ConstructQGSC(); 
   void ConstructQGSC_HP();
   void ConstructQGSC_BIC();
   void ConstructQGSC_BIC_HP();
   void ConstructQGSC_BERT();
   void ConstructQGSC_BERT_HP();
   
   void ConstructFTFP(); 
   void ConstructFTFP_HP();
   void ConstructFTFP_BIC();
   void ConstructFTFP_BIC_HP();
   void ConstructFTFP_BERT();
   void ConstructFTFP_BERT_HP();
   
   void ConstructFTFC(); 
   void ConstructFTFC_HP();
   void ConstructFTFC_BIC();
   void ConstructFTFC_BIC_HP();
   void ConstructFTFC_BERT();
   void ConstructFTFC_BERT_HP();
   
   
   void ConstructElasticPhysics();
   
   

  private:
  
   //Builders
    G4NeutronBuilder theNeutrons;
    G4PiKBuilder thePiK;
    G4ProtonBuilder thePro;
    G4MiscLHEPBuilder theMiscLHEP;
    G4StoppingHadronBuilder theStoppingHadron;
    
   //BinaryCascade Builders
    MYBinaryProtonBuilder*  theBinaryProton;
    MYBinaryNeutronBuilder* theBinaryNeutron;
    
    bool NeutronHPElasticFlag;
     

};
////////////////////////////////////////////////////////////////////////////////
#endif
