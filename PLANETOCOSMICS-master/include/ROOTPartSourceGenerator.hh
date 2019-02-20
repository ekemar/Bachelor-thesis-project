#ifndef ROOTPartSourceGenerator_HH
#define ROOTPartSourceGenerator_HH
#include"globals.hh"
#include"G4ios.hh"
#ifdef USE_ANALYSIS_ROOT

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TRandom.h"
//#include "TAxis.hh"
#include "TMath.h"

#include<vector>

#endif 

namespace CLHEP {}
using namespace CLHEP;

class ROOTPartSourceGenerator
{
public:

  
  ~ROOTPartSourceGenerator();
   static ROOTPartSourceGenerator* GetInstance();
   
   
  //Public method
  //--------------
#ifdef USE_ANALYSIS_ROOT    
   void ReadFirstSourceHisto(G4String file_name, G4String path);
   void ReadSecondSourceHisto(G4String file_name, G4String path);
   std::vector<double> GeneratePrimary();
#endif  
   
   
   
 

private:
  static ROOTPartSourceGenerator* instance;
 

private:
  ROOTPartSourceGenerator(); 
#ifdef USE_ANALYSIS_ROOT   
  TH1F* theFirst1DHisto;
  TH1F* theSecond1DHisto;
  
  TH2F* the2DHisto; 
  TRandom* theRandomGenerator;
#endif  
 

    
 

};


#endif




