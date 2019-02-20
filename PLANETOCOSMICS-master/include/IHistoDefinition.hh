#ifndef IHistoDefinition_HH
#define IHistoDefinition_HH

#ifndef USE_ANALYSIS_ROOT
#include "AIDA/AIDA.h"
using namespace AIDA;
#define HISTO1D IHistogram1D
#define HISTO2D IHistogram2D
#define HISTOBASE IHistogram 
#else
//ROOT analysis
//#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TDirectory.h"
#define HISTOBASE TH1
#define HISTO1D TH1F
#define HISTO2D TH2F
#endif 
#endif
