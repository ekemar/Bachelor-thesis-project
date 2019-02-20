#include "ROOTPartSourceGenerator.hh"
#include "Randomize.hh"

ROOTPartSourceGenerator* ROOTPartSourceGenerator::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
ROOTPartSourceGenerator::ROOTPartSourceGenerator()  
{
#ifdef USE_ANALYSIS_ROOT    
  theFirst1DHisto=0;
  theSecond1DHisto=0;
  the2DHisto=0; 
  theRandomGenerator = new TRandom();
#endif   	
}
////////////////////////////////////////////////////////////////////////////////
//  
ROOTPartSourceGenerator::~ROOTPartSourceGenerator() 
{
#ifdef USE_ANALYSIS_ROOT    
  if (theFirst1DHisto) delete theFirst1DHisto;
  if (theSecond1DHisto) delete theSecond1DHisto;
  if (the2DHisto) delete the2DHisto;
#endif 
}
////////////////////////////////////////////////////////////////////////////////
//
ROOTPartSourceGenerator* ROOTPartSourceGenerator::GetInstance()
{
  if (instance == 0) instance = new ROOTPartSourceGenerator;
  return instance;
}
#ifdef USE_ANALYSIS_ROOT 
////////////////////////////////////////////////////////////////////////////////
//
void ROOTPartSourceGenerator::ReadFirstSourceHisto(G4String file_name, G4String path)
{ TFile* root_file = new TFile(file_name,"READ");
  
  if (theFirst1DHisto) delete theFirst1DHisto;
  if (the2DHisto) delete the2DHisto;
  
  G4String aName = root_file->Get(path)->ClassName();
  G4cout<<aName<<std::endl;
  if (aName == "TH1F")  theFirst1DHisto = dynamic_cast<TH1F*> (root_file->Get(path));
  if (aName == "TH2F")  the2DHisto = dynamic_cast<TH2F*> (root_file->Get(path));
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void ROOTPartSourceGenerator::ReadSecondSourceHisto(G4String file_name, G4String path)
{ TFile* root_file = new TFile(file_name,"READ");
  
  if (theSecond1DHisto) delete theSecond1DHisto;
  
  
  G4String aName = root_file->Get(path)->ClassName();
  if (aName == "TH1F")  theSecond1DHisto = dynamic_cast<TH1F*> (root_file->Get(path));
 
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector<double> ROOTPartSourceGenerator::GeneratePrimary()
{
std::vector<double> res;
res.clear();
if (theFirst1DHisto) {
	G4double var1;
	var1 = theFirst1DHisto->GetRandom();
	res.push_back(var1); 
	if (theSecond1DHisto) {
		G4double var2;
		var2 = theSecond1DHisto->GetRandom();
		res.push_back(var2); 		
	}
}
else if (the2DHisto) {
	
/*	G4double var1,var2;
	the2DHisto->GetRandom2(var1,var2);
	res.push_back(var1);
	res.push_back(var2);
	return res;
*/	
   	double* fIntegral=0;
	fIntegral = the2DHisto->GetIntegral();
   	int nbinsx = the2DHisto->GetNbinsX();
	int nbinsy = the2DHisto->GetNbinsY();
	int nbins  = nbinsx*nbinsy;
	double fEntries = the2DHisto->GetEntries();
	double integral;
   	if (fIntegral) {
      		if (fIntegral[nbins+1] != fEntries) integral = the2DHisto->ComputeIntegral();
		fIntegral = the2DHisto->GetIntegral();
   	} 
	else {
      		integral = the2DHisto->ComputeIntegral();
		fIntegral = the2DHisto->GetIntegral();
      		if (integral == 0 || fIntegral == 0) return res;
   	}
   
   	//double r1 = theRandomGenerator->Rndm();
	double r1 = HepRandom::getTheGenerator()->flat();
   	int ibin = TMath::BinarySearch(nbins,fIntegral,r1);
   	int biny = ibin/nbinsx;
   	int	binx = ibin - nbinsx*biny;
	TAxis* fXaxis=the2DHisto->GetXaxis();
	TAxis* fYaxis=the2DHisto->GetYaxis();
   	double x = fXaxis->GetBinLowEdge(binx+1)+fXaxis->GetBinWidth(binx+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   	//double y = fYaxis->GetBinLowEdge(biny+1) + fYaxis->GetBinWidth(biny+1)*theRandomGenerator->Rndm();
	double y = fYaxis->GetBinLowEdge(biny+1) + fYaxis->GetBinWidth(biny+1)*HepRandom::getTheGenerator()->flat();
	res.push_back(x);
	res.push_back(y);

}
return res;
}
#endif 
