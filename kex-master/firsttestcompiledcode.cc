#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cmath> 
#include <cstdlib>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TObject.h>
#include <TCollection.h>
#include <TString.h>
#include <TIterator.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < int >+;
#pragma link C++ class std::vector < Double_t >+;
#endif 

using namespace std;
	

int main(int argc, char* argv[]){
	
//TFile *infile = TFile::Open("/nfs/mds/project/e18sat/sim_ekemar_mds/testdata/tree1000_1472943023.root");
TFile *infile = TFile::Open("/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/scripts/scripts/testdata/treetest2.root_1495453172.root");
	//---------Read in Tree from file-------------------------
	
	TTree *properties = (TTree*) infile->Get("properties");
	if(!properties)
	{
		cout<<"file not found!"<<endl;
		return 0;
	}
	
	//create vectors for output of the tree-------------------
	
	vector<int>* primary_PDG_vec = 0;
	vector<Double_t>* primary_Energy_vec = 0;
	vector<Double_t>* primary_zenith_vec = 0;
	vector<Double_t>* primary_azimuth_vec= 0;
	vector<Double_t>* primary_latitude_vec= 0;
	vector<Double_t>* primary_longitude_vec= 0;
	
	vector<int>* secondary_PDG_vec= 0;
	vector<Double_t>* secondary_Energy_vec= 0;
	vector<Double_t>* secondary_zenith_vec= 0;
	vector<Double_t>* secondary_azimuth_vec= 0;
	vector<Double_t>* secondary_latitude_vec= 0;
	vector<Double_t>* secondary_longitude_vec= 0;
	vector<Double_t>* secondary_boundDet_vec= 0;



	
	//------Set Branch Addresses to the vectors----------------------
	
	
	properties->GetBranch("primary_PDGcode")->SetAddress(&primary_PDG_vec);
	properties->GetBranch("primary_Energy")->SetAddress(&primary_Energy_vec);
	properties->GetBranch("primary_zenith")->SetAddress(&primary_zenith_vec);
	properties->GetBranch("primary_azimuth")->SetAddress(&primary_azimuth_vec);
	properties->GetBranch("primary_latitude")->SetAddress(&primary_latitude_vec);
	properties->GetBranch("primary_longitude")->SetAddress(&primary_longitude_vec);	

	properties->GetBranch("secondary_PDGcode")->SetAddress(&secondary_PDG_vec);
	properties->GetBranch("secondary_Energy")->SetAddress(&secondary_Energy_vec);
	properties->GetBranch("secondary_zenith")->SetAddress(&secondary_zenith_vec);
	properties->GetBranch("secondary_azimuth")->SetAddress(&secondary_azimuth_vec);
	properties->GetBranch("secondary_latitude")->SetAddress(&secondary_latitude_vec);
	properties->GetBranch("secondary_longitude")->SetAddress(&secondary_longitude_vec);
	properties->GetBranch("secondary_BoundaryDetector")->SetAddress(&secondary_boundDet_vec);


	//data analysis: loop through tree via TTree->GetEntry(a) and do your analysis. Example: Plot secondary proton energy distribution
	
	TH1D *secondary_energy = new TH1D("secondary_energy","secondary_energy",1000,0,1E4);

	Double_t	energy_secondary=0;
	Int_t 		PDG_secondary=0;
	Int_t		antiproton_counter = 0;
	Float_t		total_particles = 0;
	Float_t        	quota_particles_events = 0;
	Float_t		eventnumber = 0;
	Float_t		percent_of_particles_event = 0;
	Float_t		total_Events;


	for (float eventnumber =  0 ; eventnumber < properties->GetEntriesFast();++eventnumber){
		properties->GetEntry(eventnumber);
		antiproton_counter = 0;
		
		for(auto pdg_index = secondary_PDG_vec->begin(); pdg_index != secondary_PDG_vec->end();++pdg_index)
		{
			if(*pdg_index == -2212)
			{
				secondary_energy->Fill(secondary_Energy_vec->at(pdg_index - secondary_PDG_vec->begin()));
					antiproton_counter++;
					//multi_antiproton->Fill(antiproton_counter);
				
			}
		}
		total_particles += antiproton_counter;	
		total_Events = eventnumber;

		if(antiproton_counter > 0){ 
			cout<< " No of antiprotons: " << antiproton_counter <<" eventnumber: " << eventnumber << endl;
		}
			//return eventnumber;
	}

	
	//Write histogram to file
	TFile* outfile = new TFile("output.root","RECREATE");
	secondary_energy->Write();
	outfile->Write();
	outfile->Close();
	
	//Some calulations

	quota_particles_events = float (total_particles / total_Events);
	percent_of_particles_event = quota_particles_events * 100;

	// printing stuff
	cout<< " Total no of antiprotons: " << total_particles << endl;
	cout<< "Eventnumbers : " << total_Events << endl;
	cout<< "Quota : " << quota_particles_events<< endl;
	cout<< "Percent of events who gives hits: " << percent_of_particles_event<< "%"<< endl;

	//Clean up
	properties->ResetBranchAddresses();
	infile->Close();
	return 0;
	
}
