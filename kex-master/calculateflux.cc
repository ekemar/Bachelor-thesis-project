//Liselotts code for calculating the flux from antiprotons
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
#include <TApplication.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < int >+;
#pragma link C++ class std::vector < Double_t >+;
#endif 

using namespace std;
	

int main(int argc, char* argv[]){
	TApplication theApp("App", &argc, argv);
	//Set up flux diagram	

//creating canvas 
	TCanvas *histogram = new TCanvas("histogram","antiproton hits",1000, 7000);	
	
//definitions for readin files from folder

	TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/testdata/";
	const char * prefix = "tree";
	TString savepath = "1histo7.root";
	const char *dirname;
	const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/testdata/";
	const char *suffix = ".root";
	TFile* outfile;
	int file_counter = 0;
	TTree *properties;
//Getting a list from folder	

	TSystemDirectory dir(directory_name, directory_name);
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TIter next(files);
	TString filename;
	TString full_path;
	TFile *infile;


	if (!files) {
		cout << "!files"<<endl;
		return 1;
	}

	// Iterate over the folders with files
	TH1D *secondary_energy = new TH1D("secondary_energy","secondary_energy",1000,0,7E3);
	while ((file = (TSystemFile*)next())){
		filename = file->GetName();
		if(!filename.BeginsWith(prefix) || !filename.EndsWith(suffix)){
			continue;
		}
		cout<< "Reading "<<filename<<endl;
		full_path = path+ filename;
	
	//---------Read in Tree from file-------------------------
		infile = TFile::Open(full_path);
		properties = (TTree*) infile->Get("properties");
		if(!properties)
		{
			cout<<"file not found!"<<endl;
			return 1;
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

		Double_t	energy_secondary=0;
		Int_t 		PDG_secondary=0;
		Float_t		total_particles = 0;
		Float_t        	quota_particles_events = 0;
		Float_t		eventnumber = 0;


		for (float eventnumber =  0 ; eventnumber < properties->GetEntriesFast();++eventnumber){
			properties->GetEntry(eventnumber);
			
			for(auto pdg_index = secondary_PDG_vec->begin(); pdg_index != secondary_PDG_vec->end();++pdg_index)
			{
				if(*pdg_index == -2212)
				{
					secondary_energy->Fill(secondary_Energy_vec->at(pdg_index - secondary_PDG_vec->begin()));
					
				}
			}
		}

		



	}
	
		//Create 1 histogram with all protons
		outfile = new TFile(savepath,"UPDATE");
		secondary_energy->Write();
		//histogram->Draw();
		secondary_energy->Draw();
		histogram->Update();
		theApp.Run();
		outfile->Write();
		outfile->Close();


		//Clean up
		properties->ResetBranchAddresses();
		infile->Close();
		return 0;

}
