{
	#include <vector>

	TFile *infile = TFile::Open("/nfs/hicran/project/e18sat/SIM/new_planetocosmics/planetocosmics/planetocosmics/examples/earth_testcomplete5_1469690213.root");
	
	//---------Read in Tree from file-------------------------
	
	TTree *properties = (TTree*) infile.Get("properties");
	
	//create vectors for output of the tree-------------------
	
	vector<int>* primary_PDG_vec;
	vector<Double_t>* primary_Energy_vec;
	vector<Double_t>* primary_zenith_vec;
	vector<Double_t>* primary_azimuth_vec;
	vector<Double_t>* primary_posY_vec;
	vector<Double_t>* primary_posX_vec;
	
	vector<int>* secondary_PDG_vec;
	vector<Double_t>* secondary_Energy_vec;
	vector<Double_t>* secondary_zenith_vec;
	vector<Double_t>* secondary_azimuth_vec;
	vector<Double_t>* secondary_posY_vec;
	vector<Double_t>* secondary_posX_vec;
	vector<Double_t>* secondary_boundDet_vec;



	
	//------Set Branch Addresses to the vectors----------------------
	
	
	properties->SetBranchAddress("primary_PDGcode",&primary_PDG_vec);
	properties->SetBranchAddress("primary_Energy",&primary_Energy_vec);
	properties->SetBranchAddress("primary_zenith",&primary_zenith_vec);
	properties->SetBranchAddress("primary_azimuth",&primary_azimuth_vec);
	properties->SetBranchAddress("primary_posY",&primary_posY_vec);
	properties->SetBranchAddress("primary_posX",&primary_posX_vec);
	
	properties->SetBranchAddress("secondary_PDGcode",	&secondary_PDG_vec);
	properties->SetBranchAddress("secondary_Energy",	&secondary_Energy_vec);
	properties->SetBranchAddress("secondary_zenith",	&secondary_zenith_vec);
	properties->SetBranchAddress("secondary_azimuth",	&secondary_azimuth_vec);
	properties->SetBranchAddress("secondary_posY",		&secondary_posY_vec);
	properties->SetBranchAddress("secondary_posX",		&secondary_posX_vec);
	properties->SetBranchAddress("secondary_BoundaryDetector", &secondary_boundDet_vec);
	
	
	
	
	//data analysis: loop through tree via TTree->GetEntry(a) and do your analysis. Example: Muon creation for different primary energies
	
	TH1D *multi_mu = new TH1D("multi_mu","multiplicity of muons",100,0,1E6);
	
	Double_t	energy_primary;
	Int_t 		PDG_secondary;
	
	
	for (int a = 0 ; a < properties->GetEntriesFast();++a){
		Int_t muon_counter=0;
		properties->GetEntry(a);
		energy_primary = primary_Energy_vec.at(0);
		
		for (int b = 0; b < secondary_PDG_vec.size();++b){
			if(secondary_PDG_vec.at(b) == 13 || secondary_PDG_vec.at(b) == -13){
				++muon_counter;
			} 
			
		}
		
		if(muon_counter > 0){
			multi_mu->Fill(energy_primary,muon_counter);
		}
		
		
		cout<<"EventNO: "<<a<<" no of muons: "<<muon_counter<<" energy of prim: "<<energy_primary<<endl;
	}
	
	multi_mu->Draw();
	
}
