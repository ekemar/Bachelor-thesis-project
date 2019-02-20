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
 vector<Double_t> get_data_from_file(int interesting_col_number,string filename, int rows, int columns)
 {
      string line;
      int row,col;
      std::vector<Double_t> particle_flux;
      std::vector<Double_t> lower_energy;
      std::vector<Double_t> higher_energy;
      std::vector<Double_t> all_energies;

      float my_array[rows][columns];
      ifstream pFile (filename);
      if (pFile.is_open())
      {
          row=0;
          while(!pFile.eof())
          {
              getline(pFile, line);
              stringstream ss(line);
              col=0;
              while(ss >> my_array[row][col])
              {
                  //cout<<"col "<<col;
                 switch(interesting_col_number)
                 {
                  case 0:
                          if(col ==0)
                          {
                          particle_flux.push_back(my_array[row][0]);

                          }
                 case 1:
                         switch(col)
                         {
                                  case  1: lower_energy.push_back(my_array[row][1]);
                                  break;
                                  case  2:
                                  higher_energy.push_back(my_array[row][2]);
                                  break;
                         }
                  col++;
              }
         }
              row++;
          }
           pFile.close();
         }
         switch(interesting_col_number)
         {
                 case 0:
                 return particle_flux;
                 break;
                 case 1:
                 all_energies = lower_energy;
                 all_energies.push_back(higher_energy.back());
                 return all_energies;


         }

}


TH1D *histogram_setupp(TString name, TString title, int number_bins, int start_value, int end_value)
{
	TH1D *pointer = new TH1D(name, title, number_bins, start_value, end_value);
	return pointer;
}
TH1D *detector_check(int max_value, int min_value, int iterator,TH1D *pointer, vector<Double_t> *detector_vector) 
{

	if(min_value < detector_vector->at(iterator) && detector_vector->at(iterator) <max_value)
	{
		pointer->Fill(detector_vector->at(iterator));
	}
	return pointer;
}

double calculated_norm(double solidangle, double area_unit_converstion, double integrated_primary_flux, double sum_of_entries, double geometrical_factor)
{
	double normfactor = 0;
	normfactor = (1/solidangle)*area_unit_converstion*(integrated_primary_flux/sum_of_entries) * geometrical_factor; 
	return normfactor;
}

double calculate_geometrical_factor(double startalt, double detector_altitude, double detector_latitude_max, double detector_latitude_min, double detector_longitude_max, double detector_longitude_min)
{	
	const double pi = 3.1416;
	double earth_radii = 6317;
	double geometrical_factor = 0;
	geometrical_factor = 4 * pi*pow(6371+startalt,2)/(pow(6371+detector_altitude,2)*(detector_longitude_max*0.017453-detector_longitude_min*0.017453)*(sin(detector_latitude_max*0.017453)-sin(detector_latitude_min*0.017453)));
	return geometrical_factor;

}


TCanvas *create_canvas( TH1D *TH1D_pointer, const char *name, const char *title, Int_t width, Int_t height)
{
	TCanvas *pointer = new TCanvas(name, title, width, height);
	//pointer->Divide(numbers);
	TH1D_pointer->Draw();
	pointer->Update();
	return pointer;
}

TGraph  *create_graph(Int_t bin_numbers, vector<Double_t> some_vector)
{
	Double_t *some_array = some_vector.data();
	TGraph *pointer = new TGraph(bin_numbers, x-axis, yaxis);
	return pointer;
}

Double_t *
set_axis_stuff(double amount_bins)
{
	double yaxis = TH1D_diagram
	
}
float  calculate_ratios(TH1D *pointer1, TH1D *pointer2)
{

		cout<< "Entries antiprotons" << pointer1->GetEntries()<<endl;
		cout<<"Entries protons" <<pointer2->GetEntries()<<endl;
		float ratio = float (pointer1->GetEntries()/pointer2->GetEntries());
		return ratio;
}
int main(int argc, char* argv[]){
	TApplication theApp("App", &argc, argv);
	int total_events = 0;
	Float_t all_Events = 0;
	Float_t	total_Events = 0;
// Some definitons for stringtypes
	TString diagram_name = "Flux";
//Definding stuff with constant value
	double startalt=70000; //Altitude of the primary sphere
	double detector_altitude1=25.2842;
	double detector_altitude2=25.2842;
	double longitud_max= 180;  //Range of laltitude and longitude of the detector
	double longitud_min= -180;
	double latitud_max= -75;
	double latitud_min= -85;
	int numbers_bins = 20;
	int histogram_window = 1;
	int flux_window = 2;
	float proton_counter = 0;
	float antiproton_counter = 0;
	//Analytical factors
	
	//Some tabel values
	double the_integrated_primary_flux= 0.504709; 
	double area_conversation_factor = 1000;
	double the_solidangle=0.628319;
	//Define considered solid angle 
	//Definding binding range
	vector<Double_t> BESS2007_particle_flux= get_data_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6);
 	vector<Double_t> BESS2007_energybins=get_data_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6);
	//Set up flux diagram	

//creating canvas 
	//TCanvas *graph_windows = new TCanvas("histogram","antiproton hits",1000, 10000);	
	//graph_windows->Divide(histogram_window, flux_window);	
//definitions for readin files from folder

	//TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary/";
	TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/1moretestratio/";
	const char * prefix = "tree";
	TString savepath = "1histo9.root";
	const char *dirname;
	//const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary";
	const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/1moretestratio";
	const char *suffix = ".root";
	TFile* outfile;
	TTree *properties;
//Getting a list from folder	
		
	TSystemDirectory dir(directory_name, directory_name);
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TIter next(files);
	TString filename;
	TString full_path;
	TFile *infile;
	cout<<"yay we can read in files"<<endl;

	if (!files) {
		cout << "!files"<<endl;
		return 1;
	}

	// Iterate over the folders with files

	/*TH1D *secondary_energy = new TH1D("secondary_energy","secondary_energy",1000,0,7E3);
	//TH1D *secondary_energy_antiprotons = histogram_setupp("energy_antiprotons", "Hits antiprotons", 1000, 0, 7E3);
	TH1D *secondary_energy_protons = histogram_setupp("energy_protons", "Hits protons", 1000, 0, 7E3);
	//TH1D *secondary_energy_protons = new TH1D("secondary_energy_protons","secondary_energy_protons",100000,0,7E3);
	//TH1D *secondary_energy_antiprotons = new TH1D("secondary_energy_antiprotons","secondary_energy_antiprotons",1000,0,7E3);
	TH1D *primary_energy_protons = histogram_setupp("secondary_energy_antiprotons","secondary_energy_antiprotons",1000,0,7E3);
	TH1D *boundary_detector = histogram_setupp("boundary_detector","boundary_detector ",1000,0,100);
	TH1D *boundary_detector100km =histogram_setupp("boundary_detector 100km altitude","boundary_detector 100km ",1000,0,100);
	TH1D *boundary_detector43 = histogram_setupp("boundary_detectora 4.3 depth","boundary_detector 4.3 ",1000,0,100);
	TH1D *boundary_detector104 = histogram_setupp("boundary_detector 10.4 depth","boundary_detecto 10.4 ",1000,0,100);
	TH1D *boundary_detector156 = histogram_setupp("boundary_detector 15.6 depth","boundary_detecto 15.6 ",1000,0,100);
	TH1D *boundary_detector226 = histogram_setupp("boundary_detector 22.6 depth","boundary_detecto 22.6 ",1000,0,100);
	TH1D *boundary_detector50 = histogram_setupp("boundary_detector 50.0 depth","boundary_detecto 50.0 ",1000,0,100);
	TH1D *boundary_detector100 = histogram_setupp("boundary_detector 100.0 depth","boundary_detecto 100.0 ",1000,0,100);
	//TH1F *fluxene = new TH1F("fluxene",diagram_name,numbers_bins, BESS2007_binrange_low);
	*/
	total_events += properties->GetEntriesFast();
/*	anropa normfactor
	anropa canvaskapande	
	fyll histogram
	inormalisera histogramet
	
	fyll  normalisera histogram med experimentdata
	fyll normalisera histogram med simulerad data
	berakna ration mellan simulerad data ochexperiment datan
	skapa graf med simulerat data
	skapa grav med experimentl data
	
*/
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

		for (float eventnumber =  0 ; eventnumber < properties->GetEntriesFast();++eventnumber){
			properties->GetEntry(eventnumber);
			{

			
					
					
					int vector_index = 0;	
					for (int index =0; index < secondary_longitude_vec->size(); index++){
								
						//FILTERS===================================================================================================================================
					//detector area filter
					if( (latitud_min <  secondary_latitude_vec->at(index) && secondary_latitude_vec->at(index) <  latitud_max) &&
					 (longitud_min < secondary_longitude_vec->at(index) && secondary_longitude_vec->at(index))){
						if(secondary_PDG_vec->at(index) == -2212)
						{
							secondary_energy_antiprotons->Fill(secondary_Energy_vec->at(index));
							detector_check(110, 90, index, boundary_detector100km, secondary_boundDet_vec); 
							detector_check(45, 35, index, boundary_detector43, secondary_boundDet_vec); 
							detector_check(34, 30, index, boundary_detector104, secondary_boundDet_vec); 
							detector_check(30, 27, index, boundary_detector156, secondary_boundDet_vec); 
							detector_check(26, 22, index, boundary_detector226, secondary_boundDet_vec); 
							detector_check(21, 18, index, boundary_detector50, secondary_boundDet_vec); 
							detector_check(17, 0, index, boundary_detector50, secondary_boundDet_vec); 
							
						}

							
						if(secondary_PDG_vec->at(index) == 2212)
						{
							secondary_energy_protons->Fill(secondary_Energy_vec->at(index));
							
						}
					}

					}
					
				}
			
		}

		
	//Clean up
	properties->ResetBranchAddresses();
	infile->Close();

	}
	
		//Create 1 histogram with all protons
		outfile = new TFile(savepath,"RECREATE");
		//secondary_energy->Write();
		secondary_energy_protons->Write("proton second");
		primary_energy_protons->Write("proton prim");
		create_canvas(secondary_energy_protons, "secondary protons", "secondary protons", 1000, 700);
	//	graph_windows->cd(histogram_window);
		//secondary_energy_protons->Draw();
		//secondary_energy->Draw();
		//graph_windows->Update();
		secondary_energy_antiprotons->Write("antiproton");
		//graph_windows ->cd(flux_window);
		
//		secondary_energy_antiprotons->Draw();
	//	graph_windows->Update();
		float ratio_antiproton_proton = calculate_ratios(secondary_energy_antiprotons, secondary_energy_protons);
		cout<< "Ratio antiprotons-protons "<<ratio_antiproton_proton<<endl;
		cout<<"Total events " <<total_events<<endl;
		cout<<"Total secondary antiprotons: "<<secondary_energy_antiprotons->GetEntries()<<endl;
		float ratio_antiprotons_event = float(antiproton_counter/total_events);
		cout<<"ratio for antiprotons/event "<<ratio_antiprotons_event<<endl;
/*		cout<<"entries primary protons " <<primary_energy_protons->GetEntries()<<endl;
		cout<<"entries detector :   " <<boundary_detector->GetEntries()<<endl;
                cout<<"entries detector 100km:   " <<boundary_detector100km->GetEntries()<<endl;
                cout<<"entries detector 4.3:  " <<boundary_detector43->GetEntries()<<endl;
                cout<<"entries detector 10.4:  " <<boundary_detector104->GetEntries()<<endl;
                cout<<"entries detector 15.6: " <<boundary_detector156->GetEntries()<<endl;
                cout<<"entries detector 26.2: " <<boundary_detector226->GetEntries()<<endl;
                cout<<"entries detector 50.0:   " <<boundary_detector50->GetEntries()<<endl;
                cout<<"entries detector 100:  " <<boundary_detector100->GetEntries()<<endl;
*/
		outfile->Write();
		outfile->Close();
		
		cout<<"Hej"<<endl;
		theApp.Run();

		return 0;

}
