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
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TApplication.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < int >+;
#pragma link C++ class std::vector < Double_t >+;
#endif 

using namespace std;

//Section 1 reading files, setting up stuff
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

/*struct set_sim_vectors
{
	vector<int>* PDG_vec = 0;
 	vector<Double_t>* Energy_vec = 0;
	vector<Double_t>*zenith_vec = 0;
	vector<Double_t>*azimuth_vec= 0;
	vector<Double_t>*latitude_vec= 0;
	vector<Double_t>*longitude_vec = 0;
}primary[6], secondary[6]; 
*/

void set_properties()
{
	for(int n = 0; n < 6; n++)
	{
		
	}
	cout<<"\n";
}


 TH1D *filter_and_filling(TH1D *histopointer, int particle_monto_carlo_number, unsigned iterator, vector<int> *some_vector1, vector<Double_t> *some_vector2, vector<Double_t> *some_vector3)
{

	if(some_vector1->at(iterator) == particle_monto_carlo_number)
	{
		histopointer->Fill(some_vector2->at(iterator)/10000., 1./cos(some_vector3->at(iterator)));
	}
	return histopointer;

}

TH1D *filter_pointer(double latitud_mini, double latitud_maxi, double longitud_mini, double longitude_maxi, vector<Double_t>* latitude_vec, vector<Double_t>* longitude_vec,vector<int>* PDG_vector, vector<Double_t>* energy_vector,
   vector<Double_t>* zenith_vector, int particle_monto_carlo_number, TH1D *pointer1, TH1D *pointer2, TH1D *pointer3)
{

	for (unsigned index =0; index < longitude_vec->size(); index++){
		if((latitud_mini <  latitude_vec->at(index) && latitude_vec->at(index) <  latitud_maxi) &&
		 (longitud_mini < longitude_vec->at(index) && longitude_vec->at(index)< longitude_maxi)  
		&& (0.1E3 < energy_vector->at(index) && energy_vector->at(index) < 4.4E3))
		
		{

			filter_and_filling(pointer1, -particle_monto_carlo_number, index, PDG_vector, energy_vector, zenith_vector);
			filter_and_filling(pointer2, -particle_monto_carlo_number, index, PDG_vector, energy_vector, zenith_vector);
			filter_and_filling(pointer3, particle_monto_carlo_number, index, PDG_vector, energy_vector, zenith_vector);

		}
	
	}
	return 0;
}

void iterate_over_file(TTree *propertie, double latitud_mini, double latitud_maxi, double longitud_mini, double longitude_maxi, vector<Double_t>* latitude_vec, vector<Double_t>* longitude_vec,vector<int>* PDG_vector, vector<Double_t>* energy_vector, vector<Double_t>* zenith_vector, int particle_monto_carlo_number, TH1D *pointer1, TH1D *pointer2, TH1D *pointer3)
{

	for (float eventnumber =  0 ; eventnumber < propertie->GetEntriesFast();++eventnumber)
	{
		propertie->GetEntry(eventnumber);
		{
		filter_pointer(latitud_mini, latitud_maxi, longitud_mini, longitude_maxi, latitude_vec, longitude_vec,  PDG_vector, energy_vector,
		zenith_vector,particle_monto_carlo_number, pointer1, pointer2, pointer3);
		}
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

TH1D *normalized_histogram(TH1D *histogram, int bin_numbers, double norm)
{
 	double value=histogram->GetBinContent(bin_numbers+1); 
	histogram->SetBinContent(bin_numbers+1, value);
	
	histogram->Scale(norm, "width");
	return histogram;
}
//Section 2 calculating stuff for using as constants in the analysis
double calculated_norm(double solidangle, double area_unit_converstion, double integrated_primary_flux, double sum_of_entries, double area_detector)
{
	double normfactor = 0;
	normfactor = (1/solidangle)*area_unit_converstion*(integrated_primary_flux/sum_of_entries) * area_detector; 
	return normfactor;
}

double ratio_area(double startalt, double detector_altitude, double detector_latitude_max, double detector_latitude_min, double detector_longitude_max, double detector_longitude_min)
{	
	const double pi = 3.1416;
	double earth_radii = 6317;
	double geometrical_factor = 0;
	geometrical_factor = 4 * pi*pow(earth_radii+startalt,2)/(pow(earth_radii+detector_altitude,2)*(detector_longitude_max*0.017453-detector_longitude_min*0.017453)*(sin(detector_latitude_max*0.017453)-sin(detector_latitude_min*0.017453)));
	return geometrical_factor;

}



//Section 3 analysis

struct axis{Double_t axis_array[13];}xgraph, ygraph; 

void get_axis(TH1D *histogram, int bin_numbers, double norm) 
{
	Double_t an_axis[bin_numbers];	
	TH1D *the_normalized_histogram = normalized_histogram(histogram, bin_numbers, norm);
	for (int n=0; n<bin_numbers; n++)
	{
		xgraph.axis_array[n]=the_normalized_histogram->GetBinCenter(n+1); 
		ygraph.axis_array[n]=the_normalized_histogram->GetBinContent(n+1);
	}
}

Double_t *ratio_calculations(int bin_numbers, Double_t *flux_data)
{
	Double_t *ratio;
	ratio = new Double_t[bin_numbers];
	for(int n = 0; n <bin_numbers; n++)
	{
		if(0 < flux_data[n])
		{

			ratio[n] = ygraph.axis_array[n]/flux_data[n];
		}
	}
	return ratio;

}
float  calculate_ratio_between_particles(TH1D *pointer1, TH1D *pointer2)
{

		float ratio = float (pointer1->GetEntries()/pointer2->GetEntries());
		return ratio;
}

//Section 4 showing data, graphs, creating tabels 
TGraph *create_graph(TString title,int bin_numbers, Double_t *xaxis, Double_t *yaxis,Float_t color)
{
	TGraph *graph = new TGraph(bin_numbers, xaxis, yaxis);
	graph->SetTitle(title);
	graph->SetMarkerColor(color);
	graph->SetMarkerSize(0.8);
	graph->SetMarkerStyle(21);
	graph->SetLineWidth(2);
	//graph->Draw("LP");

	return graph;
}

void draw_two_graphs(TString name, TString title,TString graph_title1, TString graph_title2, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *yaxis2, Float_t color, Float_t color2) 
{
	TMultiGraph *multipointer = new TMultiGraph();
	multipointer->SetNameTitle(name, title);
	TGraph *graph_pointer = create_graph(graph_title1, bin_numbers, xaxis, yaxis, color);
	graph_pointer->Draw("LP");
	TGraph *graph_pointer2 = create_graph(graph_title2, bin_numbers, xaxis, yaxis2, color2);
	graph_pointer2->Draw("LP");
	multipointer->Add(graph_pointer);
	multipointer->Add(graph_pointer2);
	multipointer->Draw("ALP");
	multipointer->GetXaxis()->SetTitle("Kinetic energy MeV");
	multipointer->GetYaxis()->SetTitle("Antiproton flux(1/(sr s MeV m^2))" );

}

void draw_graph(TString name, TString title,TString graph_title, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Float_t color) 
{
	TMultiGraph *multipointer = new TMultiGraph();
	multipointer->SetNameTitle(name, title);
	TGraph *graph_pointer = create_graph(graph_title, bin_numbers, xaxis, yaxis, color);
	graph_pointer->Draw("P");
	multipointer->Add(graph_pointer);
	multipointer->Draw("ALP");
	multipointer->GetXaxis()->SetTitle("Kinetic energy MeV");
	multipointer->GetYaxis()->SetTitle("Ratio simulated/experimental flux");
	

}
void create_canvas(TString name, TString title, TString name2, TString title2, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *yaxis2, Double_t *yaxis3, Float_t color, Float_t color2) 
{

	TCanvas *canvas = new TCanvas("c3","plots energy",600,600,1800, 1200);
	canvas->Divide(2,1); 
	canvas->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	draw_two_graphs(name, title,"experimental flux", "simulated flux", bin_numbers, xaxis, yaxis, yaxis2, color, color2);
	gPad->BuildLegend();
	canvas->Update();
	canvas->cd(2);
	gPad->SetLogx();
	draw_graph(name2, title2,"ratio", bin_numbers, xaxis, yaxis3, color);
	gPad->BuildLegend();
	canvas->Update();
	canvas->Modified();

}


void create_and_write_file(string file_name,int bin_numbers, Double_t *indata_array1, Double_t *indata_array2)
{
	
	ofstream myfile(file_name);
	ofstream outputFile;
	
	outputFile.open(file_name);	
	
	myfile<<"Simulated particle flux;"<<"Kinetic energy"<<"\n";
	myfile<<"1/(sr s GeV m^2)"<<";"<<"GeV"<<"\n";
	

	for(int i = 0; i < bin_numbers; i++)
	{
		myfile<< indata_array1[i]<<"; "<<indata_array2[i]<<"\n";
		
	}
	outputFile.close();
}

//section 5 functions for checking stuff quick
void check_content_vector(string name, int bin_numbers, vector<Double_t> a_vector )
{

	cout<<name<<" ";
	for(int i = 0; i< bin_numbers; i++)
	 {
		cout<<a_vector.at(i)<<", ";
	}
	cout<<"\n";
}

void check_content_array(string name, int bin_numbers, Double_t *an_array )
{

	cout<<name<<" ";
	for(int i = 0; i< bin_numbers; i++)
	 {
		cout<<an_array[i]<<", ";
	}
	cout<<"\n";
}

Double_t *converted_array(int bin_numbers, Double_t *start_array)
{
	Double_t *modify_array;
	modify_array = new Double_t[bin_numbers];
	for(int i = 0; i < bin_numbers; i++)
	{
		modify_array[i] = 1000 * start_array[i];	
	}	
	return modify_array;
}
//Section 6 functions for filter stuff

int main(int argc, char* argv[]){
	TApplication theApp("App", &argc, argv);
	axis ax;
	Float_t all_Events = 0;
	Float_t	total_events = 0;
// Some definitons for stringtypes
	TString diagram_name = "Flux";
//Definding stuff with constant value
	double startaltitude=70000; //Altitude of the primary sphere
	double detector_altitude=25.2842;
	double detector_altitude2=25.2842;
	double longitud_max= 180;  //Range of laltitude and longitude of the detector
	double longitud_min= -180;
	double latitud_max= -75;
	double latitud_min= -85;
	//Some tabel values
	int nbins = 13;
	double the_integrated_primary_flux= 1.304709; 
	double area_conversation_factor = 10000;
	double the_solidangle=0.628319;
	//Define considered solid angle 
	//Definding binding range
	vector<Double_t> BESS_data_particle_flux_vector= get_data_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
 	vector<Double_t> BESS_data_energybins_vector=get_data_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
	cout<<"Bin numbers: "<<BESS_data_energybins_vector.size()<<endl;
	Double_t *BESS_data_particle_flux= BESS_data_particle_flux_vector.data();
	Double_t *BESS_data_energybins=BESS_data_energybins_vector.data();

	//TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary/";
	TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/testbins2/";
	const char * prefix = "tree";
	TString savepath = "manipulationstuff.root";
	const char *dirname;
	//const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary";
	const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/testbins2";
	const char *suffix = ".root";
	double norm_factor;
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

	if (!files) {
		cout << "!files"<<endl;
		return 1;
	}

	// Iterate over the folders with files
	TH1D *flux_histogram = new TH1D("fluxene","calculated flux", nbins, BESS_data_energybins);
	TH1D *secondary_energy = new TH1D("secondary_energy","secondary_energy",1000,0,1E3);
	TH1D *secondary_energy_antiprotons = histogram_setupp("energy_antiprotons", "Hits antiprotons", 1000, 0, 1E3);
	TH1D *secondary_energy_protons = histogram_setupp("energy_protons", "Hits protons", 1000, 5.6, 10E6);
	
	
	
	while ((file = (TSystemFile*)next())){
		filename = file->GetName();
		if(!filename.BeginsWith(prefix) || !filename.EndsWith(suffix)){
			continue;
		}
		//cout<< "Reading "<<filename<<endl;
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


		total_events += properties->GetEntriesFast();


		double detector_area = ratio_area(startaltitude, detector_altitude, latitud_max,latitud_min, longitud_max, longitud_min);
		norm_factor = calculated_norm(the_solidangle ,area_conversation_factor, the_integrated_primary_flux,total_events, detector_area);

		iterate_over_file(properties, latitud_min, latitud_max, longitud_min, longitud_max, secondary_latitude_vec, secondary_longitude_vec, secondary_PDG_vec, secondary_Energy_vec,
		 secondary_zenith_vec, 2212, flux_histogram, secondary_energy_antiprotons,secondary_energy_protons);

		
	//Clean up
	properties->ResetBranchAddresses();
	infile->Close();

	}
	
 		
				
	outfile = new TFile(savepath,"RECREATE");
	secondary_energy_protons->Write("proton second");
	flux_histogram->Write("test");
	secondary_energy_antiprotons->Write("antiproton");

	check_content_vector("BESS energy bins", nbins+1, BESS_data_energybins_vector);
	check_content_vector("BESS particle flux", nbins, BESS_data_particle_flux_vector);
	Double_t *MeV_array = converted_array(nbins-1, BESS_data_particle_flux);
	
//	set_properties();
	get_axis(flux_histogram,nbins, norm_factor);
	check_content_array("XGraph", nbins,  xgraph.axis_array);
	check_content_array("YGraph", nbins +1,  ygraph.axis_array);
	Double_t *firstratio = ratio_calculations(nbins-1, MeV_array);
	create_canvas("Flux", "flux", "Ratio", "Ratio" , nbins, xgraph.axis_array, MeV_array, ygraph.axis_array, firstratio, 4, 3);
	create_and_write_file("test.cvs", nbins, ygraph.axis_array, xgraph.axis_array);
	memset(BESS_data_energybins, 0, 14*sizeof(*BESS_data_energybins));
	memset(BESS_data_particle_flux, 0, 13*sizeof(*BESS_data_particle_flux));
	memset(BESS_data_particle_flux, 0, 13*sizeof(*MeV_array));
	memset(firstratio, 0, 13*sizeof(*firstratio));
	outfile->Write();
	outfile->Close();
	
	cout<<"run done"<<endl;
	theApp.Run();

	return 0;

}
