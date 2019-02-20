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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TApplication.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < int >+;
#pragma link C++ class std::vector < Double_t >+;
#endif 

using namespace std;

//Section 1 reading files, setting up stuff
 vector<Double_t> get_one_vector_from_file(int interesting_col_number,string filename, int rows, int columns)
 {
      string line;
      int row,col;
      std::vector<Double_t> particle_flux;
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
                          if(col == interesting_col_number)
                          {
                          particle_flux.push_back(my_array[row][interesting_col_number]);

                          }
                  col++;
              }
         }
              row++;
          }
           pFile.close();
         
	return particle_flux;
}

vector<Double_t> modifiy_vector(vector<Double_t> vector1, vector<Double_t> vector2)
{
	vector<Double_t> modified_vector;	
	modified_vector = vector1;
	modified_vector.push_back(vector2.back());
	return modified_vector;


}

 TH1D *filter_and_filling(TH1D *histopointer, int particle_monto_carlo_number, unsigned iterator, vector<int> *some_vector1, vector<Double_t> *some_vector2, vector<Double_t> *some_vector3)
{
	if(some_vector1->at(iterator) == particle_monto_carlo_number)
	{
		histopointer->Fill(some_vector2->at(iterator)/1000., 1./cos(some_vector3->at(iterator)));
	}
	return histopointer;

}

TH1D *filter_pointer(double latitud_mini, double latitud_maxi, double longitud_mini, double longitude_maxi, vector<Double_t>* latitude_vec, vector<Double_t>* longitude_vec,vector<int>* PDG_vector, vector<Double_t>* energy_vector,
   vector<Double_t>* zenith_vector, vector<Double_t>* detector_vector, double detector_min, double detector_max, int particle_monto_carlo_number, TH1D *pointer1, TH1D *pointer2, TH1D *pointer3)
{
	vector<Double_t>* vector_with_ones = new vector<Double_t> (longitude_vec->size());
	for(unsigned n =0; n < longitude_vec->size(); n++)
	{
		vector_with_ones->at(n) = 1;
	}
	for (unsigned index =0; index < longitude_vec->size(); index++)
		{
		filter_and_filling(pointer2, -particle_monto_carlo_number, index, PDG_vector, energy_vector, vector_with_ones);
		filter_and_filling(pointer3, particle_monto_carlo_number, index, PDG_vector, energy_vector, vector_with_ones);
		if((latitud_mini <  latitude_vec->at(index) && latitude_vec->at(index) <  latitud_maxi) &&
		 (longitud_mini < longitude_vec->at(index) && longitude_vec->at(index)< longitude_maxi)  
		&& (detector_min < detector_vector->at(index) && detector_vector->at(index) < detector_max)
	 	&&(0.1E3 < energy_vector->at(index) && energy_vector->at(index) < 4.4E3)
		 && zenith_vector->at(index)<(0.451027))
		
		{

			filter_and_filling(pointer1, -particle_monto_carlo_number, index, PDG_vector, energy_vector, zenith_vector);

		}
	
	}
	return 0;
}

void iterate_over_file(TTree *propertie, double latitud_mini, double latitud_maxi, double longitud_mini, double longitude_maxi, vector<Double_t>* latitude_vec, vector<Double_t>* longitude_vec,vector<int>* PDG_vector, vector<Double_t>* energy_vector, vector<Double_t>* zenith_vector, vector<Double_t>* detector_vector, double detector_min, double detector_max, int particle_monto_carlo_number, TH1D *pointer1, TH1D *pointer2, TH1D *pointer3)
{

	for (float eventnumber =  0 ; eventnumber < propertie->GetEntriesFast();++eventnumber)
	{
		propertie->GetEntry(eventnumber);
		{
		filter_pointer(latitud_mini, latitud_maxi, longitud_mini, longitude_maxi, latitude_vec, longitude_vec,  PDG_vector, energy_vector,
		zenith_vector,  detector_vector, detector_min, detector_max, particle_monto_carlo_number, pointer1, pointer2, pointer3);
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



double integrated_differencial_energy(vector<Double_t> flux_vector)
{

	const double pi = 3.1416;
	double integrated_diffrensial_energy = 0;
	for(unsigned i =0; i<flux_vector.size(); i++)
	{	
		integrated_diffrensial_energy +=flux_vector[i];		
	}
	return pi * integrated_diffrensial_energy;
;
}



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

struct axis{Double_t axis_array[13];}xgraph, ygraph, bin_error; 

void get_axis(TH1D *histogram, int bin_numbers, double norm) 
{
	Double_t an_axis[bin_numbers];	
	TH1D *the_normalized_histogram = normalized_histogram(histogram, bin_numbers, norm);
	for (int n=0; n<bin_numbers; n++)
	{
		xgraph.axis_array[n]=the_normalized_histogram->GetBinCenter(n+1); 
		ygraph.axis_array[n]=the_normalized_histogram->GetBinContent(n+1);
		bin_error.axis_array[n] = the_normalized_histogram->GetBinError(n+1);
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

Double_t *ratio_error(int bin_numbers, Double_t *flux_data_error)
{
	Double_t *ratio_error_array = new Double_t[bin_numbers];
	for(auto n =0; n <bin_numbers; n++)
	{
	ratio_error_array[n] = bin_error.axis_array[n]/flux_data_error[n];
	}
	return ratio_error_array;
}
//Section 4 showing data, graphs, creating tabels 
TGraphAsymmErrors *create_error_graph(int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *yerrorhigh, Double_t *yerrorlow, Float_t color)
{
	Double_t xerrorlow[bin_numbers];
	Double_t xerrorhigh[bin_numbers];
	auto  *graph = new TGraphAsymmErrors(bin_numbers, xaxis, yaxis, xerrorlow,xerrorhigh, yerrorlow, yerrorhigh);
	graph->SetMarkerColor(color); 
	graph->SetMarkerSize(3);
	graph->SetMarkerStyle(22);
	graph->SetLineWidth(2);
	graph->SetLineColor(color);
	return graph;
}

struct statisticdata
{
	Double_t variabel;
}mean, rms;
void get_statistic_data(TGraphAsymmErrors *pointer)
{
	mean.variabel = pointer->GetMean(2);
	rms.variabel= pointer->GetRMS(2);
}
void create_graph(TString title, TString name, int bin_numbers, Double_t *xaxis, Double_t *yaxis,Float_t color)
{
	TGraph *graph = new TGraph(bin_numbers, xaxis, yaxis);
	graph->SetNameTitle(name, title);
	graph->SetMarkerColor(color);
	graph->SetMarkerSize(0.8);
	graph->SetMarkerStyle(21);
	graph->SetLineWidth(3);
	graph->Draw("ALP");

}

void create_canvas(TString name, TString title, TString name2, TString title2, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *yaxis2, Double_t *yerrorlow_stat, Double_t *yerrorhigh_stat, Float_t color, Float_t color2) 
{
	Double_t *ratio_error_high = ratio_error(bin_numbers, yerrorhigh_stat);
	Double_t *ratio_error_low = ratio_error(bin_numbers, yerrorlow_stat);


	Double_t *bin_error_array = new Double_t[bin_numbers];
	for(unsigned n = 0; n<bin_numbers; n++)
	{
		bin_error_array[n] = bin_error.axis_array[n];
	}
	TCanvas *canvas = new TCanvas("c3","plots energy",600,600,1800, 1200);
	canvas->Divide(2,1); 
	canvas->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	TMultiGraph *multipointer = new TMultiGraph();
	multipointer->SetNameTitle(name, title);
	TGraphAsymmErrors *pointer = create_error_graph(bin_numbers, xaxis, yaxis, yerrorlow_stat, yerrorhigh_stat,color);
	pointer ->Draw("P");
	pointer->SetTitle("experimental flux");

	TGraphAsymmErrors *pointer2 = create_error_graph(bin_numbers, xgraph.axis_array, ygraph.axis_array, bin_error.axis_array, bin_error.axis_array, color2);
	pointer2->Draw("P");
	pointer2->SetTitle("Simulated flux");
	multipointer->Add(pointer);
	multipointer->Add(pointer2);
	multipointer->Draw("ALP");
	multipointer->GetXaxis()->SetTitle("Kinetic energy GeV");
	multipointer->GetYaxis()->SetTitle("Antiproton flux(1/(sr s GeV m^2))" );
	canvas->Update();
	gPad->BuildLegend();
	canvas->Update();
	canvas->cd(2);
	gPad->SetLogx();
	TGraphAsymmErrors *ratio_pointer = create_error_graph(bin_numbers, xaxis, yaxis2, ratio_error_low, ratio_error_high, 1);
	ratio_pointer->Draw("ALP");
	ratio_pointer->SetTitle("ratio simulated experimental flux");
	gPad->BuildLegend();
	canvas->Update();
	canvas->Modified();

}


void create_and_write_file(string file_name,string witch_run, int bin_numbers, Double_t *indata_array1, Double_t *indata_array2, Double_t *indata_array3)
{
	
	ofstream myfile(file_name);
	ofstream outputFile;
	
	outputFile.open(file_name);	
	
	myfile<<witch_run<<" "<<"Kinetic energy"<<" "<<"positive stastistical error"<<"\n";
	myfile<<"(1/(sr s GeV m^2))"<<" "<<"GeV"<<" "<<"GeV"<<" "<<"GeV"<<"\n";
	

	for(int i = 0; i < bin_numbers; i++)
	{
		myfile<< indata_array1[i]<<" "<<indata_array2[i]<<" "<<indata_array3[i]<<"\n";
		
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

Double_t *converted_array(int bin_numbers, Double_t *start_array, Double_t conversion_factor)
{
	Double_t *modify_array;
	modify_array = new Double_t[bin_numbers];
	for(int i = 0; i < bin_numbers; i++)
	{
		modify_array[i] =  conversion_factor *start_array[i];	
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
	double detector_altitude=37;
	double longitud_max= 180;  //Range of laltitude and longitude of the detector
	double longitud_min= -180;
	double latitud_max= -75;
	double latitud_min= -85;
	//Some tabel values
	int nbins = 13;
	double area_conversation_factor = 10000;
	double the_solidangle=0.628319;
	//Define considered solid angle 
	//Definding binding range
	vector<Double_t> BESS_data_particle_flux_vector= get_one_vector_from_file(3,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
 	vector<Double_t> BESS_data_energybins_low_vector=get_one_vector_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
 	vector<Double_t> BESS_data_energybins_high_vector=get_one_vector_from_file(2,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
 	vector<Double_t> BESS_data_stat_error_low_vector=get_one_vector_from_file(4,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
 	vector<Double_t> BESS_data_stat_error_high_vector=get_one_vector_from_file(5,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/BESS2004.dat",13,6);
	vector<Double_t> BESS_bins_vector= modifiy_vector(BESS_data_energybins_low_vector, BESS_data_energybins_high_vector);
 	vector<Double_t> LIS_vector=get_one_vector_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/clusterSimScripts/spectrumtest550.dat",68,2);

	double integrated_LIS =  integrated_differencial_energy(LIS_vector);
	Double_t *BESS_data_particle_flux= BESS_data_particle_flux_vector.data();
	Double_t *BESS_bins=BESS_bins_vector.data();
	Double_t *BESS_stat_error_low = BESS_data_stat_error_low_vector.data();
	Double_t *BESS_stat_error_high = BESS_data_stat_error_high_vector.data();
	//TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary/";
	TString path = "/nfs/mds/project/e18sat/sim_ekemar_mds/supershorttest550/";
	const char * prefix = "tree";
	TString savepath = "manipulationstuff.root";
	const char *dirname;
	//const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/antiprotononsasprimary";
	const char *directory_name = "/nfs/mds/project/e18sat/sim_ekemar_mds/supershorttest550/";
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
	TH1D *flux_histogram = new TH1D("fluxene","calculated flux", nbins, BESS_bins);
	TH1D *secondary_energy = new TH1D("secondary_energy","secondary_energy",1000,0,1E3);
	TH1D *secondary_energy_antiprotons = histogram_setupp("energy_antiprotons", "Hits antiprotons", 1000, 0, 5E2);
	TH1D *secondary_energy_protons = histogram_setupp("energy_protons", "Hits protons", 1000, 5.6, 10E2);
	
	
	
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


		//data analysis: 


		total_events += properties->GetEntriesFast();


		double detector_area = ratio_area(startaltitude, detector_altitude, latitud_max,latitud_min, longitud_max, longitud_min);
		norm_factor = calculated_norm(the_solidangle ,area_conversation_factor, integrated_LIS, total_events, detector_area);

		iterate_over_file(properties, latitud_min, latitud_max, longitud_min, longitud_max, secondary_latitude_vec, secondary_longitude_vec, secondary_PDG_vec, secondary_Energy_vec,
		 secondary_zenith_vec, secondary_boundDet_vec, 90, 105, 2212, flux_histogram, secondary_energy_antiprotons,secondary_energy_protons);

		
	//Clean up
	properties->ResetBranchAddresses();
	infile->Close();

	}
	
 		
				
	outfile = new TFile(savepath,"RECREATE");
	secondary_energy_protons->Write("proton second");
	flux_histogram->Write("test");
	secondary_energy_antiprotons->Write("antiproton");

	check_content_vector("BESS energy bins", nbins+1, BESS_bins_vector);
	check_content_vector("BESS particle flux", nbins, BESS_data_particle_flux_vector);
	
	get_axis(flux_histogram,nbins, norm_factor);
	check_content_array("XGraph", nbins,  xgraph.axis_array);
	check_content_array("YGraph", nbins +1,  ygraph.axis_array);
	check_content_array("Bin errors: ", nbins,  bin_error.axis_array);
	Double_t *firstratio = ratio_calculations(nbins-1, BESS_data_particle_flux);
	TGraphAsymmErrors * graph = create_error_graph(nbins, xgraph.axis_array, ygraph.axis_array,  bin_error.axis_array, bin_error.axis_array, 3);
//	get_statistic_data(graph);
	cout<<"Mean: "<<mean.variabel<<endl;
	cout<<"RMS: "<<rms.variabel<<endl;
	cout<<"Integrated diffrential energy: "<<integrated_LIS<<endl;
	create_canvas("Flux", "flux", "Ratio", "Ratio" , nbins, xgraph.axis_array, BESS_data_particle_flux, firstratio, BESS_stat_error_low, BESS_stat_error_high, 4, 3);


	Double_t *ratio = ratio_calculations(nbins, BESS_data_particle_flux);
	Double_t *the_ratio_error = ratio_error(nbins, BESS_stat_error_low);

	string save_file = "/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/theisis/datafiles/550test.dat";
	string save_file2 = "/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/theisis/datafiles/ratio550test.dat";
	string save_file_nofilter = "/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/theisis/datafiles/550testnofilter.dat";
	string save_file2_nofilter = "/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/theisis/datafiles/ratio550testnofilter.dat";

	create_and_write_file(save_file, "secundary  flux test",nbins, ygraph.axis_array, xgraph.axis_array, bin_error.axis_array);
	create_and_write_file(save_file2, "ratio",nbins, ratio, xgraph.axis_array,the_ratio_error);

	memset(BESS_bins, 0, 14*sizeof(*BESS_bins));
	memset(BESS_data_particle_flux, 0, 13*sizeof(*BESS_data_particle_flux));
	memset(firstratio, 0, 13*sizeof(*firstratio));
	memset(bin_error.axis_array, 0, 13*sizeof(*bin_error.axis_array));
	outfile->Write();
	outfile->Close();
	
	cout<<"run done"<<endl;
	theApp.Run();

	return 0;

}
