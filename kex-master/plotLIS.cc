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
#include <TApplication.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < int >+;
#pragma link C++ class std::vector < Double_t >+;
#endif 

using namespace std;

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



void check_content_array(string name, int bin_numbers, Double_t *an_array )
{

	cout<<name<<" ";
	for(int i = 0; i< bin_numbers; i++)
	 {
		cout<<an_array[i]<<", ";
	}
	cout<<"\n";
}
void check_content_vector(string name, int bin_numbers, vector<Double_t> a_vector )
{

	cout<<name<<" ";
	for(int i = 0; i< bin_numbers; i++)
	 {
		cout<<a_vector.at(i)<<", ";
	}
	cout<<"\n";
}


TGraph *create_graph(TString title,int bin_numbers, Double_t *xaxis, Double_t *yaxis,Float_t color, int marker)
{
	TGraph *graph = new TGraph(bin_numbers, xaxis, yaxis);
	graph->SetTitle(title);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerSize(0.8);
	graph->SetMarkerStyle(marker);
	graph->SetLineWidth(2);
	graph->Draw("LP");

	return graph;
}


	

void draw_graph(TString name, TString title,TString graph_title, TString graph_title2, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *xaxis2, Double_t *yaxis2, Float_t color, Float_t color2, int marker1, int marker2)
{
	TMultiGraph *multipointer = new TMultiGraph();
	multipointer->SetNameTitle(name, title);
	TGraph *graph_pointer = create_graph(graph_title, bin_numbers, xaxis, yaxis, color, marker1);
	TGraph *graph_pointer2 = create_graph(graph_title2, bin_numbers, xaxis2, yaxis2, color2, marker2);
	graph_pointer->Draw("P");
	graph_pointer2->Draw("P");
	multipointer->Add(graph_pointer);
	multipointer->Add(graph_pointer2);
	multipointer->Draw("ALP");
  	multipointer->GetYaxis()->SetTitleOffset(1.5);
  	//multipointer->GetYaxis()->SetTitleSize(0.9);
	multipointer->GetXaxis()->SetTitle("Kinetic energy (MeV)");
	multipointer->GetYaxis()->SetTitle("Differential energy spectrum, (1/(sr s MeV m^2)) ");
	

}
void create_canvas(TString name, TString title, int bin_numbers, Double_t *xaxis, Double_t *yaxis, Double_t *xaxis2, Double_t *yaxis2, Float_t color, Float_t color2, int marker1, int marker2) 
{

	TCanvas *canvas = new TCanvas("c3","plots energy",900,900,900, 900);
	canvas->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	draw_graph(name, title,"LIS solar modulation 615 MeV", "LIS solar modulation 550 MeV", bin_numbers, xaxis, yaxis,xaxis2, yaxis2,  color, color2, marker1, marker2);
	gPad->BuildLegend();
	canvas->Update();
	canvas->Modified();

}

Double_t *calculate_differenses(int bin_numbers, Double_t *array1, Double_t *array2)
{
	Double_t *differens_vector = new Double_t[bin_numbers];	
	for(auto n = 0; n <bin_numbers; n++)
	{
		differens_vector[n] = array1[n] - array2[n];	
	}
	return differens_vector;
}

void create_and_write_file(string file_name,int bin_numbers, Double_t *indata_array1)
{
	
	ofstream myfile(file_name);
	ofstream outputFile;
	
	outputFile.open(file_name);	
	
	myfile<<"Differences LIS  "<<"\n";
	myfile<<"(1/(sr s MeV m^2)) "<<"\n";
	

	for(int i = 0; i < bin_numbers; i++)
	{
		myfile<< indata_array1[i]; 
		
	}
	outputFile.close();
}

int main(int argc, char* argv[])
{
	TApplication theApp("App", &argc, argv);
	TFile *outfile;
	int nbins = 62;
	vector<Double_t> energy_vector_550= get_one_vector_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/clusterSimScripts/spectrumtest550.dat",67,2);
 	vector<Double_t> LIS_vector_550=get_one_vector_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/clusterSimScripts/spectrumtest550.dat",67,2);
	vector<Double_t> complete_energy_vector_550 = energy_vector_550;
	complete_energy_vector_550.push_back(1.1601e+08);
	Double_t * energy_550 = new Double_t[energy_vector_550.size()];
	Double_t * LIS_550 = new Double_t[LIS_vector_550.size()];
	energy_550 = energy_vector_550.data();
	LIS_550 = LIS_vector_550.data();



	vector<Double_t> energy_vector_615= get_one_vector_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/clusterSimScripts/spectrumtest615.dat",67,2);
 	vector<Double_t> LIS_vector_615=get_one_vector_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/clusterSimScripts/spectrumtest615.dat",67,2);
	vector<Double_t> complete_energy_vector_615 = energy_vector_615;
	complete_energy_vector_615.push_back(1.1601e+080);
	Double_t * energy_615 = new Double_t[energy_vector_615.size()];
	Double_t * LIS_615 = new Double_t[LIS_vector_615.size()];
	energy_615 = energy_vector_615.data();
	LIS_615 = LIS_vector_615.data();
	Double_t *diff_array = new Double_t[nbins];
	diff_array = calculate_differenses(nbins, LIS_550,LIS_615);
	check_content_array("Differenses LIS: ", nbins, diff_array);
	string file = "/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/theisis/datafiles/differensesLIS.csv";
	create_and_write_file(file, nbins, diff_array);
	outfile = new TFile("test.root","Update");
	create_canvas("LIS ", "Local interstellar flux ", nbins, energy_615, LIS_615,energy_550, LIS_550, 2, 3, 22, 23);
	memset(energy_615, 0, 62*sizeof(*energy_615));
	memset(energy_550, 0, 62*sizeof(*energy_550));
	memset(LIS_615, 0, 62*sizeof(*LIS_615));
	memset(LIS_550, 0, 62*sizeof(*LIS_550));
	outfile->Write();
	outfile->Close();
	theApp.Run();
}
