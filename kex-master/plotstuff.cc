//Liselotts code for calculating the flux from antiprotons
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
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

vector<double> get_energy_bins_from_file(string filename,int rows, int columns)
{

     string line;
     int row,col;
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

                 switch(col)
                 {
                 case  1:
                 lower_energy.push_back(my_array[row][1]);
		 break;
                 case  2:
                 higher_energy.push_back(my_array[row][2]);
		 break;
                 }
                 col++;
             }
             row++;
         }
          pFile.close();

	}
	return all_energies;
}

/*TCanvas plot_stuff(TString name, TString Title, Double_t x-axis, vector<Double_t> y-axis, Int_t bin_numbers)
{
	TCanvas *name = new TCanvas(name, Title, x-axis, y-axis.size()):
	TGraph *graph = new TGraph(bin_numbers, x-axis, y-axis);
	graph->Draw();
		
}*/
int main(int argc, char* argv[]){
	TApplication theApp("App", &argc, argv);
	vector<Double_t> BESS2007_particle_flux_vector=get_data_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6);
	vector<Double_t> BESS2007_energybins_vector=get_data_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6);
	Double_t *BESS2007_particle_flux= BESS2007_particle_flux_vector.data();
	Double_t *BESS2007_energybins=BESS2007_energybins_vector.data();
	cout<<"Particle flux"<<endl;
	for(int i = 0; i <BESS2007_particle_flux_vector.size();i++)
	{
		cout<<BESS2007_particle_flux_vector.at(i)<<" ";
	}
	cout<<"\n"<<endl;
	cout<<"Particle energy bins"<<endl;
	for(int i = 0; i <BESS2007_energybins_vector.size();i++)
	{
		cout<<BESS2007_energybins_vector.at(i)<<" ";
	}
	
	TCanvas *flux_graph = new TCanvas("flux_graph", "particle flux", 600, 600, 700, 500);
	Int_t bin_numbers = 29;
	Double_t x[29];
	//TGraph *plotfromfile = new TGraph("/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat", "%lg %lg");
	for(int i =0; i <bin_numbers; i++)
	{
		x[i] = 0.1*i;	
	}
	//TGraph *plotfromfile = new TGraph(bin_numbers, BESS2007_energybins, BESS2007_particle_flux);
	TGraph *plotfromfile = new TGraph(bin_numbers, x, BESS2007_particle_flux);
	plotfromfile->Draw();
	cout<<"Hej"<<endl;
	memset(BESS2007_energybins, 0, 30*sizeof(*BESS2007_energybins));
	memset(BESS2007_particle_flux, 0, 29*sizeof(*BESS2007_particle_flux));
	theApp.Run();
}
