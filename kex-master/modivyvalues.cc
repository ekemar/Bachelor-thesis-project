//Code for recaluclate given data values
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
 vector<Double_t> get_data_from_file(int interesting_col_number,string filename, int rows, int columns)
 {
      string line;
      int row,col;
      std::vector<Double_t> particle_flux;
      std::vector<Double_t> energy;

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
			if(col == 1)
			{
                    	  energy.push_back(my_array[row][1]);
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
                 case 1:
                 return energy;


         }

}

void check_content_vector(string name, vector<Double_t> a_vector )
{

	cout<<name<<endl;
	for(int i = 0; i< a_vector.size(); i++)
	 {
		cout<<a_vector.at(i)<<", ";
	}
	cout<<"\n";
}

vector<Double_t> convert_units(vector<Double_t> my_vector, long double converter_unit)
{
	vector<Double_t> modified_vector (my_vector.size());
	for(int i = 0; i<my_vector.size(); i++)
	{
		modified_vector.at(i) = my_vector.at(i) * converter_unit;

	}
	return modified_vector;
}

void read_vector_to_file(string file_name, vector<Double_t> indata_vector1, vector<Double_t> indata_vector2)
{
	ofstream myfile(file_name);
	ofstream outputFile;
	
	outputFile.open(file_name);	
	


	for(int i = 0; i < indata_vector1.size(); i++)
	{
		myfile<< indata_vector1.at(i)<<" "<<indata_vector2.at(i)<<"\n";
		
	}
	outputFile.close();
}
int main(int argc, char* argv[])
{
	long double elementary_charge = 1 / 1.6021766208e-19;
	vector<Double_t> energy_vector = get_data_from_file(0,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6); 
	cout<<"Vector size: "<<energy_vector.size()<<endl;
	vector<Double_t> flux_vector = get_data_from_file(1,"/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/galpropmodified.dat",29,6); 
	cout<<"Vector size flux: "<<flux_vector.size()<<endl;
	vector<Double_t> energy_GeV_vector = convert_units(flux_vector, elementary_charge);
	vector<Double_t> GeV_vector = convert_units(energy_vector, elementary_charge);
	//
	check_content_vector("Experimental flux", flux_vector);
	check_content_vector("Flux converted to GeV", GeV_vector);
	check_content_vector("Energy bins", energy_vector);
	check_content_vector("Energy bins GeV", energy_GeV_vector);

	read_vector_to_file("/nfs/freenas/tuph/e18/project/e18sat/ekemar_sim/DATA/GeVgalpropmodified.dat", energy_GeV_vector, GeV_vector);

	fill(flux_vector.begin(), flux_vector.end(), 0);
	fill(GeV_vector.begin(), GeV_vector.end(), 0);
	fill(energy_vector.begin(), energy_vector.end(), 0);
	fill(energy_GeV_vector.begin(), energy_GeV_vector.end(), 0);
}
