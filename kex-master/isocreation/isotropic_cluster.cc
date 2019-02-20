/*
-This code generates PLANETOCOSMICS macros for an isotropic flux inside a spherical starting-positions sphere.
-Description in isotropicstart_macgenerator_README.txt
-Author: Andrea Meraner (Sep. 2016) andrea.meraner@tum.de
*/
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
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>




using namespace std;
	

int main(int argc, char* argv[]){
	
		if(argc != 4){
			cout<<argc<<endl;
			cout<<"usage: <startid> <number of macros to be created> <baseclassfile>"<<endl;
			exit(0);
		}

//Define the section of the sphere (in latitude and longitude degree limits) in which the starting positions will be generated. Values for the entire sphere: -90, 90, -180, 180.
	int minlat=-90; //[-90, 90)
	int maxlat=90; //(-90,90]
	int minlon=-180; //[-180, 180)
	int maxlon=180; //(-180,180]

//Define the altitude -which is not the radius- of the starting-positions sphere in km.
	double altitude=70000; //A higher value gives a better magnetic field simulation but requires higher computation time, as less particles will hit the Earth. 70000 km proved to be a good compromise.

	
//Define starting positions properties
	int startpos=1000; //Number of starting positions per macro
	int part=1; //Number of particles generated at each starting position. The angle distribution will be a cosine-law, as implemented in PLANETOCOSMICS.

//Define if macro should be set for a cluster run
	bool cluster_run = true;
	
//Define the name and number of macros to be created
	int macnumb= atoi(argv[2]) ;
	int macnumbOffset = atoi(argv[1]);
	string macroname="startTest_";
	
//Define if macro contains the tree generator command
	bool tree_run=true;

//Define the input file
	ifstream filein(argv[3]); //Source file

//Define file for saving the maps
//	TFile* mapsFile= new TFile("mapstest1.root","RECREATE"); 

//Define distribution
	int nbinlat=(maxlat-minlat);
	int nbinlon=(maxlon-minlon);
	TH2D* latlondist = new TH2D("latlon","lat and long distribution of starting positions",nbinlon,minlon,maxlon,nbinlat,minlat,maxlat);

	TH2D* latlonpos = (TH2D*)latlondist->Clone(); //Clone empty map
	latlonpos->SetNameTitle("starting_positions","lat and long coordinates of starting positions");

//Fill distribution map
	for (int x=minlon; x<=maxlon;x++){
		for (int y=minlat; y<=maxlat; y++){
			double c=3.14159/180;
			latlondist->Fill(x,y,0.5*cos(y*c)); //Cosine distribution to regard the latitude-longitude coordinates
		}
	}
	cout<<"Distribution set"<<endl;
//	latlondist->Write();
	
	ostringstream sstream; //stringstream to convert doubles in strings
	string latStr, lonStr, altStr, partStr, startStr, startPart, lineTemp;
	double lat,lon;
	
	
//----------------------Macro loop start------------------------------------------------------------------------------------
	for (int macn=macnumbOffset; macn<=macnumb; macn++){
	
		sstream << macn;
		string macnStr=sstream.str();
		sstream.str( std::string() );
		sstream.clear(); 
		
		string currentfileout= macroname + macnStr + ".g4mac";
		ofstream fileout(currentfileout.c_str()); //Numbered output file
		cout<<"New macro created: "<<currentfileout<<endl;
		
		
			//Seek to starting positions part in macro and copy test until then
				string lineSeek="#StartingPositions";
	
				bool found = false;

				while(getline(filein,lineTemp))
				{
					if(lineTemp == lineSeek){
						found = true;
					}
					lineTemp += "\n";
					fileout << lineTemp;
					if(found) break;
				}

			//Compute starting positions and add lines to macro
				for (int i=0; i<startpos; i++) {
					latlondist->GetRandom2(lon, lat); //get lat lon according to distribution
//					latlonpos->Fill(lon, lat); //fill starting positions map
					//cout<<lat<<endl;
					//cout<<lon<<endl;
					//cout<<altitude<<endl;
					//cout<<endl;
	
					sstream << lat;
					latStr=sstream.str();
					sstream.str( std::string() );
					sstream.clear();
		
					sstream << lon;
					lonStr=sstream.str();
					sstream.str( std::string() );
					sstream.clear();
		
					sstream << altitude;
					altStr=sstream.str();
					sstream.str( std::string() );
					sstream.clear();
		
					sstream << part;
					partStr=sstream.str();
					sstream.str( std::string() );
					sstream.clear();
		
					//Add lines to macro
					startStr= "/PLANETOCOS/SOURCE/SetPosition PLA " + altStr + " km " + latStr + " " + lonStr + " degree \n"; 
					fileout << startStr;
					startPart="/run/beamOn "+ partStr + "\n \n";
					fileout << startPart;
	
				}
			
			//Add tree save command
			if (tree_run) { 
				string treecommand="#SaveTree \n/PLANETOCOS/ANALYSIS/COMPLETE_EVENT/WriteFluxRoot \n \n";
				fileout << treecommand;
			}
	
			//Seek to final save part in macro and copy remaining text (only if not in cluster mode)
				lineSeek="#Save";
				while(getline(filein,lineTemp)){
					if(lineTemp == lineSeek){
						lineTemp += "\n";
					 	fileout << lineTemp; 
					 	
					 	if (!cluster_run){
							while(getline(filein,lineTemp)){
								lineTemp += "\n";
								fileout << lineTemp;  
							}
						}
					}
				}
				
		cout<<"Creation of macro ended"<<endl;
		filein.clear();
		filein.seekg(0,ios::beg); //Return at the start of the input file for the next iteration
		
		fileout.close();
		
	}
//-------------Macro loop ends-------------------------------------------------------------------------------------------------------------
	
	
	
	//latlonpos->Write(); //save starting positions map

	
	//TCanvas *c1 = new TCanvas("maps","distribution maps",50,50,1600, 600);
	//c1->Divide(2,1);
	//c1->cd(1);
    //latlondist->Draw("colz");
    //c1->cd(2);
    //latlonpos->Draw("colz"); 
	return 1;
}
