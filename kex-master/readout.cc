/*
I 'm cleaning andreas code so I understand what happens
*/

{
	#include <vector>
	#include <iostream>
	#include <fstream>
	#include <string>
	#include <sstream>
	#include <math.h>
	#include <cmath> 
	using namespace std;
	

	//Define data path=====================================================================
	TString path="/nfs/mds/project/e18sat/sim_meraner_mds/magfieldstart/data8_2/";	
	
	//Define prefix of data that shall be considered
	string prefix="tree";
	
	//Define parameters for normalisation==================================================
	
	double startalt=70000; //Altitude of the primary sphere
	double detector_altitude1=25.2842;
	double detector_altitude2=25.2842;
	double detector_longitude_max1=-92;  //Range of latitude and longitude of the detector
	double detector_longitude_min1=-116;
	double detector_latitude_max1=39;
	double detector_latitude_min1=31;
	
	
	
	//Analytical factors
	double geomfac1=4*3.1416*pow(6371+startalt,2)/(pow(6371+detector_altitude1,2)*(detector_longitude_max1*0.017453-detector_longitude_min1*0.017453)*(sin(detector_latitude_max1*0.017453)-sin(detector_latitude_min1*0.017453)));
	double geomfac2=4*3.1416*pow(6371+startalt,2)/(pow(6371+detector_altitude2,2)*(detector_longitude_max1*0.017453-detector_longitude_min1*0.017453)*(sin(detector_latitude_max1*0.017453)-sin(detector_latitude_min1*0.017453)));
	
	//Define integrated primary flux
	double integrated_primary_flux= 0.504709; 
	
	//Define considered solid angle 
	double solidangle=0.628319;
	
	//Define names of histograms==========================================================
	TString name1="Flux of protons at 26.4 g/cm2";
	TString name2="Flux of muons at 26.4 g/cm2";
	
	
	//Plots and saving setup==============================================================
	
	//Define if plots should be saved
	bool saveplots=true;

	//Set limit to files to be analysed (set 0 for all files)
	int filelimit =25;
	
	//Define file name where the extracted data will be saved
	TString savepath="ProtonMuonQGSP2001high_26-4_26-4_921163931_12812couts.root";
	
	
	//Define frequency of live update
	const Int_t kUPDATE =20;	
	
    //Define range of map plots
	int minlat=-90;
	int maxlat=90;
	int minlon=-180;
	int maxlon=180;
	
	
	//Canvas setup========================================================================
	
	
	

	TCanvas *c3 = new TCanvas("c3","plots energy",300,300,900, 600);
	
	TCanvas *c4 = new TCanvas("c4","ratio",300,300,800, 600);
	
	
	TCanvas *c32 = new TCanvas("c32","plots energy2",300,300,900, 600);
	
	TCanvas *c42 = new TCanvas("c42","ratio2",300,300,800, 600);
	
	
	
	
	//Histograms setup==================================================================
	
	

	//BESS data binning ranges
	Double_t bessx[31] = {0.215,0.251,0.293,0.342,0.398,0.464,0.541,0.631,0.736,0.858,1.000,1.17,1.36,1.58,1.85,2.15,2.51,2.93,3.42,3.98,4.64,5.41,6.31,7.36,8.58,10.0,11.7,13.6,15.8,18.5,21.5}; //BESS2000 TOA Flux energy axis (Shikaze et.al)
	Double_t bessxd[32]={0.185,0.215,0.251,0.293,0.341,0.398,0.464,0.541,0.631,0.736,0.858,1.000,1.17,1.36,1.58,1.85,2.15,2.51,2.93,3.41,3.98,4.64,5.41,6.31,7.36,8.58,10.0,11.7,13.6,15.8,18.5,21.5}; //BESS2000 atmo flux energy axis (Shikaze et.al)

	Double_t bessx01prot[21]={0.46,0.54,0.63,0.74,0.8,1.00,1.17,1.36,1.58,1.85,2.15,2.51,2.93,3.41,3.98,4.64,5.41,6.31,7.36,8.58,10.00}; //BESS2001 proton flux at 26.4 g/cm2 energy axis (Abe et.al)
	Double_t bessxmuon[21]={0.50, 0.58, 0.67, 0.78, 0.90, 1.05, 1.21, 1.41, 1.63, 1.90, 2.20, 2.55, 2.96, 3.44, 3.99, 4.63, 5.38, 6.24, 7.25, 8.41,9.76}; //BESS2001 muonflux at 26.4 g/cm2 momentum axis (Abe et.al)
	
	int nbins=20; //x axis vector size -1
	int nbins2=20;
		
	TH1F *fluxene = new TH1F("fluxene",name1,nbins,bessx01prot);
	TH1F *fluxenecounts = new TH1F("fluxenecounts","Counts1", nbins, bessx01prot);
	TH1F *fluxene2 = new TH1F("fluxene2",name2,nbins2,bessxmuon);
	TH1F *fluxenecounts2 = new TH1F("fluxenecounts2","Counts2", nbins2, bessxmuon);
	
	
	//TGraph setup====================================================================
	
	Double_t xgraph[nbins], ygraph[nbins], rat[nbins], rat2[nbins2], xgraph2[nbins2], ygraph2[nbins2]; //X axis initialization
	
	TMultiGraph *fluxes= new TMultiGraph();
	fluxes->SetNameTitle("fluxes","BESS Flux & PCsim Flux");
	
	TMultiGraph *fluxes2= new TMultiGraph();
	fluxes2->SetNameTitle("fluxes2","BESS Flux & PCsim Flux 2");
	
	Double_t ybessgraph[30]={169,192,206,213,232,233,243,246,248,249,239,234,212,198,177,155,134,114,95.5,77.5,62.4,49.7,37.4,27.0,20.4,14.2,10.6,7.22,5.15,3.60};//BESS2000 TOA Flux
	Double_t ybessgraphd[31]={789,786,792,646,560,504,461,436,358,366,351,302,273,329,170,256,162,139,98.9,103,64.6,66.6,34.8,29.1,25.9,19.9,12.2,8.58,5.82,4.42,2.44};//BESS2000 atmo flux

	Double_t ybessmuon[20]={57, 50.7, 31.5, 30.8, 22.6, 19.3, 17.4, 13.5, 10.3, 8.31, 5.83, 4.42, 3.97, 2.21, 1.83, 1.15, 0.840, 0.405, 0.404, 0.301}; //BESS2001 muon flux at 26.4 g/cm2
	Double_t ybess01prot25km[20]={178.0,141.0,124.0,103.0,83.0,74.9,71.7,58.5,50.9,41.5,40.3,45.2,68.2,91.2,75.4,63.3, 47.9, 34.7, 26.8, 17.6}; //BESS2001 proton flux at 26.4 g/cm2

	
	//Initialize vectors for output of the tree ============================================
	vector<int>* primary_PDG_vec;
	vector<Double_t>* primary_Energy_vec;
	vector<Double_t>* primary_zenith_vec;
	vector<Double_t>* primary_azimuth_vec;
	vector<Double_t>* primary_latitude_vec;
	vector<Double_t>* primary_longitude_vec;

	vector<int>* secondary_PDG_vec;
	vector<Double_t>* secondary_Energy_vec;
	vector<Double_t>* secondary_zenith_vec;
	vector<Double_t>* secondary_azimuth_vec;
	vector<Double_t>* secondary_latitude_vec;
	vector<Double_t>* secondary_longitude_vec;
	vector<Double_t>* secondary_boundDet_vec;
	
	TString fname;
	long long int sumentries=0;
	int sumentries1=0;
	int sumentries2=0;
	double normfac1,normfac2,mom;
	bool stop=false;
	
	
//ANALYSIS START==============================================================================================================================================

//Opening folder and reading files
	TSystemDirectory dir(path,path);
    TList *files = dir.GetListOfFiles();
 	if (files) {
  		files.Sort();
        TSystemFile *file;
        TIter next(files);
        int count=1;
        
        while (true) {
        	
//Analyze data and update histograms during execution=====================================================================
				if (count%kUPDATE == 0 || stop) {
				
					cout<<"Canvas update - Data after "<<count<<" files"<<endl;
					cout<<"Primaries: "<<sumentries<<endl;
					
					//First analysis+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					cout<<name1<<": "<<endl;
					cout<<"Entries: "<<sumentries1<<endl;
					normfac1= (1/solidangle)*10000*((integrated_primary_flux/sumentries)*(geomfac1));  //The factor 10000 accounts for the area unit conversion (from cm2 to m2)
					cout<<"Normalisation factor: "<<normfac1<<endl;
				
	
					c3->cd();  //Change to energy canvas
	
					gPad->SetLogx();  //Set logarithmic axes
					gPad->SetLogy();
				
					TH1F* fluxenecopy = (TH1F*)fluxene->Clone(); //Create a working copy of the current histogram and modify only on this copy. This is needed for the live data show. The original histogram shall not be modified. 
					fluxenecopy->SetNameTitle("fluxenecopy", "Flux");
					
				    //Fill the working histogram
					for (int n=0; n<nbins; n++){
	
						double val=fluxene->GetBinContent(n+1); 
						fluxenecopy->SetBinContent(n+1, val);
					}
				
					fluxenecopy->Scale(normfac1, "width");	//Normalization of the histogram. The bin contents are divided by the width of the bin to achieve an x-axis specific value. 
				    
				    //Fill the graph arrays with the normalized histogram data	
					for (int n=0; n<nbins; n++){
						xgraph[n]=fluxenecopy->GetBinCenter(n+1); 
						ygraph[n]=fluxenecopy->GetBinContent(n+1);
						rat[n]=ygraph[n]/ybess01prot25km[n]; //Save the ratio values between simulation and experimental data.
					}
				
					fluss= new TGraph(nbins,xgraph,ygraph); //create graph with simulation data
					fluss->SetMarkerColor(4);
					fluss->SetMarkerSize(0.8);
					fluss->SetMarkerStyle(21);
	
					bessfluss= new TGraph(nbins,xgraph,ybess01prot25km); //create graph with experimental data
					bessfluss->SetMarkerColor(3);
					bessfluss->SetMarkerSize(0.8);
					bessfluss->SetMarkerStyle(21);
				
					fluxes->Add(fluss); //add graphs to multipgraph
					fluxes->Add(bessfluss);
					fluxes->Draw("ALP");
				
					c3->Update();
			
					c4->cd(); //Change to the ratio canvas
					gPad->SetLogx();
			
					ratio=new TGraph(nbins, xgraph, rat);
					ratio->SetNameTitle("ratio", "error ratio");
					ratio->SetMarkerColor(4);
					ratio->SetMarkerSize(0.8);
					ratio->SetMarkerStyle(21);
					ratio->Draw("ALP");
				
					c4->Update();
				
					//Second analysis (same steps as above)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					cout<<name2<<": "<<endl;
					cout<<"Entries: "<<sumentries2<<endl;
					normfac2= (1/(solidangle))*10000*((integrated_primary_flux/sumentries)*(geomfac2));
					cout<<"Normalisation factor: "<<normfac2<<endl;
				
				
	
	
	
					c32->cd();
	
					gPad->SetLogx();
					gPad->SetLogy();
				
					TH1F* fluxene2copy = (TH1F*)fluxene2->Clone();
					fluxene2copy->SetNameTitle("fluxene2copy", "Flux ene 2");
				
					for (int n=0; n<nbins2; n++){
	
						double val2=fluxene2->GetBinContent(n+1); 
						fluxene2copy->SetBinContent(n+1, val2);
					}
				
					fluxene2copy->Scale(normfac2, "width");	
								    	
					for (int n=0; n<nbins2; n++){
	
						xgraph2[n]=fluxene2copy->GetBinCenter(n+1); 
		
						ygraph2[n]=fluxene2copy->GetBinContent(n+1); 
					
						rat2[n]=ygraph2[n]/ybessmuon[n];
					}
				
				
					fluss2= new TGraph(nbins2,xgraph2,ygraph2);
					fluss2->SetMarkerColor(4);
					fluss2->SetMarkerSize(0.8);
					fluss2->SetMarkerStyle(21);
	
					bessfluss2= new TGraph(nbins2,xgraph2,ybessmuon);
					bessfluss2->SetMarkerColor(3);
					bessfluss2->SetMarkerSize(0.8);
					bessfluss2->SetMarkerStyle(21);
				
					fluxes2->Add(fluss2);
					fluxes2->Add(bessfluss2);
					fluxes2->Draw("ALP");
				
					c32->Update();
			
					c42->cd();
					gPad->SetLogx();
					

					ratio2=new TGraph(nbins2, xgraph2, rat2);
					ratio2->SetNameTitle("ratio2", "error ratio2");
					ratio2->SetMarkerColor(4);
					ratio2->SetMarkerSize(0.8);
					ratio2->SetMarkerStyle(21);
					ratio2->Draw("ALP");
				
					c42->Update();
				

					
					if (stop){ //Save plots if the filelimit has been reached or the file list has ended
						if (saveplots){
							TFile* mapsFile= new TFile(savepath,"UPDATE");
						
							fluxes->Write();
							fluxenecounts->Write();
							ratio->Write();
							fluxes2->Write();
							fluxenecounts2->Write();
							ratio2->Write();
							
							mapsFile->Close();
						}
						cout<<"Analysis ended!"<<endl;
						break; //Stop analysis
					}
						
					//Delete working histograms and graphs before the next iteration
					delete fluss;
					delete bessfluss;
					delete fluxenecopy;
					delete fluss2;
					delete bessfluss2;
					delete fluxene2copy;
					delete ratio;
					delete ratio2; 
					
				}




//DATA FILTERING and EXTRACTION=========================================================================================
        
				if ((file=(TSystemFile*)next()) ) {  
				
				    fname = file->GetName();
				    if (fname.BeginsWith(prefix.c_str())){  //prefix check
				    	
				
						if (filelimit!=0 && (count)==filelimit) {  //filelimit check
							cout<<"Analysis stopping at file number "<<filelimit<<endl;
							stop=true;
						}
						
						
				    	count++;
				    	
				    	
			   			cout<<"Reading: "<<fname<<endl;
				    	fname=path+fname;
				 		TFile *infile = TFile::Open(fname);
				 		
						//Read in tree from file
						TTree *properties = (TTree*)infile.Get("properties");

						//Set Branch Addresses to the vectors
						properties->SetBranchAddress("primary_PDGcode",&primary_PDG_vec);
						properties->SetBranchAddress("primary_Energy",&primary_Energy_vec);
						properties->SetBranchAddress("primary_zenith",&primary_zenith_vec);
						properties->SetBranchAddress("primary_azimuth",&primary_azimuth_vec);
						properties->SetBranchAddress("primary_latitude",&primary_latitude_vec);
						properties->SetBranchAddress("primary_longitude",&primary_longitude_vec);
	
						properties->SetBranchAddress("secondary_PDGcode",	&secondary_PDG_vec);
						properties->SetBranchAddress("secondary_Energy",	&secondary_Energy_vec);
						properties->SetBranchAddress("secondary_zenith",	&secondary_zenith_vec);
						properties->SetBranchAddress("secondary_azimuth",	&secondary_azimuth_vec);
						properties->SetBranchAddress("secondary_latitude",		&secondary_latitude_vec);
						properties->SetBranchAddress("secondary_longitude",		&secondary_longitude_vec);
						properties->SetBranchAddress("secondary_BoundaryDetector", &secondary_boundDet_vec);
	
	
						sumentries+=properties->GetEntries();

						for (int a = 0 ; a<properties->GetEntries();++a){
							properties->GetEntry(a);
						
							//loop over secondaries
							for (int i=0; i<secondary_longitude_vec.size(); i++){
							
							//FILTERS===================================================================================================================================
								//detector area filter
								if ((detector_longitude_min1<secondary_longitude_vec.at(i) && secondary_longitude_vec.at(i)<detector_longitude_max1)) {
										if ((detector_latitude_min1<secondary_latitude_vec.at(i) && secondary_latitude_vec.at(i)<detector_latitude_max1)){
									
												//detector altitude filter
												if (secondary_boundDet_vec.at(i)<25.3 && secondary_boundDet_vec.at(i)>25.2){
												
													//incident angle filter
													if(secondary_zenith_vec.at(i)<(0.451027)){
															
																//First analysis++++++++++++++++++++++++++++++++++++++++++++++++++++
																//particle type filter 
																if (secondary_PDG_vec.at(i)==2212){ //only protons
															
																	//energy filter
																	if (460<=secondary_Energy_vec.at(i) && secondary_Energy_vec.at(i)<=10000) {
																			sumentries1++;
																		fluxene->Fill(secondary_Energy_vec.at(i)/1000, 1/cos(secondary_zenith_vec.at(i))); //cosine weight for planar geometry
																			fluxenecounts->Fill(secondary_Energy_vec.at(i)/1000); //MeV to GeV unit conversion
																		}
																	}
																
																//Second analysis+++++++++++++++++++++++++++++++++++++++++++++++++++
																//particle type filter
																if (secondary_PDG_vec.at(i)==13){ //only muons
															
																	//relativistic momentum calculation
																	mom=sqrt(pow(secondary_Energy_vec.at(i),2)+2*secondary_Energy_vec.at(i)*105.658);
																
																	//momentum filter
																	if (500<=mom && mom<=9760) {
																			sumentries2++;
																			fluxene2->Fill(mom/1000, 1/cos(secondary_zenith_vec.at(i)));
																			fluxenecounts2->Fill(mom/1000);
																		}
																	}
													}										
											}
										}
									}
						
							}
					
						}
			
					infile->Close();

					}
				 }
				 
			else {  //stop if no file is remaining (list is completely analyzed)
				stop=true;
			}
		
		}
	}
	
	
}
