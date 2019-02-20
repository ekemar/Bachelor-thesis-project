#include <iostream>
#include <fstream>
#include <vector>

#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

namespace isotropic
{
double GetPi();

double randomUniformSinDisValue(double lowerBound, double upperBound, TRandom3 * random);

double randomUniformCosDisValue(double lowerBound, double upperBound, TRandom3 * random);

vector<double> CalcIntersection(double RE, double HD, double HA, double theta, double phi, double alpha, double beta);

vector<double> TransformPolar(vector <double> vec);

vector<double> CalcZenithAzimuth(double RE, double HD, double HA, double theta, double phi, double alpha, double beta);

void CalculateStartingPointBottomTop(double RE, double HD, double HA, 
				     double low_theta, double up_theta, 
				     double low_phi, double up_phi, 
				     int eventNumber, char * filename);
				     
void CompareDistributions(double RE, double H1, char * filename1, double H2, char * filename2);

vector<double> GetRandomStartingPosition(double RE, double HD, double HA, 
				     double low_theta, double up_theta, 
				     double low_phi, double up_phi, 
				     TRandom3 * random);
}		  
