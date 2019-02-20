#include <isotropic.hh>

double isotropic::GetPi() 
{
return 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491;
}

double isotropic::randomUniformSinDisValue(double lowerBound, double upperBound, TRandom3 * random)
{
lowerBound = sin(lowerBound * GetPi()/180);
upperBound = sin(upperBound * GetPi()/180);

return asin(random->Uniform(lowerBound, upperBound)) * 180/GetPi();
}

double isotropic::randomUniformCosDisValue(double lowerBound, double upperBound, TRandom3 * random)
{
lowerBound = cos(lowerBound * GetPi()/180);
upperBound = cos(upperBound * GetPi()/180);

return acos(random->Uniform(lowerBound, upperBound)) * 180/GetPi();
}

vector<double> isotropic::CalcIntersection(double RE, double HD, double HA, double theta, double phi, double alpha, double beta)
{
theta = theta * GetPi()/180;
phi = phi * GetPi()/180;
alpha = alpha * GetPi()/180;
beta = beta * GetPi()/180;

double intersection_x = (HD + RE) * cos(phi) * cos(theta) + (-0.1e1 * RE * cos(alpha) - 0.1e1 * HD * cos(alpha) + sqrt(RE * RE * pow(cos(alpha), 0.2e1) + 0.2e1 * RE * pow(cos(alpha), 0.2e1) * HD + HD * HD * pow(cos(alpha), 0.2e1) + HA * HA - 0.1e1 * HD * HD - 0.2e1 * HD * RE + 0.2e1 * RE * HA)) * (-0.1e1 * sin(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * cos(phi) * sin(theta) * cos(alpha) + cos(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1)))) / sin(theta);

double intersection_y = (HD + RE) * sin(phi) * cos(theta) + (-0.1e1 * HD * cos(alpha) - 0.1e1 * RE * cos(alpha) + sqrt(HD * HD * pow(cos(alpha), 0.2e1) + 0.2e1 * HD * pow(cos(alpha), 0.2e1) * RE + RE * RE * pow(cos(alpha), 0.2e1) + 0.2e1 * RE * HA - 0.2e1 * HD * RE + HA * HA - 0.1e1 * HD * HD)) * (sin(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(theta), 0.2e1)) * (-0.1e1 + pow(cos(alpha), 0.2e1))) + cos(theta) * sin(phi) * sin(theta) * cos(alpha) + cos(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(theta), 0.2e1)) * (-0.1e1 + pow(cos(alpha), 0.2e1)))) / sin(theta);

double intersection_z = (HD + RE) * sin(theta) - 0.1e1 * (-0.1e1 * HD * cos(alpha) - 0.1e1 * RE * cos(alpha) + sqrt(HD * HD * pow(cos(alpha), 0.2e1) + 0.2e1 * HD * pow(cos(alpha), 0.2e1) * RE + RE * RE * pow(cos(alpha), 0.2e1) + 0.2e1 * RE * HA - 0.2e1 * HD * RE + HA * HA - 0.1e1 * HD * HD)) * (-0.1e1 * cos(alpha) + cos(beta) * cos(theta) * sqrt((-0.1e1 + pow(cos(theta), 0.2e1)) * (-0.1e1 + pow(cos(alpha), 0.2e1))) + pow(cos(theta), 0.2e1) * cos(alpha)) / sin(theta);

vector<double> cartesic;
cartesic.push_back(intersection_x);
cartesic.push_back(intersection_y);
cartesic.push_back(intersection_z);

vector<double> polar = TransformPolar(cartesic);

return polar;
}

vector<double> isotropic::TransformPolar(vector <double> vec)
{
vector<double> result;

double R = sqrt(pow(vec[0], 0.2e1) + pow(vec[1], 0.2e1) + pow(vec[2], 0.2e1));

double Phi;
if(vec[1] < 0) Phi = 0.2e1 * GetPi() - acos(vec[0] * pow((pow(vec[0], 2) + pow(vec[1], 2)), -0.1e1 / 0.2e1)); 
else Phi = acos(vec[0] * pow((pow(vec[0], 2) + pow(vec[1], 2)), -0.1e1 / 0.2e1));

double Theta = atan(vec[2] * pow(vec[0] * vec[0] + vec[1] * vec[1], -0.1e1 / 0.2e1));

result.push_back(R);
result.push_back(Phi*180/GetPi());
result.push_back(Theta*180/GetPi());

return result;
}

vector<double> isotropic::CalcZenithAzimuth(double RE, double HD, double HA, double theta, double phi, double alpha, double beta)
{
vector<double> intersection = CalcIntersection(RE, HD, HA, theta, phi, alpha, beta);

theta = theta * GetPi()/180;
phi = phi * GetPi()/180;
alpha = alpha * GetPi()/180;
beta = beta * GetPi()/180;

double intersection_r = intersection[0];
double intersection_phi = GetPi()/180*intersection[1];
double intersection_theta = GetPi()/180*intersection[2];

vector<double> dir_rot(3);

dir_rot[0] = -0.1e1 * sin(intersection_theta) * cos(intersection_phi) * (-0.1e1 * sin(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * cos(phi) * sin(theta) * cos(alpha)) / sin(theta) - 0.1e1 * sin(intersection_theta) * sin(intersection_phi) * (cos(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + sin(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * sin(phi) * sin(theta) * cos(alpha)) / sin(theta) - 0.1e1 * cos(intersection_theta) * (pow(cos(theta), 0.2e1) * cos(alpha) + cos(beta) * cos(theta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) - 0.1e1 * cos(alpha)) / sin(theta);

dir_rot[1] = -0.1e1 * sin(intersection_phi) * (-0.1e1 * sin(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * cos(phi) * sin(theta) * cos(alpha)) / sin(theta) + cos(intersection_phi) * (cos(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + sin(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * sin(phi) * sin(theta) * cos(alpha)) / sin(theta);

dir_rot[2] = -0.1e1 * cos(intersection_theta) * cos(intersection_phi) * (-0.1e1 * sin(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * cos(phi) * sin(theta) * cos(alpha)) / sin(theta) - 0.1e1 * cos(intersection_theta) * sin(intersection_phi) * (cos(phi) * sin(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + sin(phi) * sin(theta) * cos(beta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) + cos(theta) * sin(phi) * sin(theta) * cos(alpha)) / sin(theta) + sin(intersection_theta) * (pow(cos(theta), 0.2e1) * cos(alpha) + cos(beta) * cos(theta) * sqrt((-0.1e1 + pow(cos(alpha), 0.2e1)) * (-0.1e1 + pow(cos(theta), 0.2e1))) - 0.1e1 * cos(alpha)) / sin(theta);

vector<double> r_az_zen = TransformPolar(dir_rot);

vector<double> result;

result.push_back(intersection_r);

result.push_back(intersection_phi*180/GetPi());

result.push_back(intersection_theta*180/GetPi());

double az;
if(180 - r_az_zen[1] < 0) az = 360 + 180 -  r_az_zen[1];
else az = 180 - r_az_zen[1];
if(result[2] < 0) 
	{
	if(az < 180) az += 180;
	else az -= 180;
	}
result.push_back(az);
	
result.push_back(90 + r_az_zen[2]);

return result;
}

void isotropic::CalculateStartingPointBottomTop(double RE, double HD, double HA, 
				     		double low_theta, double up_theta, 
				     		double low_phi, double up_phi, 
				     		int eventNumber, char * filename)
{
ofstream file(filename);
file.precision(20);

double downLat = low_theta;
double upLat = up_theta;

TRandom3 * random = new TRandom3();
random->SetSeed(0);

int ctr = 1;
while(ctr <= eventNumber) 
	{
	double randomLat = randomUniformSinDisValue(downLat, upLat, random);;
	
	double randomLong = random->Uniform(low_phi, up_phi);
	
	double randomBeta = random->Uniform(0,360);	  
	
	double randomAlpha = randomUniformCosDisValue(0, 90, random);
	
	vector<double> StartingPoint = CalcZenithAzimuth(RE, HD, HA, randomLat, randomLong, randomAlpha, randomBeta);
	
	vector <double> det_vector(3);
	det_vector[0] = (RE+HD)*cos(randomLong*GetPi()/180)*cos(randomLat*GetPi()/180);
	det_vector[1] = (RE+HD)*sin(randomLong*GetPi()/180)*cos(randomLat*GetPi()/180);
	det_vector[2] = (RE+HD)*sin(randomLat*GetPi()/180);
	
	vector <double> atmo_vector(3);
	atmo_vector[0] = (RE+HA)*cos(StartingPoint[1]*GetPi()/180)*cos(StartingPoint[2]*GetPi()/180);
	atmo_vector[1] = (RE+HA)*sin(StartingPoint[1]*GetPi()/180)*cos(StartingPoint[2]*GetPi()/180);
	atmo_vector[2] = (RE+HA)*sin(StartingPoint[2]*GetPi()/180);
	
	double distance = sqrt(pow(atmo_vector[0] - det_vector[0],2) + pow(atmo_vector[1] - det_vector[1],2)+ pow(atmo_vector[2] - det_vector[2],2));

	if(!isnan(StartingPoint[2]))	file<<StartingPoint[2]<<"\t"<<StartingPoint[1]<<"\t"<<StartingPoint[4]<<"\t"<<StartingPoint[3]<<"\t"<<StartingPoint[0]<<"\t"<<distance<<"\t"<<1<<endl;
		
	ctr++;	
	}
file.close();

delete random;			random = 0;
}

void isotropic::CompareDistributions(double RE, double H1, char * filename1, double H2, char * filename2)
{
ifstream file1(filename1);

ofstream file2(filename2);
file2.precision(20);

while(1)
	{
	if(!file1.good()) break;
	else 
		{
		double theta_1, phi_1, alpha_1, beta_1, r_1, distance_1, DetectorHit;
		
		file1>>theta_1>>phi_1>>alpha_1>>beta_1>>r_1>>distance_1>>DetectorHit;
		
		vector<double> StartingPoint = CalcZenithAzimuth(RE, H1, H2, theta_1, phi_1, alpha_1, beta_1);		
	
		vector <double> H1_vector(3);
		H1_vector[0] = (RE+H1)*cos(phi_1*GetPi()/180)*cos(theta_1*GetPi()/180);
		H1_vector[1] = (RE+H1)*sin(phi_1*GetPi()/180)*cos(theta_1*GetPi()/180);
		H1_vector[2] = (RE+H1)*sin(theta_1*GetPi()/180);
	
		vector <double> H2_vector(3);
		H2_vector[0] = (RE+H2)*cos(StartingPoint[1]*GetPi()/180)*cos(StartingPoint[2]*GetPi()/180);
		H2_vector[1] = (RE+H2)*sin(StartingPoint[1]*GetPi()/180)*cos(StartingPoint[2]*GetPi()/180);
		H2_vector[2] = (RE+H2)*sin(StartingPoint[2]*GetPi()/180);
	
		double distance = sqrt(pow(H2_vector[0] - H1_vector[0],2) + pow(H2_vector[1] - H1_vector[1],2)+ pow(H2_vector[2] - H1_vector[2],2));

		if(!isnan(StartingPoint[2])) file2<<StartingPoint[2]<<"\t"<<StartingPoint[1]<<"\t"<<StartingPoint[4]<<"\t"<<StartingPoint[3]<<"\t"<<StartingPoint[0]<<"\t"<<distance<<"\t"<<1<<endl;
		}
	}

file1.close();
file2.close();
}

vector<double> isotropic::GetRandomStartingPosition(double RE, double HD, double HA, 
				     		double low_theta, double up_theta, 
				     		double low_phi, double up_phi, 
				     		TRandom3 * random)
{
vector <double> lat_long_alpha_beta;

double randomLat = randomUniformSinDisValue(low_theta, up_theta, random);
double randomLong = random->Uniform(low_phi, up_phi);
if(randomLong > 360) randomLong -= 360;
double randomAlpha = randomUniformCosDisValue(0, 90, random);
double randomBeta = random->Uniform(0,360);	  
	
vector<double> StartingPoint = CalcZenithAzimuth(RE, HD, HA, randomLat, randomLong, randomAlpha, randomBeta);
	
lat_long_alpha_beta.push_back(StartingPoint[2]);
lat_long_alpha_beta.push_back(StartingPoint[1]);
lat_long_alpha_beta.push_back(StartingPoint[4]);
lat_long_alpha_beta.push_back(StartingPoint[3]);
	
return lat_long_alpha_beta;	
}


