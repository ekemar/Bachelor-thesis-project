#include "EarthMagneticField.hh"
#include "PlanetMagneticField.hh"
#include "EarthMagneticFieldMessenger.hh" 
#include "globals.hh"
#include "earth_magneto_fsubroutine_def.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4FieldManager.hh"
#include "SpaceCoordinatePlanet.hh"
#include "PlanetUnits.hh"
#include <time.h>

////////////////////////////////////////////////////////////////////////////////
//
EarthMagneticField::EarthMagneticField() : PlanetMagneticField("Earth")
{ 
  ListOfInternalFieldModels.push_back("IGRF");
  ListOfInternalFieldModels.push_back("IGRFTiltedDipole");
  ListOfInternalFieldModels.push_back("IGRFEccentricTiltedDipole");
  ListOfExternalFieldModels.push_back("TSY89");
  ListOfExternalFieldModels.push_back("TSY96");
  ListOfExternalFieldModels.push_back("TSY2001");
  ListOfExternalFieldModels.push_back("TSY2004");
  ListOfMagnetopauseModels.push_back("TSY89");
  ListOfMagnetopauseModels.push_back("TSY96");
  ListOfMagnetopauseModels.push_back("TSY2001");
  ListOfMagnetopauseModels.push_back("TSY2004");
  
  HasAGlobalField =true;
 
 //messenger
  theFieldMessenger = new EarthMagneticFieldMessenger(this); 
  
 //coefficients for igrf
  hh=new float[105];
  gg=new float[105];
  rec=new float[105];

  //Igrf table
  if (getenv("IGRF_TABLE")){
  	ReadIgrfTable(getenv("IGRF_TABLE"));
  }
  else {
   G4cout<<"The environment variable IGRF_TABLE has not been defined"<<std::endl;
  }

 TiltAngle = -1000;
  
 //Tsyganenko Parameters
 
  Pdyn=2.;
  Dst=0.;
  Imfy=0.;
  Imfz=1.*nanotesla;
  G1=1.;
  G2=0.;
  nm_gauss=13;
  nm_gauss=13;
  
  W1=0;
  W2=0;
  W3=0;
  W4=0;
  W5=0;
  W6=0;
  
  times_of_data.clear();
  Pdyn_data.clear();
  Dst_data.clear();
  Imf_gsm_y_data.clear();
  Imf_gsm_z_data.clear();
  g1_tsy_param_data.clear();
  g2_tsy_param_data.clear();
  w1_tsy_param_data.clear();
  w2_tsy_param_data.clear();
  w3_tsy_param_data.clear();
  w4_tsy_param_data.clear();
  w5_tsy_param_data.clear();
  w6_tsy_param_data.clear();
  n_tsy_data=0;
  iopt =1;

  //Field model  
     
  SetInternalField("IGRF");
  SetExternalField("NOFIELD");
  SetMagnetopauseModel("SPHERE");
  RadiusMagnetosphere = 25.* Re;
 
  External=false;
  Internal=true;
  
  iopt=1;
  
  ConsiderDipoleShift=false;  
  
 //Set the start date and compute IGRF coefficient 

  StartDate = DateAndTime(2000,1,1);
  SetTimeOfB(0.);
}
////////////////////////////////////////////////////////////////////////
//
EarthMagneticField::~EarthMagneticField()
{ 
}			       
/////////////////////////////////////////////////////////////////////////////////
//
/*G4ThreeVector EarthMagneticField::FindPositionOnMagnetopause(G4double theta, 
                                                              G4double xgsm, 
							      G4double precision) const
{ G4ThreeVector pos1 = G4ThreeVector(xgsm,0.,0.);
 
  SpaceCoordinatePlanet* theCoordinateConvertor
                            = SpaceCoordinatePlanet::GetInstance();
  theCoordinateConvertor->SetSystemInAndOut("GSM","GEO");
 
  if  (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos1)))
                                            return G4ThreeVector (0.,0.,0.);
 
  G4ThreeVector pos2 =G4ThreeVector(xgsm,40.*re*std::sin(theta),40.*re*std::cos(theta));
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
 	pos3= (pos1+pos2)/2.;
	if (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos3))) pos2=pos3;
	else pos1=pos3;
  }
  return (pos1+pos2)/2.;
}
*/
////////////////////////////////////////////////////////////////////////////////
//
/*G4ThreeVector EarthMagneticField::FindStandOffPosition(G4double precision ) const
{ G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.)*Re;
  G4ThreeVector pos2 = G4ThreeVector(40.,0.,0.)*Re;
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
  	pos3= (pos1+pos2)/2.;
	if (OutsideMagnetosphere(pos3)) pos2=pos3;
	else pos1=pos3;/usr/lib/gcc-lib/i586-suse-linux/3.3.3/specs
  }
  return pos1;	
}
*/
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::SetIopt(  G4int val)
{
  iopt=val;
  if (val <1) iopt=1;
  if (val >7) iopt=7;
  return ;
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::SetInternalField(G4String aString)
{  
   if (aString ==  "IGRFTiltedDipole"){
   	
	SetTiltedDipoleParameterFromIGRF();
	SetMotherInternalField("DIPOLE");
	PInternalField= PMotherInternalField;
	return;
	
   }
   else if (aString ==  "IGRFEccentricTiltedDipole"){
   	
	SetEccentricDipoleParameterFromIGRF();
	SetMotherInternalField("DIPOLE");
	PInternalField= PMotherInternalField;
	return;
	
   }
   else if (!SetMotherInternalField(aString)){
	Internal=true;   
  	if (aString == "IGRFC++")
  		PInternalField=&EarthMagneticField::GetIGRF;
  	else if (aString == "IGRFC++DOUBLE")
  		PInternalField=&EarthMagneticField::GetIGRFDOUBLE;	 
  	else if (aString == "IGRFFORTRAN")
  		PInternalField=&EarthMagneticField::GetIGRFFortran;
  	else if (aString == "IGRFFORTRAN1")
  		PInternalField=&EarthMagneticField::GetIGRFFortran1;
	
	else if (aString == "IGRF"){
  		PInternalField=&EarthMagneticField::GetIGRFFortran1;
		
	}
  	else Internal=false;  
  }
  else {
  	PInternalField= PMotherInternalField;
  }	   
}
///////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::SetExternalField(  G4String aString)
{External=true;
 if (aString == "TSY89"){
  	PExternalField=&EarthMagneticField::GetTSY89;
  	SetMagnetopauseModel("TSY89");
 }
 else if (aString == "TSYBOB89"){ 
  	PExternalField=&EarthMagneticField::GetTSYBOB89;
  	SetMagnetopauseModel("TSY89");
 } 
 else if (aString == "TSY96"){
  	PExternalField=&EarthMagneticField::GetTSY96;
  	SetMagnetopauseModel("TSY96");
 }
 else if (aString == "TSY2001"){
  	PExternalField=&EarthMagneticField::GetTSY2001;
  	SetMagnetopauseModel("TSY2001");
 }
 else if (aString == "TSY2004"){
  	PExternalField=&EarthMagneticField::GetTSY2004;
  	SetMagnetopauseModel("TSY2004");
 } 
 else  {External=false;
 	SetMagnetopauseModel("SPHERE");
	
 }	
  
}

////////////////////////////////////////////////////////////////////////////////
//
void  EarthMagneticField::SetMagnetopauseModel( G4String aString) 
{ if (!SetMotherMagnetopauseModel( aString)){
  	Magnetopause =true;
	if (aString == "TSY2004")
 		POutsideMagnetosphere
                  	=&EarthMagneticField::TSY2004OutsideMagnetosphere;
			
	if (aString == "TSY2001")
 		POutsideMagnetosphere
                  	=&EarthMagneticField::TSY2001OutsideMagnetosphere;
		  
 	else if (aString == "TSY96")
     		POutsideMagnetosphere
                  	=&EarthMagneticField::TSY96OutsideMagnetosphere; 		  
		  
 	else  if (aString == "TSY89")
     		POutsideMagnetosphere
                  	=&EarthMagneticField::TSY89OutsideMagnetosphere; 
	    
 	else	Magnetopause =false;
 		
  }
  else 	POutsideMagnetosphere = PMotherOutsideMagnetosphere;	   	    

}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetInternalField(G4ThreeVector aVector) const
{ return (this->*PInternalField)(aVector);  
}
///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetExternalField(G4ThreeVector aVector) const
{ return (this->*PExternalField)(aVector); 
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ComputeBfieldParametersAtReferenceTime()
{ ComputeIgrfCoefAccordingToTime();
  ComputeTSY2001Parameter(TimeOfB);



}

////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ReadTSY2001Parameter(G4String nameFile)
{  G4cout<<"ok"<<std::endl;
  times_of_data.clear();
  Pdyn_data.clear();
  Dst_data.clear();
  Imf_gsm_y_data.clear();
  Imf_gsm_z_data.clear();
  g1_tsy_param_data.clear();
  g2_tsy_param_data.clear();
  n_tsy_data=0;
  G4cout<<"ok1"<<std::endl;
  std::fstream File_Input(nameFile, (std::ios::binary | std::ios::in));
  File_Input>>n_tsy_data;
  G4cout<<"ok2"<<std::endl;
 G4cout<<n_tsy_data<<std::endl;
  if (n_tsy_data <= 0){
  	G4cout<<"the file does not exist or is not well written"<<G4endl;
  } 
  else{ G4int nyear;
   	G4int nmonth;
   	G4int nday;
   	G4int nhour;
   	G4int nminute;
   	G4int nsecond=0;
   	File_Input>>nyear>>nmonth>>nday>>nhour>>nminute;
	G4double a,b,c,d,e,f,g;
 	G4cout<<"ok3"<<std::endl;
   	for (G4int i=0;i<n_tsy_data;i++){
   		File_Input>>a>>b>>c>>d>>e>>f>>g;
		G4cout<<a<<std::endl;
		times_of_data.push_back(a);
	        Pdyn_data.push_back(b);
              	Dst_data.push_back(c);
              	Imf_gsm_y_data.push_back(d);
                Imf_gsm_z_data.push_back(e);
              	g1_tsy_param_data.push_back(f);
              	g2_tsy_param_data.push_back(g);
   	}
	SetStartDate(nyear, nmonth, nday, nhour, nminute, nsecond);

  } 
  File_Input.close(); 
  double val =TimeOfB;
  SetTimeOfB(val);
}

////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ReadTSY2004Parameter(G4String nameFile)
{ //to be done
	
	/* G4cout<<"ok"<<std::endl;
  times_of_data.clear();
  Pdyn_data.clear();
  Dst_data.clear();
  Imf_gsm_y_data.clear();
  Imf_gsm_z_data.clear();
  g1_tsy_param_data.clear();
  g2_tsy_param_data.clear();
  n_tsy_data=0;
  G4cout<<"ok1"<<std::endl;
  std::fstream File_Input(nameFile, (std::ios::binary | std::ios::in));
  File_Input>>n_tsy_data;
  G4cout<<"ok2"<<std::endl;
 G4cout<<n_tsy_data<<std::endl;
  if (n_tsy_data <= 0){
  	G4cout<<"the file does not exist or is not well written"<<G4endl;
  } 
  else{ G4int nyear;
   	G4int nmonth;
   	G4int nday;
   	G4int nhour;
   	G4int nminute;
   	G4int nsecond=0;
   	File_Input>>nyear>>nmonth>>nday>>nhour>>nminute;
	G4double a,b,c,d,e,f,g;
 	G4cout<<"ok3"<<std::endl;
   	for (G4int i=0;i<n_tsy_data;i++){
   		File_Input>>a>>b>>c>>d>>e>>f>>g;
		G4cout<<a<<std::endl;
		times_of_data.push_back(a);
	        Pdyn_data.push_back(b);
              	Dst_data.push_back(c);
              	Imf_gsm_y_data.push_back(d);
                Imf_gsm_z_data.push_back(e);
              	g1_tsy_param_data.push_back(f);
              	g2_tsy_param_data.push_back(g);
   	}
	SetStartDate(nyear, nmonth, nday, nhour, nminute, nsecond);

  } 
  File_Input.close(); 
  double val =TimeOfB;
  SetTimeOfB(val);*/
}

////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ReadIgrfTable(G4String nameFile)
{ n_tsy_data=0;
  std::fstream File_Input(nameFile,  std::ios::in);
 //clear vectors
  igrf_year.clear();
  h_nm.clear();
  g_nm.clear();
  dh_nm.clear();
  dg_nm.clear();
  hh_nm.clear();
  gg_nm.clear();
  hh_nm1.clear();
  gg_nm1.clear();
  
  recurrence_coef.clear();
 
 

 //Process the first line 
  char ch;
  std::stringstream astream;
  File_Input.get(ch);
  while (ch != '\n') {
  	astream<<ch;
   	File_Input.get(ch);
  }
  G4String str1;
  astream>>str1>>str1>>str1;
  if (str1!="m"){
  	G4cout<<"The environment variable IGRF_TABLE is not defined correctly" <<std::endl;
  	G4cout<<"The program will be interrupted"<<std::endl;
	exit(0);
	
  }
  G4double year=0.;
  while (!astream.eof() && year >= 0.){
 	year=-999999.;
  	astream>>year;
  	if (year > 0){ 
     		igrf_year.push_back(year);
     	}
  }

  //number of column
 
  G4int nyear=igrf_year.size();
  nmax=0;
  G4int n,m;

  // initialisation h_nm,g_nm;  
  
  for (int i=0; i <nyear ; i++){
   	h_nm.push_back(std::vector<float> ());
    	g_nm.push_back(std::vector<float> ());
  }

  while (!File_Input.eof()){
  	str1="bibab";
    	n=-999;
    	m=-999;
    	File_Input>>str1>>n>>m;

    	if (str1 != "bibab" && n != -999 && m != -999) {

     		if (n>nmax){
      			nmax=n;
       			for (int i=0; i <nyear ; i++){
        			h_nm[i].insert(h_nm[i].end(),nmax+1,0.);
	 			g_nm[i].insert(g_nm[i].end(),nmax+1,0.);
			}
       			dg_nm.insert(dg_nm.end(),nmax+1,0.);
       			dh_nm.insert(dh_nm.end(),nmax+1,0.);
       			hh_nm.insert(hh_nm.end(),nmax+1,0.);
       			gg_nm.insert(gg_nm.end(),nmax+1,0.);
			hh_nm1.insert(hh_nm1.end(),nmax+1,0.);
       			gg_nm1.insert(gg_nm1.end(),nmax+1,0.);	    
      		}
     		G4int index = ((n+2)  * (n-1)) /2 + m;

     		G4double value;
     		int i_start=0;
     		if (n >10) {
     			i_start = 20;
			for (int i=0;i<i_start;i++){	
				h_nm[i][index]  = 0.;
				g_nm[i][index]  = 0.;
			}
     		}				
     
     		for (int i= i_start; i <nyear ; i++){
      			File_Input>>value;

       			if (str1 == "h") h_nm[i][index] = value;
       			else if (str1 == "g") g_nm[i][index] = value;
      		}
     		if (n<9) {//cautiona the dh and dg parameter are pro year 
     			File_Input>>value;
     			if (str1 == "h") dh_nm[index] = value;
     			else if (str1 == "g") dg_nm[index] = value;
     			
     		}
     		else {
     			dh_nm[index] = 0.;
			dg_nm[index] = 0.;
     		}		
    	}
  }
 
 //Here we calculate the coefficient for the recurrence equation used to compute 
 //the associated Legendre function in the Gauss normalsation. 
 // In this normalisation the equation of recurrence is
 // P_n_m(x) = x * p_n-1_m -(n+m)*(n-m)* p_n-2_m / ((2n-3)*(2n-1)) (1)
 // The coefficient that we calculate here is (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1))
 // The qeaution (1) is a transformation of the classical  equation of recurrence of lengendre assocaited function
 // (n-m) * P_n_m = x  * p_n-1_m *(2n-1) -(n + m -1) * p_m_n-2  (2)
 // by using the Gauss normalisation (2) transform in (1)
 // (1) is better for computing purpose because you have just to compute 
 //the coefficient  (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) at initialisation an multiplied it by
 //   p_m_n-2  when computinhg p_n_m
 // In  case of (2) you have to compute three vector of cooefficient at initilaisation you have to
 // make three multiplications instead of one 


   for (n=1;n<nmax + 1;n++)
  	for (int m= 0; m < n+1; m++) 
   		recurrence_coef.push_back(
       			 float ((n-1+m)*(n-1-m)) / float ((2*n-1)*(2*n-3)) );
   
   for (n=1;n<=nmax + 1;n++)
  	for (int m= 1; m < n+1; m++){
		int mn=n*(n-1)/2+m;
   		igrfcc_.rec[mn-1] = float ((n-m)*(n+m-2)) / float ((2*n-1)*(2*n-3));
	
	};
			 		 
 
 
 File_Input.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ComputeIgrfCoefAccordingToTime()
{
  G4int nyear = igrf_year.size();
  G4int n=(nmax+3)*nmax/2;
  
  //G4cout<<"number of year is "<<nyear<<endl;
 
  if (ReferenceDate.year < igrf_year[0]  || 
                     ReferenceDate.year> igrf_year[nyear-1] + 5 ){
     	G4cout<<"The year should be in the period "<<int(igrf_year[0])
                                   <<"-"<<int(igrf_year[nyear-1]) + 5<<std::endl;
      
      	G4cout<<"The coefficient will be set according to year "<<
       	                                   int(igrf_year[nyear-1])<<std::endl;
      	for (int i=0;i<n;i++){
       		gg_nm[i]=g_nm[nyear-1][i];
        	hh_nm[i]=h_nm[nyear-1][i];
		gg_nm1[i]=g_nm[nyear-1][i];
        	hh_nm1[i]=h_nm[nyear-1][i];}
  }
  else {
     	G4int index=0;
     
      	while (ReferenceDate.year >= (igrf_year)[index] && index < nyear) index++;
      
     	// G4cout <<"Selected index"<<index<<std::endl;
      	if (index<nyear){
        	DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
	 	DateAndTime Date2=DateAndTime(int (igrf_year[index]), 1,1);
	 	G4double delta_day12=Date2.DifferenceInDays(Date1);
	 	G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
	 	G4double f2= delta_day / delta_day12;
         	G4double f1=1.-f2;
       		for (int i=0;i<n;i++){
          		gg_nm[i]=g_nm[index-1][i]*f1 + g_nm[index][i]*f2;
           		hh_nm[i]=h_nm[index-1][i]*f1 + h_nm[index][i]*f2;
			hh_nm1[i]=hh_nm[i];
			gg_nm1[i]=gg_nm[i];
		}
      	}
      	else {
        	DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
	 	G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
	 	G4double dt=delta_day / 365.25;
         	for (int i=0;i<n;i++){
          		gg_nm[i]=g_nm[nyear-1][i] + dt*dg_nm[i];
           		hh_nm[i]=h_nm[nyear-1][i] + dt*dh_nm[i];
			hh_nm1[i]=hh_nm[i];
			gg_nm1[i]=gg_nm[i];
		}
        }
  }
     
  //The coefficient are  Schmidt quasi normalised
  // For using the recuurence relation in Gauss normalisation
  // we have to convert the coefficient in Gauss normalisation
  // the g_nm and h_nm coefficient become
  float s=1.;
  for (n=1;n<nmax+1;n++){
  	G4int mn = ((n + 1) * n) / 2  -1;
       	s*=float(2*n-1)/float(n);
       	hh_nm[mn] *= s;
       	gg_nm[mn] *= s;
       	G4double p=s;
       	for (int m=1;m<n+1;m++){
        	float aa=1.;
	 	if (m == 1) aa=2.;
	 	p *= std::sqrt(aa * float(n - m + 1) / float (n + m) );
	 	int nmn =  mn + m ;
	 	hh_nm[nmn] *= p;
         	gg_nm[nmn] *= p;
		hh_nm1[nmn]=hh_nm[nmn];
		gg_nm1[nmn]=gg_nm[nmn];
		
	}
  }
  igrfcc_.hh[0]=0.;
  igrfcc_.gg[0]=0.;
  for (unsigned int i=0;i<hh_nm.size();i++){
  	hh[i]=hh_nm[i];
	gg[i]=gg_nm[i];
	igrfcc_.hh[i+1]=hh_nm[i];
	igrfcc_.gg[i+1]=gg_nm[i];
  } 				   			     
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::SetTiltedDipoleParameterFromIGRF()
{ComputeIgrfCoefAccordingToTime();
 SetTiltedDipoleParameterFromGAUSSCoefficients();
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::SetEccentricDipoleParameterFromIGRF()
{ComputeIgrfCoefAccordingToTime();
 SetEccentredDipoleParameterFromGAUSSCoefficients();
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::PrintStormParameterTSY2001()
{ G4cout<<"iopt "<<iopt<<G4endl;
  G4cout<<"Pdyn "<<Pdyn<<G4endl;
  G4cout<<"TiltAngle "<<TiltAngle<<G4endl;
  G4cout<<"Dst "<<Dst/nanotesla<<G4endl;
  G4cout<<"Imfy "<<Imfy/nanotesla<<G4endl;
  G4cout<<"Imfz "<<Imfz/nanotesla<<G4endl;
  G4cout<<"G1 "<<G1<<G4endl;
  G4cout<<"G2 "<<G2<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::PrintStormParameterTSY2004()
{
  G4cout<<"iopt "<<iopt<<G4endl;
  G4cout<<"Pdyn "<<Pdyn<<G4endl;
  G4cout<<"TiltAngle "<<TiltAngle<<G4endl;
  G4cout<<"Dst "<<Dst/nanotesla<<G4endl;
  G4cout<<"Imfy "<<Imfy/nanotesla<<G4endl;
  G4cout<<"Imfz "<<Imfz/nanotesla<<G4endl;
  G4cout<<"W1 "<<W1<<G4endl;
  G4cout<<"W2 "<<W2<<G4endl;
  G4cout<<"W3 "<<W3<<G4endl;
  G4cout<<"W4 "<<W4<<G4endl;
  G4cout<<"W5 "<<W5<<G4endl;
  G4cout<<"W6 "<<W6<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetIGRFFortran(G4ThreeVector pos) const 
{ G4double r=pos.mag()/re;
  G4double t=pos.theta();
  G4double p=pos.phi();
  float Br=1.;
  float Bp=1.;
  float Bt=1.;
  G4int nm=nm_gauss;
  //year_igrf=1980;
   
 /* G4double ndays_per_y  = DateAndTime(ReferenceDate.year+1,1,1)
                                  .DifferenceInDays(DateAndTime(ReferenceDate.year,1,1));
  */ 
 /* G4double fraction_of_year =
                     ReferenceDate.DifferenceInDays(
		                    DateAndTime(ReferenceDate.year,1,1))/ndays_per_y;
  float yea=float(ReferenceDate.year + fraction_of_year)  ;*/
  float r1 =float(r);
  float t1 =float(t);
  float p1 =float(p);

   
  if (nm > 10) nm=10;
  
  
  igrf_to_cc_(&r1,&t1,&p1,&Br,&Bt,&Bp); 
   
  G4double Bunit=nanotesla;
  G4double st=sin(t);
  G4double ct=cos(t);
  G4double sp=sin(p);
  G4double cp=cos(p);
  G4double Be=double(Br)*st+double(Bt)*ct;
  G4double Bx=Bunit*(Be*cp-double(Bp)*sp);
  G4double By=Bunit*(Be*sp+double(Bp)*cp);
  G4double Bz=Bunit*(double(Br)*ct-double(Bt)*st);
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetIGRFFortran1(G4ThreeVector pos) const 
{ G4double r=pos.mag()/rplanet;
  G4double t=pos.theta();
  G4double p=pos.phi();
  float Br=1.;
  float Bp=1.;
  float Bt=1.;
  G4int old_nm_gauss =nm_gauss;
  /*G4cout<<nm_gauss<<std::endl;
  G4cout<<pos/rplanet<<std::endl;*/
  if (ReferenceDate.year <2000. && nm_gauss >10) old_nm_gauss =10; 
  igrfcc_.nm=old_nm_gauss;
 
   
 
  
 
  float r1 =float(r);
  float t1 =float(t);
  float p1 =float(p);
  
  
  igrf_to_cc_(&r1,&t1,&p1,&Br,&Bt,&Bp);
   
   
  G4double Bunit=nanotesla;
  G4double st=sin(t);
  G4double ct=cos(t);
  G4double sp=sin(p);
  G4double cp=cos(p);
  G4double Be=double(Br)*st+double(Bt)*ct;
  G4double Bx=Bunit*(Be*cp-double(Bp)*sp);
  G4double By=Bunit*(Be*sp+double(Bp)*cp);
  G4double Bz=Bunit*(double(Br)*ct-double(Bt)*st);
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetIGRF(G4ThreeVector pos) const 
{ 
  // This method computes the igrf magnetic field
  // For this purpose associated legendre function are calculated by the recurrence equation
  // P_n_m(x) = x * p_n-1_m(x) - recurrence_coef(n,m) * p_n-2_m(x)  (1)
  // where recurrence_coef = (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) and x= cos_theta
  // By deriving  equation (1) by theta we obtain
  // dp_n_m/dtheta = -sin_theta * p_n-1_m + cos_theta * dp_n-1_m/dtheta  -  
  //                                                        recurrence_coef * dp_n-2_m/dtheta
  // and dp_n_m/dtheta can  also be computed by recurrence. 
  // we should also now p_m_m = (sin (theta))^m


  float r=float(pos.mag()/re);
  float theta=float(pos.theta());
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  float phi=float(pos.phi());
   
  float pp = float(1.) / r;
  float p =  pp * pp;
 
   
   
   
  float cos_phi = std::cos(phi);
  float two_cos_phi = std::cos(phi) * 2.;
  float sin_phi = std::sin(phi);
  float cos_theta = std::cos(theta);
  float sin_theta = std::sin(theta);
   
   
   
  // a_n = (1 /r) ^ (n+1)
  // da_n_dth = -(a_n1) * n+1 
 
  std::vector<float> a;
  std::vector<float> da_dr;
  G4int n;
  G4int old_nm_gauss =nm_gauss;
  if (ReferenceDate.year <2000. && nm_gauss >10) old_nm_gauss =10; 	
   
   
  for ( n=1;n<old_nm_gauss+1;n++){
    	p *= pp;
     	a.push_back(p);
     	da_dr.push_back( -p * float (n+1) );
  }
    
  // case where m=0, p_1_0 =1.
  //---------------------------
    
  // n=1 
  float p_nm=cos_theta;
  float dp_nm_dth=-sin_theta;
    
  float Br=-da_dr[0]*gg_nm[0]*p_nm;
  float Bt=-a[0]*gg_nm[0]*dp_nm_dth;
  float Bp=0.;
    
  float p_nm1=p_nm; //corresponds to P for n-1 and m
  float p_nm2=1.;   //corresponds to P for n-2 and m
  
  float dp_nm1_dth=dp_nm_dth;
  float dp_nm2_dth=0.;
    
  // n>2
    
 for (n = 2 ; n < old_nm_gauss +1 ; n++ ){
    	G4int ii= ((n+2) * (n-1) )/2; 
      	p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
      	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                       - recurrence_coef[ii]*dp_nm2_dth;
      	Br-=da_dr[n-1]*gg_nm[ii]*p_nm;
      	Bt-=a[n-1]*gg_nm[ii]*dp_nm_dth; 
      	p_nm2=p_nm1;
      	p_nm1=p_nm;
      	dp_nm2_dth=dp_nm1_dth;
      	dp_nm1_dth=dp_nm_dth;
  } 
    
  // continue for m>0
  //---------  
    
  float cos_mphi2=cos_phi;
  float cos_mphi1=1.;
  float cos_mphi=cos_phi;
  float sin_mphi2=-sin_phi;
  float sin_mphi1=0.;
  float sin_mphi=sin_phi;
    
    
  float p_mm=1.;
  float dp_mm_dth=1.; 
  float Bpm=0.;
  float Bp_hh,Bp_gg,Br_hh,Br_gg,Bt_hh,Bt_gg;
  for (int m = 1 ; m < old_nm_gauss +1 ; m++ ){
    	//Initialisation
     	//------------------
     
       	Bp_hh=0.;
       	Bp_gg=0.;
       	Br_hh=0.;
       	Br_gg=0.;
       	Bt_hh=0.;
       	Bt_gg=0.;
      
     	//compute cos_mphi and sin_mphi with
     	// recurrence formula from numerical recipies in C++ p 184 
      	cos_mphi= two_cos_phi*cos_mphi1 -cos_mphi2;
      	sin_mphi= two_cos_phi*sin_mphi1 -sin_mphi2;
      
      
      	cos_mphi2=cos_mphi1;
      	sin_mphi2=sin_mphi1;
      	cos_mphi1=cos_mphi;
      	sin_mphi1=sin_mphi;
    
      	dp_mm_dth= float(m) * p_mm *cos_theta;
      	p_mm *= sin_theta;
      	Bpm=0;
     	// case where n=m
      	n=m;
      	G4int ii= ((n+3) * (n) )/2 -1; 
      	p_nm = p_mm;
      	dp_nm_dth = dp_mm_dth;
	
	float gg=gg_nm[ii];
	float hh=hh_nm[ii];
	float aa=a[n-1];
	float da=da_dr[n-1]; 
      
      	float ggg = gg*p_nm;
      	float hhh = hh*p_nm;
      
      	Br_gg-=da*ggg;
      	Br_hh-=da*hhh;
        
	float w=aa*dp_nm_dth;
      	Bt_gg-=gg*w;
      	Bt_hh-=hh*w;
      
      	//Bt-=a[n-1]*gh*dp_nm_dth;
      
      	//Bpm-=a[n-1]*dgh_dphi*p_nm; 
      
      	Bp_gg-=aa*ggg;
      	Bp_hh-=aa*hhh; 
      
      
      	p_nm1=p_nm; //corresponds to P for n-1 and m
      	p_nm2=0.;  //corresponds to P for n-2 and m
      	dp_nm1_dth=dp_nm_dth;
     	 dp_nm2_dth=0.;
        
     
        for ( n = m + 1 ; n < old_nm_gauss +1 ; n++ ){
        	G4int ii= ((n+2) * (n-1) )/2 + m;
		float rc= recurrence_coef[ii];
         	p_nm = cos_theta *p_nm1 - rc*p_nm2;
         	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                                 - rc*dp_nm2_dth;
	 
		
		gg=gg_nm[ii];
		hh=hh_nm[ii];
	 
	 
	 	float ggg = gg*p_nm;
         	float hhh = hh*p_nm;
		da=da_dr[n-1];
		aa=a[n-1];
      
         	Br_gg-=da*ggg;
        	Br_hh-=da*hhh;
                
		w=aa*dp_nm_dth;
         	Bt_gg-=w*gg;
         	Bt_hh-=w*hh;
      
         	Bp_gg-=aa*ggg;
         	Bp_hh-=aa*hhh; 

	 	p_nm2=p_nm1;
	 	p_nm1=p_nm;
	 	dp_nm2_dth=dp_nm1_dth;
	 	dp_nm1_dth=dp_nm_dth;
	}
	Br+=Br_gg*cos_mphi + Br_hh*sin_mphi;
	Bt+=Bt_gg*cos_mphi + Bt_hh*sin_mphi;
	Bp+=float(m)*(Bp_hh*cos_mphi - Bp_gg*sin_mphi);

  } 
  
  if (sin_theta >= 1.e-8) Bp/=sin_theta; 
  else if (cos_theta<0.) Bp=-Bp;
   
  G4double Be= Br*sin_theta+ Bt*cos_theta;
  G4double Bx= Be*cos_phi- Bp*sin_phi;
  G4double By= Be*sin_phi+ Bp*cos_phi;
  G4double Bz= Br*cos_theta - Bt*sin_theta;
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
    
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz)*nanotesla;
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetIGRFDOUBLE(G4ThreeVector pos) const 
{ 
  // This method computes the igrf magnetic field
  // For this purpose associated legendre function are calculated by the recurrence equation
  // P_n_m(x) = x * p_n-1_m(x) - recurrence_coef(n,m) * p_n-2_m(x)  (1)
  // where recurrence_coef = (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) and x= cos_theta
  // By deriving  equation (1) by theta we obtain
  // dp_n_m/dtheta = -sin_theta * p_n-1_m + cos_theta * dp_n-1_m/dtheta  -  
  //                                                        recurrence_coef * dp_n-2_m/dtheta
  // and dp_n_m/dtheta can  also be computed by recurrence. 
  // we should also now p_m_m = (sin (theta))^m


  double r=double(pos.mag()/re);
  double theta=double(pos.theta());
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  double phi=double(pos.phi());
   
  double pp = double(1.) / r;
  double p =  pp * pp;
 
   
   
   
  double cos_phi = std::cos(phi);
  double two_cos_phi = std::cos(phi) * 2.;
  double sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);
   
   
   
  // a_n = (1 /r) ^ (n+1)
  // da_n_dth = -(a_n1) * n+1 
 
  std::vector<double> a;
  std::vector<double> da_dr;
  G4int n;
  G4int old_nm_gauss =nm_gauss;
  if (ReferenceDate.year <2000. && nm_gauss >10) old_nm_gauss =10; 	
   
   
  for ( n=1;n<old_nm_gauss+1;n++){
    	p *= pp;
     	a.push_back(p);
     	da_dr.push_back( -p * double (n+1) );
  }
    
  // case where m=0, p_1_0 =1.
  //---------------------------
    
  // n=1 
  double p_nm=cos_theta;
  double dp_nm_dth=-sin_theta;
    
  double Br=-da_dr[0]*gg_nm1[0]*p_nm;
  double Bt=-a[0]*gg_nm1[0]*dp_nm_dth;
  double Bp=0.;
    
  double p_nm1=p_nm; //corresponds to P for n-1 and m
  double p_nm2=1.;   //corresponds to P for n-2 and m
  
  double dp_nm1_dth=dp_nm_dth;
  double dp_nm2_dth=0.;
    
  // n>2
    
 for (n = 2 ; n < old_nm_gauss +1 ; n++ ){
    	G4int ii= ((n+2) * (n-1) )/2; 
      	p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
      	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                       - recurrence_coef[ii]*dp_nm2_dth;
      	Br-=da_dr[n-1]*gg_nm1[ii]*p_nm;
      	Bt-=a[n-1]*gg_nm1[ii]*dp_nm_dth; 
      	p_nm2=p_nm1;
      	p_nm1=p_nm;
      	dp_nm2_dth=dp_nm1_dth;
      	dp_nm1_dth=dp_nm_dth;
  } 
    
  // continue for m>0
  //---------  
    
  double cos_mphi2=cos_phi;
  double cos_mphi1=1.;
  double cos_mphi=cos_phi;
  double sin_mphi2=-sin_phi;
  double sin_mphi1=0.;
  double sin_mphi=sin_phi;
    
    
  double p_mm=1.;
  double dp_mm_dth=1.; 
  double Bpm=0.;
  double Bp_hh,Bp_gg,Br_hh,Br_gg,Bt_hh,Bt_gg;
  for (int m = 1 ; m < old_nm_gauss +1 ; m++ ){
    	//Initialisation
     	//------------------
     
       	Bp_hh=0.;
       	Bp_gg=0.;
       	Br_hh=0.;
       	Br_gg=0.;
       	Bt_hh=0.;
       	Bt_gg=0.;
      
     	//compute cos_mphi and sin_mphi with
     	// recurrence formula from numerical recipies in C++ p 184 
      	cos_mphi= two_cos_phi*cos_mphi1 -cos_mphi2;
      	sin_mphi= two_cos_phi*sin_mphi1 -sin_mphi2;
      
      
      	cos_mphi2=cos_mphi1;
      	sin_mphi2=sin_mphi1;
      	cos_mphi1=cos_mphi;
      	sin_mphi1=sin_mphi;
    
      	dp_mm_dth= double(m) * p_mm *cos_theta;
      	p_mm *= sin_theta;
      	Bpm=0;
     	// case where n=m
      	n=m;
      	G4int ii= ((n+3) * (n) )/2 -1; 
      	p_nm = p_mm;
      	dp_nm_dth = dp_mm_dth;
	
	double gg=gg_nm1[ii];
	double hh=hh_nm1[ii];
	double aa=a[n-1];
	double da=da_dr[n-1]; 
      
      	double ggg = gg*p_nm;
      	double hhh = hh*p_nm;
      
      	Br_gg-=da*ggg;
      	Br_hh-=da*hhh;
        
	double w=aa*dp_nm_dth;
      	Bt_gg-=gg*w;
      	Bt_hh-=hh*w;
      
      	//Bt-=a[n-1]*gh*dp_nm_dth;
      
      	//Bpm-=a[n-1]*dgh_dphi*p_nm; 
      
      	Bp_gg-=aa*ggg;
      	Bp_hh-=aa*hhh; 
      
      
      	p_nm1=p_nm; //corresponds to P for n-1 and m
      	p_nm2=0.;  //corresponds to P for n-2 and m
      	dp_nm1_dth=dp_nm_dth;
     	 dp_nm2_dth=0.;
        
     
        for ( n = m + 1 ; n < old_nm_gauss +1 ; n++ ){
        	G4int ii= ((n+2) * (n-1) )/2 + m;
		double rc= recurrence_coef[ii];
         	p_nm = cos_theta *p_nm1 - rc*p_nm2;
         	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                                 - rc*dp_nm2_dth;
	 
		
		gg=gg_nm1[ii];
		hh=hh_nm1[ii];
	 
	 
	 	double ggg = gg*p_nm;
         	double hhh = hh*p_nm;
		da=da_dr[n-1];
		aa=a[n-1];
      
         	Br_gg-=da*ggg;
        	Br_hh-=da*hhh;
                
		w=aa*dp_nm_dth;
         	Bt_gg-=w*gg;
         	Bt_hh-=w*hh;
      
         	Bp_gg-=aa*ggg;
         	Bp_hh-=aa*hhh; 

	 	p_nm2=p_nm1;
	 	p_nm1=p_nm;
	 	dp_nm2_dth=dp_nm1_dth;
	 	dp_nm1_dth=dp_nm_dth;
	}
	Br+=Br_gg*cos_mphi + Br_hh*sin_mphi;
	Bt+=Bt_gg*cos_mphi + Bt_hh*sin_mphi;
	Bp+=double(m)*(Bp_hh*cos_mphi - Bp_gg*sin_mphi);

  } 
  
  if (sin_theta >= 1.e-8) Bp/=sin_theta; 
  else if (cos_theta<0.) Bp=-Bp;
   
  G4double Be= Br*sin_theta+ Bt*cos_theta;
  G4double Bx= Be*cos_phi- Bp*sin_phi;
  G4double By= Be*sin_phi+ Bp*cos_phi;
  G4double Bz= Br*cos_theta - Bt*sin_theta;
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
    
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz)*nanotesla;
  return Bfield;
}

////////////////////////////////////////////////////////////////////////////////
/////////////// External field models

G4ThreeVector EarthMagneticField::GetTSY89( G4ThreeVector pos) const
{  bdip_.b0dip= DipoleB0/nanotesla;
   SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

   G4ThreeVector pos_geo=pos;
   if (ConsiderDipoleShift) pos_geo-= DipoleShift;  
   //GEOtoGSM
   G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
   
   float xgsm=float(pos_gsm.x()/re);
   float ygsm=float(pos_gsm.y()/re);
   float zgsm=float(pos_gsm.z()/re);
  
 
   
   float parmod[10]; 
   float Bxgsm,Bygsm,Bzgsm;
   
   float ps=float(theConvertor->GetTiltAngle());
   if(TiltAngle != -1000) ps = TiltAngle; //PVD
   
   G4int ipar=iopt;
   t89c_(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                      &Bxgsm, &Bygsm, &Bzgsm );
   
   G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
   Bgsm *= nanotesla;
   
   //GSMtoGEO
   return theConvertor->TransformPSMinPLA(Bgsm);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector EarthMagneticField::GetTSYBOB89( G4ThreeVector pos) const
{  bdip_.b0dip= DipoleB0/nanotesla;
   SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

   G4ThreeVector pos_geo=pos;
   if (ConsiderDipoleShift) pos_geo-= DipoleShift;  
   //GEOtoGSM
   G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
   
   float xgsm=float(pos_gsm.x()/re);
   float ygsm=float(pos_gsm.y()/re);
   float zgsm=float(pos_gsm.z()/re);
  
 
   
   float parmod[10]; 
   parmod[1]=Dst;
   float Bxgsm,Bygsm,Bzgsm;
   float ps=float(theConvertor->GetTiltAngle());
   if(TiltAngle != -1000) ps = TiltAngle; //PVD
   G4int ipar=iopt;
   t89c_boberg_(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                      &Bxgsm, &Bygsm, &Bzgsm );
   
   G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
   Bgsm *= nanotesla;
   
   //GSMtoGEO
   return theConvertor->TransformPSMinPLA(Bgsm);
}
////////////////////////////////////////////////////////////////////////////
G4ThreeVector EarthMagneticField::GetTSY96(G4ThreeVector pos) const
{ bdip_.b0dip= DipoleB0/nanotesla;
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

  G4ThreeVector pos_geo=pos;
  if (ConsiderDipoleShift) pos_geo-= DipoleShift;  
  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
  

  float parmod[10]; 
  parmod[0]=Pdyn;
  parmod[1]=Dst/nanotesla; 
  parmod[2]=Imfy/nanotesla;
  parmod[3]=Imfz/nanotesla;
   
   
  float Bxgsm,Bygsm,Bzgsm;
  
  float ps=float(theConvertor->GetTiltAngle());
  if(TiltAngle != -1000) ps = TiltAngle; //PVD
  G4int ipar=iopt;
  
  mcos_t96_01_(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                        &Bxgsm, &Bygsm, &Bzgsm );

  G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
   
  //GSMtoGEO
  return theConvertor->TransformPSMinPLA(Bgsm);

}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector
   EarthMagneticField::GetTSY2001( G4ThreeVector pos) const
{bdip_.b0dip= DipoleB0/nanotesla; 
 SpaceCoordinatePlanet* theConvertor =
                             SpaceCoordinatePlanet::GetInstance(); 

   
 G4ThreeVector pos_geo=pos;
 
 // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
 //GEOtoGSM
 G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
 
 float xgsm=float(pos_gsm.x()/re);
 float ygsm=float(pos_gsm.y()/re);
 float zgsm=float(pos_gsm.z()/re);
  
// call tsy2001 fortran subroutine   
   
 float parmod[10]; 
 parmod[0]=Pdyn;
 parmod[1]=Dst/nanotesla; 
 parmod[2]=Imfy/nanotesla;
 parmod[3]=Imfz/nanotesla;
 parmod[4]=G1;
 parmod[5]=G2;
   
 float Bxgsm,Bygsm,Bzgsm;
 float ps=float(theConvertor->GetTiltAngle());
 if(TiltAngle != -1000) ps = TiltAngle; //PVD
   
 G4int ipar=iopt;
 t01_01_(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                        &Bxgsm, &Bygsm, &Bzgsm );
   
 G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
 Bgsm *= nanotesla;
   
 //GSMtoGEO
 return theConvertor->TransformPSMinPLA(Bgsm);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector
   EarthMagneticField::GetTSY2004( G4ThreeVector pos) const
{bdip_.b0dip= DipoleB0/nanotesla; 
 SpaceCoordinatePlanet* theConvertor =
                             SpaceCoordinatePlanet::GetInstance(); 

   
 G4ThreeVector pos_geo=pos;
 
 // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
 //GEOtoGSM
 G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
 
 float xgsm=float(pos_gsm.x()/re);
 float ygsm=float(pos_gsm.y()/re);
 float zgsm=float(pos_gsm.z()/re);
  
// call tsy2004 fortran subroutine   
   
 float parmod[10]; 
 parmod[0]=Pdyn;
 parmod[1]=Dst/nanotesla; 
 parmod[2]=Imfy/nanotesla;
 parmod[3]=Imfz/nanotesla;
 parmod[4]=W1;
 parmod[5]=W2;
 parmod[6]=W3;
 parmod[7]=W4;
 parmod[8]=W5;
 parmod[9]=W6;
   
 float Bxgsm,Bygsm,Bzgsm;
 float ps=float(theConvertor->GetTiltAngle());
 if(TiltAngle != -1000) ps = TiltAngle; //PVD
   
 G4int ipar=iopt;
 
 t04_s_(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                        &Bxgsm, &Bygsm, &Bzgsm );
   
 G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
 Bgsm *= nanotesla;
   
 //GSMtoGEO
 return theConvertor->TransformPSMinPLA(Bgsm);
}
////////////////////////////////////////////////////////////////////////////////
//-------------------------------Outside magnetosphere
bool EarthMagneticField::OutsideMagnetosphere( G4ThreeVector pos) const
{ return  (this->*POutsideMagnetosphere)(pos);
}
////////////////////////////////////////////////////////////////////////////////
//
bool EarthMagneticField::IGRFOutsideMagnetosphere( G4ThreeVector pos)  const
{ G4double r=pos.mag()/re;
  bool res=false;
  if (r > 25.) res =true;
  return res;
}
////////////////////////////////////////////////////////////////////////////////
//
bool EarthMagneticField::TSY89OutsideMagnetosphere( G4ThreeVector pos) const
{ //model from Kobel PhD "Zu die magnetosphärischen Effekten der kosmischen
  //  Strahlung"
  G4double a=-0.0545;
  G4double b_iopt[]={11.7,11.1,10.8,10.4,10.4,10.2,10.2};
  G4double f_iopt[]={20.,15.,6.67,10.,5.,6.,6.};
  
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 
  G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos)/re;
  G4double rho2 = pos_gsm.y()*pos_gsm.y() + pos_gsm.z()*pos_gsm.z();
  G4bool res = (rho2>900 || pos_gsm.x() <-60);
  if (!res){ 
  	double  k_par = theConvertor->GetTiltAngle()/f_iopt[iopt-1];
  	G4double cos_k=std::cos(k_par);
  	G4double sin_k=std::sin(k_par);
  	G4RotationMatrix mat = G4RotationMatrix(G4ThreeVector(cos_k,0,sin_k),
  					  G4ThreeVector(0,1,0), 
					  G4ThreeVector(-sin_k,0,cos_k));
  	G4ThreeVector pos_rot = mat*pos_gsm;
	G4double rho_rot= pos_rot.y()*pos_rot.y() + pos_rot.z()*pos_rot.z();
	res = (pos_rot.x() > a*rho_rot +b_iopt[iopt-1]); 
  }
 						  
  return res;
}
///////////////////////////////////////////////////////////////////////////////
//
bool EarthMagneticField::TSY96OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko96 code
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

  G4ThreeVector pos_geo=pos;
  
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  
  // magnetopause model
  G4double am0=70.;
  G4double s0=1.08;
  G4double x00=5.48;
  G4double dsig=.005;
  G4double xappa=std::pow(Pdyn/2.,0.14);
  G4double x0=x00/xappa;
  G4double am=am0/xappa;
  G4double rho2=ygsm*ygsm+zgsm*zgsm;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
///////////////////////////////////////////////////////////////////////////////
//
bool EarthMagneticField::TSY2001OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko2001 code
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

  G4ThreeVector pos_geo=pos;
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  G4double r= pos_gsm.mag()/re;
  
  G4double xss=xgsm;
  G4double zss=zgsm;
  
  G4double rh0=8.94335;
  G4double dd;
  G4double rh2= -5.2;
  
  float sps = std::sin(theConvertor->GetTiltAngle());
  if(TiltAngle != -1000) sps = std::sin(TiltAngle);//PVD

  //change of coordinate due to the tail warping
  do{ G4double xsold=xss;
      G4double zsold=zss;
      G4double rh=rh0+rh2*(zss/r)*(zss/r);
      G4double sinpsas=sps/std::pow(1.+std::pow(r/rh,3.),1./3.);
      G4double cospsas=std::sqrt(1. - sinpsas*sinpsas);
      zss=xgsm*sinpsas+zgsm*cospsas;
      xss=xgsm*cospsas-zgsm*sinpsas;
      dd=std::abs(xss-xsold)+std::abs(zss-zsold);
  } while (dd > 1.e-6);
  
  // magnetopause model
  G4double a0_a=34.586;
  G4double s0=1.196;
  G4double a0_x0=3.4397;
  G4double dsig=.003;
  G4double xappa=std::pow(Pdyn/2.,0.15821);
  G4double x0=a0_x0/xappa;
  G4double am=a0_a/xappa;
  G4double rho2=ygsm*ygsm+zss*zss;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}


///////////////////////////////////////////////////////////////////////////////
//
bool EarthMagneticField::TSY2004OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko2004 code	
	
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 

  G4ThreeVector pos_geo=pos;
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformPLAinPSM(pos_geo);  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  G4double r= pos_gsm.mag()/re;
  
  G4double xss=xgsm;
  G4double zss=zgsm;
  
  G4double rh0=7.5;
  G4double dd;
  G4double rh2= -5.2;
  
  float sps = std::sin(theConvertor->GetTiltAngle());
  if(TiltAngle != -1000) sps = std::sin(TiltAngle);//PVD
  
  //change of coordinate due to the tail warping
  do{ G4double xsold=xss;
      G4double zsold=zss;
      G4double rh=rh0+rh2*(zss/r)*(zss/r);
      G4double sinpsas=sps/std::pow(1.+std::pow(r/rh,3.),1./3.);
      G4double cospsas=std::sqrt(1. - sinpsas*sinpsas);
      zss=xgsm*sinpsas+zgsm*cospsas;
      xss=xgsm*cospsas-zgsm*sinpsas;
      dd=std::abs(xss-xsold)+std::abs(zss-zsold);
  } while (dd > 1.e-6);
  
  // magnetopause model
  G4double a0_a=34.586;
  G4double s0=1.196;
  G4double a0_x0=3.4397;
  G4double dsig=.005;
  G4double xatpa=std::pow(Pdyn/2.,0.152759);
  G4double x0=a0_x0/xatpa;
  G4double am=a0_a/xatpa;
  G4double rho2=ygsm*ygsm+zss*zss;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));

  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
  
}

////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ComputeTSY2001Parameter(G4double t)
{ G4int index =-1;
  if (n_tsy_data>0){
  	for (G4int i=0;i<n_tsy_data-1;i++){
   		if ( (t >= times_of_data[i]) && (t <= times_of_data[i+1])){
       			index=i;
        		i=n_tsy_data-1;
		}
	}
  	if (index > -1){
    		G4double dx=(t-times_of_data[index])
                	/(times_of_data[index+1]-times_of_data[index]);
     		Pdyn=Pdyn_data[index] + dx * (Pdyn_data[index+1]-Pdyn_data[index]);
     		Dst=Dst_data[index] + dx * (Dst_data[index+1]-Dst_data[index]);
     		Dst*=nanotesla;
		Imfy=Imf_gsm_y_data[index] + dx * 
                      (Imf_gsm_y_data[index+1]-Imf_gsm_y_data[index]);
		Imfy*=nanotesla;      
     		Imfz=Imf_gsm_z_data[index] + dx * 
                      (Imf_gsm_z_data[index+1]-Imf_gsm_z_data[index]);
     		Imfz*=nanotesla;
		G1=g1_tsy_param_data[index] + dx * 
                       (g1_tsy_param_data[index+1]-g1_tsy_param_data[index]);
     		G2=g2_tsy_param_data[index] + dx * 
                      (g2_tsy_param_data[index+1]-g2_tsy_param_data[index]);     
  	}
  }  		
}

////////////////////////////////////////////////////////////////////////////////
//
void EarthMagneticField::ComputeTSY2004Parameter(G4double t)
{ //to be done
	
	/*G4int index =-1;
  if (n_tsy_data>0){
  	for (G4int i=0;i<n_tsy_data-1;i++){
   		if ( (t >= times_of_data[i]) && (t <= times_of_data[i+1])){
       			index=i;
        		i=n_tsy_data-1;
		}
	}
  	if (index > -1){
    		G4double dx=(t-times_of_data[index])
                	/(times_of_data[index+1]-times_of_data[index]);
     		Pdyn=Pdyn_data[index] + dx * (Pdyn_data[index+1]-Pdyn_data[index]);
     		Dst=Dst_data[index] + dx * (Dst_data[index+1]-Dst_data[index]);
     		Dst*=nanotesla;
		Imfy=Imf_gsm_y_data[index] + dx * 
                      (Imf_gsm_y_data[index+1]-Imf_gsm_y_data[index]);
		Imfy*=nanotesla;      
     		Imfz=Imf_gsm_z_data[index] + dx * 
                      (Imf_gsm_z_data[index+1]-Imf_gsm_z_data[index]);
     		Imfz*=nanotesla;
		G1=g1_tsy_param_data[index] + dx * 
                       (g1_tsy_param_data[index+1]-g1_tsy_param_data[index]);
     		G2=g2_tsy_param_data[index] + dx * 
                      (g2_tsy_param_data[index+1]-g2_tsy_param_data[index]);     
  	}
  }  */		
}

