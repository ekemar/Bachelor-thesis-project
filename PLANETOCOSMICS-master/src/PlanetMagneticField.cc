
#include "PlanetMagneticField.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4MagIntegratorStepper.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "SpaceCoordinatePlanet.hh"
#include "PlanetUnits.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4UImanager.hh"
#include"G4RunManager.hh"
#include"BlineTool.hh"
#include"PlanetEquationOfMotion.hh"
#include"PlanetManager.hh"
#include"DateAndTime.hh"
#include "PLANETOCOSPrimaryGeneratorAction.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4SolidStore.hh"
#include "G4Box.hh"
#include "MagneticShieldingTool.hh"

////////////////////////////////////////////////////////////////////////////////
//
PlanetMagneticField::PlanetMagneticField(G4String aName)
{ 
  ListOfInternalFieldModels.clear();
  ListOfInternalFieldModels.push_back("NOFIELD");
  ListOfInternalFieldModels.push_back("FLATGRID");
  ListOfInternalFieldModels.push_back("DIPOLE");
  ListOfInternalFieldModels.push_back("GAUSS");
  
  ListOfExternalFieldModels.clear();
  ListOfExternalFieldModels.push_back("NOFIELD");
  
  ListOfMagnetopauseModels.clear();
  ListOfMagnetopauseModels.push_back("NOMAGNETOPAUSE");
  ListOfMagnetopauseModels.push_back("SPHERE");
  
  
   
  //Planet name
  PlanetName = aName;
  
  //Time
  TimeOfB=0.; 
  StartDate = DateAndTime(2000,1,1);
  ReferenceDate = StartDate;
  
  //EquationOfMotion  
  fEquationOfMotion = new PlanetEquationOfMotion(this);
   
  //Define parametersfor the integration method
  fChordFinder=0;
  fStepper=0;
  rplanet=PlanetManager::GetInstance()->GetRplanet();
  RadiusMagnetosphere = 25.*rplanet;
  DefaultDeltaChord=.1*rplanet;
  DeltaChord=DefaultDeltaChord;
  DefaultDeltaIntersection= .01*rplanet;
  DefaultEpsilon= 1.e-6;
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()
	                        ->SetLargestAcceptableStep(10.*rplanet);
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->SetDetectorField(this);
  ResetIntegrationParameters(); 
  
  
  IsFlatGeometry=true;

}
////////////////////////////////////////////////////////////////////////
//
PlanetMagneticField::~PlanetMagneticField()
{ if (fEquationOfMotion) delete fEquationOfMotion;
  if (fChordFinder) delete fChordFinder;
  if (fStepper) delete fStepper;
}
///////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::GetFieldValue( const G4double yTrack[7], G4double B[3]) const 
{ 
  G4ThreeVector position = G4ThreeVector(yTrack[0],yTrack[1],yTrack[2]);
  G4ThreeVector Bfield = GetFieldValue(position); 
  B[0]=Bfield.x();
  B[1]=Bfield.y();
  B[2]=Bfield.z();
 
  return ;
  
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector PlanetMagneticField::GetFieldValue(const G4ThreeVector position) const 
{ G4ThreeVector Bfield; 
  Bfield = G4ThreeVector(0.,0.,0.);
  //G4cout<<"BfieldConst"<<BfieldConst<<std::endl;
  if (InternalModelName == "CONST") Bfield = BfieldConst;
  else {
  	G4ThreeVector PLAPosition = position;
  
  	G4bool Rotate = (IsFlatGeometry && 
                   RotatePosAndField && 
		   InternalModelName != "FLATGRID" );
  
  	if (Rotate){ 
  		PLAPosition = PLAPosition.rotateY(ZenithDirectionInPLA.theta())
                                  	 	 .rotateZ(ZenithDirectionInPLA.phi());
  				 
  		PLAPosition += WorldCentPLAPosition;
  	}
  
  	if (Internal) Bfield = GetInternalField(PLAPosition);	
  	if (External) Bfield = Bfield + GetExternalField(PLAPosition); 		
 
 
 	if (Rotate)
    		Bfield = Bfield.rotateZ(-ZenithDirectionInPLA.phi()).
  				     rotateY(-ZenithDirectionInPLA.theta());
  }				     	
  
  return Bfield;
  
 

}			       
/////////////////////////////////////////////////////////////////////////////////
//
/*G4ThreeVector PlanetMagneticField::FindPositionOnMagnetopause(G4double theta, 
                                                              G4double xgsm, 
							      G4double precision) const
{ G4ThreeVector pos1 = G4ThreeVector(xgsm,0.,0.);
 
  SpaceCoordinatePlanet* theCoordinateConvertor
                            = SpaceCoordinatePlanet::GetInstance();
  theCoordinateConvertor->SetSystemInAndOut("GSM","PLA");
 
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

////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector PlanetMagneticField::FindStandOffPosition(G4double precision ) const
{ G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.)*Re;
  G4ThreeVector pos2 = G4ThreeVector(40.,0.,0.)*Re;
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
  	pos3= (pos1+pos2)/2.;
	if (OutsideMagnetosphere(pos3)) pos2=pos3;
	else pos1=pos3;
  }
  return pos1;	
}
*/
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetTimeOfB(G4double val)
{ 
  G4cout<<val<<std::endl;
  TimeOfB=val;
  G4double JulDate = StartDate.JulianDate() + ( val / (24. * 3600.) );
  ReferenceDate = DateAndTime(JulDate);
  ComputeBfieldParametersAtReferenceTime();
  SpaceCoordinatePlanet* theCoordinateConvertor 
  				= SpaceCoordinatePlanet::GetInstance();
  theCoordinateConvertor->SetReferenceDate(ReferenceDate);
  if (HasAGlobalField){
  	SetTiltedDipoleParameterFromGAUSSCoefficients();
	G4ThreeVector DipoleAxis = G4ThreeVector(1.,0.,0.);
	DipoleAxis.setTheta(DipoleTheta);
	DipoleAxis.setPhi(DipolePhi);
	theCoordinateConvertor->SetDipoleAxisInPLA(DipoleAxis);
        
  }
 
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetStartDate(DateAndTime aDate)
{
  // Set the start date
  StartDate = aDate;
  double val =TimeOfB;
  SetTimeOfB(val);
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetStartDate(G4int year, G4int month,G4int day,
	                               G4int hour, G4int minute,G4int second )
{
  // Set the start date
  StartDate = DateAndTime(year,month,day,hour,minute,second);
  double val =TimeOfB;
  SetTimeOfB(val);
}		
////////////////////////////////////////////////////////////////////////////////
//
bool PlanetMagneticField::SetMotherInternalField(  G4String aString)
{ Internal=true;
  
  InternalModelName = aString;
  
  if (aString =="DIPOLE") {
  	PMotherInternalField=&PlanetMagneticField::GetPLADipole;
  } 
  else if (aString == "CONST"){
  	PMotherInternalField=&PlanetMagneticField::GetBConst;
        ExternalModelName = "NOFIELD";
  }	
  else if (aString == "FLATGRID"){
  	PMotherInternalField=&PlanetMagneticField::GetBInterpolFlat;
	ExternalModelName = "NOFIELD";
  }	
  
  else Internal=false;
  return Internal;     
}


////////////////////////////////////////////////////////////////////////////////
//
bool  PlanetMagneticField::SetMotherMagnetopauseModel( G4String aString) 
{ 
  if  (aString == "SPHERE") {
       PMotherOutsideMagnetosphere
                 = &PlanetMagneticField::OutsideSphericalMagnetosphere;
       return true;
       Magnetopause =true;	
  }
  return false;	  	    

}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetDipoleB0(G4double aVal)
{ DipoleB0=aVal;
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetDipoleAxis(G4double Theta,G4double Phi)
{ DipoleTheta=Theta;
  DipolePhi=Phi;
  G4ThreeVector dipole_axis_in_PLA;
  dipole_axis_in_PLA.setRThetaPhi(1.,Theta, Phi);
  SpaceCoordinatePlanet::GetInstance()->SetDipoleAxisInPLA(dipole_axis_in_PLA); 

}


////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetEpsilon(G4double aVal)
{ //precision for G4Integration methods
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMinimumEpsilonStep(aVal);
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMaximumEpsilonStep(aVal*1.00001);
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetDeltaChord(G4double aVal)
{ if (fChordFinder) fChordFinder->SetDeltaChord(aVal);
  DeltaChord=aVal;
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetDeltaIntersection(G4double aVal)
{ ///////Delta intersection for G4Integration method
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                                                 ->SetDeltaIntersection(aVal);

}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ResetIntegrationParameters()
{ SetStepper("CashKarpRKF45");
  SetDeltaChord(DefaultDeltaChord);
  SetDeltaIntersection(DefaultDeltaIntersection);
  SetEpsilon(DefaultEpsilon);
}

////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ReadGaussCoefficients(G4String name_file)
{std::fstream File_Input(name_file,  std::ios::in);
 
 //clear vectors
  hh_nm.clear();
  gg_nm.clear();
  recurrence_coef.clear();
 
  nmax =0; 
  nm_gauss=0;
  G4int n,m;
  
  // initialisation h_nm,g_nm;  
  
  
  
   
  while (!File_Input.eof()){
    	n=-999;
    	m=-999;
    	File_Input>>n>>m;
	//G4cout<<n<<'\t'<<m<<std::endl;;
    	if (n != -999 && m != -999) {
    		
     		if (n>nmax){
      			nmax=n;
       			hh_nm.insert(hh_nm.end(),nmax+1,0.);
	 		gg_nm.insert(gg_nm.end(),nmax+1,0.); 
      		}
     		G4int index = ((n+2)  * (n-1)) /2 + m;
     		G4double value;
		File_Input>>value;
		gg_nm[index] = value;
		//G4cout<<value<<'\t';
		if (m>0) {
			File_Input>>value;
			hh_nm[index] = value;
			//G4cout<<value<<'\t';
		}
		else 	hh_nm[index]=0;
		//G4cout<<std::endl; 
		
     				
    	}
  }
  nm_gauss=nmax;
 // G4cout<<"nmax "<<nmax<<std::endl;
 
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
       			 double ((n-1+m)*(n-1-m)) / double ((2*n-1)*(2*n-3)) );
 File_Input.close();
 
 //The coefficient are  Schmidt quasi normalised
 // For using the recuurence relation in Gauss normalisation
 // we have to convert the coefficient in Gauss normalisation
 // the g_nm and h_nm coefficient become
  G4double s=1.;
  for (n=1;n<nmax+1;n++){
  	G4int mn = ((n + 1) * n) / 2  -1;
       	s*=double(2*n-1)/double(n);
       	hh_nm[mn] *= s;
       	gg_nm[mn] *= s;
       	G4double p=s;
       	for (int m=1;m<n+1;m++){
        	G4double aa=1.;
	 	if (m == 1) aa=2.;
	 	p *= std::sqrt(aa * double(n - m + 1) / double (n + m) );
	 	int nmn =  mn + m ;
	 	hh_nm[nmn] *= p;
         	gg_nm[nmn] *= p;
	}
  }  
 
 
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector PlanetMagneticField::GetGAUSS(G4ThreeVector pos) const 
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


  float r=pos.mag()/rplanet;
  float theta=pos.theta();
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  float phi=pos.phi();
   
  float pp = 1. / r;
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
 
  for ( n=1;n<nm_gauss+1;n++){
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
    
 for (n = 2 ; n < nm_gauss +1 ; n++ ){
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
  for (int m = 1 ; m < nm_gauss +1 ; m++ ){
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
      
      	float ggg = gg_nm[ii]*p_nm;
      	float hhh = hh_nm[ii]*p_nm;
      
      	Br_gg-=da_dr[n-1]*ggg;
      	Br_hh-=da_dr[n-1]*hhh;
      
      	Bt_gg-=a[n-1]*gg_nm[ii]*dp_nm_dth;
      	Bt_hh-=a[n-1]*hh_nm[ii]*dp_nm_dth;
      
      	//Bt-=a[n-1]*gh*dp_nm_dth;
      
      	//Bpm-=a[n-1]*dgh_dphi*p_nm; 
      
      	Bp_gg-=a[n-1]*ggg;
      	Bp_hh-=a[n-1]*hhh; 
      
      
      	p_nm1=p_nm; //corresponds to P for n-1 and m
      	p_nm2=0.;  //corresponds to P for n-2 and m
      	dp_nm1_dth=dp_nm_dth;
     	 dp_nm2_dth=0.;
        
     
        for ( n = m + 1 ; n < nm_gauss +1 ; n++ ){
        	G4int ii= ((n+2) * (n-1) )/2 + m; 
         	p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
         	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                                 - recurrence_coef[ii]*dp_nm2_dth;
	 
		 /*gh = gg_nm[ii]*cos_mphi + hh_nm[ii]*sin_mphi;
         	dgh_dphi = hh_nm[ii]*cos_mphi - gg_nm[ii]*sin_mphi;
        
	 	Br-=da_dr[n-1]*gh*p_nm;
         	Bt-=a[n-1]*gh*dp_nm_dth;
         	Bpm-=a[n-1]*dgh_dphi*p_nm;*/ 
	 
	 
	 	float ggg = gg_nm[ii]*p_nm;
         	float hhh = hh_nm[ii]*p_nm;
      
         	Br_gg-=da_dr[n-1]*ggg;
        	 Br_hh-=da_dr[n-1]*hhh;
      
         	Bt_gg-=a[n-1]*gg_nm[ii]*dp_nm_dth;
         	Bt_hh-=a[n-1]*hh_nm[ii]*dp_nm_dth;
      
         	Bp_gg-=a[n-1]*ggg;
         	Bp_hh-=a[n-1]*hhh; 

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
void PlanetMagneticField::PrintBfield(G4ThreeVector pla_pos) const
{ SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 
		   
  G4ThreeVector Bfield =GetFieldValue(pla_pos);
  Bfield /= nanotesla; 
  G4cout<<"Bfield at primary position" <<std::endl;
  if (IsFlatGeometry){
  	G4cout<<"Bfield : "<<Bfield<<" "<<Bfield.mag()<<std::endl;
  	return;
  }	
  
  
 
 
  
  G4cout<<"Bfield PLA: "<<Bfield<<" "<<Bfield.mag()<<std::endl;
  G4double altitude, longitude, latitude;
  theConvertor->ComputePLAGCoordinatesFromPLAPosition(pla_pos,
 						     altitude,
						     longitude,latitude);
  //Northward, eastward, dowmward component
  //---------------------------------------
 
  G4ThreeVector vertical = G4ThreeVector(1.,0.,0.);
  vertical.setRThetaPhi(1.,90.*degree-latitude,longitude);
  G4ThreeVector east_dir =G4ThreeVector(-std::sin(longitude),std::cos(longitude),0. );
  G4ThreeVector north_dir = -east_dir.cross(vertical);
  G4double Bdownward = -Bfield.dot(vertical);
  G4double Bnorthward = Bfield.dot(north_dir);
  G4double Beastward = Bfield.dot(east_dir);
 
  G4cout<<"Bfield north, east, down: "<<Bnorthward
                       <<" "<<Beastward<<" "<<Bdownward<<std::endl;
		       
  if (theConvertor->GetPlanetIsMagnetic()) {		       
  	G4ThreeVector Bpsm = theConvertor->Transform(Bfield, "PLA", "PSM");
  	G4ThreeVector Bmag = theConvertor->Transform(Bfield, "PLA", "PMAG");
  	G4ThreeVector Bsmag = theConvertor->Transform(Bfield, "PLA", "PSMAG");		       
  	G4cout<<"Bfield PSM: "<<Bpsm<<" "<<Bpsm.mag()<<std::endl;
  	G4cout<<"Bfield SMAG: "<<Bsmag<<" "<<Bsmag.mag()<<std::endl;
  	G4cout<<"Bfield MAG: "<<Bmag<<" "<<Bmag.mag()<<std::endl; 
  }	
  
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ComputeBfieldAtDifferentPosition(G4String cosys_pos, 
				G4double alt0, G4double dAlt, G4int nAlt,
				G4double lat0, G4double dlat, G4int nlat,
				G4double long0, G4double dlong, G4int nlong,
				G4String file_output) const
{ clock_t clock1,clock2;
  if (IsFlatGeometry){
  	G4cout<<"This command is not available in the case of flat geometry"<<std::endl;
  	return;
  }
  clock1=clock();
  SpaceCoordinatePlanet* theConvertor =
                   SpaceCoordinatePlanet::GetInstance(); 
  std::ofstream OutputFile;
  OutputFile.open(file_output, (std::ios::binary | std::ios::out));
  OutputFile<<"System of coordinate for position: "<<cosys_pos<<std::endl;
  OutputFile<<"altitude [km]"<<'\t'<<"latitude "<<'\t'
				  <<"longitude "<<'\t'
				  <<"Bx"<<'\t'<<'\t'<<"By"<<'\t'<<'\t'
				  <<"Bz"<<'\t'<<'\t'<<"B"<<'\t'<<'\t'
				  <<"Br"<<'\t'<<'\t'<<"Bth"<<'\t'<<'\t'<<"Bphi"<<'\t'<<'\t'
				  <<"Bvert"<<'\t'<<'\t'<<"Bnorth"<<std::endl;
 /* OutputFile.precision(2);
  OutputFile.setf(std::ios::fixed); */
  OutputFile<<std::setiosflags(std::ios::scientific);
  OutputFile<<std::setprecision(4);
  G4double Bmax = 0.;	
  G4double alt_max=0.;
  G4double lat_max = 0.;
  G4double lon_max =0.;
  for (double  altitude=alt0;altitude<alt0 + double(nAlt)*dAlt; altitude+=dAlt){
  	for (double  latitude=lat0; latitude< lat0 + double(nlat)*dlat; 
							latitude+=dlat){
		for (double  longitude=long0; longitude< long0 + double(nlong)*dlong;
							 longitude+=dlong){
			G4ThreeVector pla_pos;
			if (cosys_pos != "PLAG"){
				G4ThreeVector pos=G4ThreeVector(0.,0.,1.);
				pos.setRThetaPhi(rplanet+altitude,
						  90*degree-latitude,
						  longitude);
				pla_pos = theConvertor->Transform(pos, cosys_pos,"PLA"); 
	           	}
			else {
				theConvertor->ComputePLAPositionFromPLAG
	            			(altitude, latitude, longitude, pla_pos);
			}
			G4ThreeVector Bfield =GetFieldValue(pla_pos);
  			Bfield /= nanotesla;
			G4double plag_alt,plag_long,plag_lat;
			//G4cout<<"Test1"<<std::endl;
			theConvertor->ComputePLAGCoordinatesFromPLAPosition
								(pla_pos,
 						     		 plag_alt,
						     		 plag_long,
								 plag_lat);
			//G4cout<<"Test2"<<std::endl;
			G4ThreeVector vertical = G4ThreeVector(1.,0.,0.);
  			//G4cout<<(90.*degree-plag_lat)/degree<<std::endl;
			vertical.setRThetaPhi(1.,
					      90.*degree-plag_lat,
					      plag_long);
  			G4ThreeVector east_dir =G4ThreeVector(-std::sin(longitude),std::cos(longitude),0. );
  			G4ThreeVector north_dir = -east_dir.cross(vertical);
  			G4double Bdown = -Bfield.dot(vertical);
  			G4double Bnorth = Bfield.dot(north_dir);
  			G4double Beast = Bfield.dot(east_dir);
			
			vertical = pla_pos/pla_pos.mag();
			G4double Br= Bfield.dot(vertical);
			G4ThreeVector theta_dir = east_dir.cross(vertical);
			G4double Bth = Bfield.dot(theta_dir);
			
			
			OutputFile<<altitude/km<<'\t'
				  <<latitude/degree<<'\t'
				  <<longitude/degree<<'\t'
				  <<Bfield.x()<<'\t'<<Bfield.y()<<'\t'
				  <<Bfield.z()<<'\t'<<Bfield.mag()<<'\t'
				  <<Br<<'\t'<<Bth<<'\t'<<Beast<<'\t'
				  <<Bdown<<'\t'<<Bnorth<<std::endl;
			if (Bfield.mag() > Bmax) {
				Bmax= Bfield.mag();
				lat_max =latitude/degree;
				lon_max =longitude/degree;
				alt_max =altitude/km; 
			}	              	
		}
	}
  
  }
  OutputFile<<alt_max<<'\t'<<lat_max<<'\t'<<lon_max<<'\t'<<Bmax<<std::endl;
  OutputFile.close();
   
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the computation: "<<tclock<<" s"<<std::endl;
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ComputeBfieldAtDifferentPosition(
				G4double alt0, G4double dAlt, G4int nAlt,
				G4double x0,G4double  dx, G4int nX,
				G4double y0,G4double  dy, G4int nY,
				G4String file_output) const
{ clock_t clock1,clock2;
  
  clock1=clock();
  if (!IsFlatGeometry){
  	G4cout<<"This command is not available in the case of spherical geometry"<<std::endl;
  	return;
  }
  std::ofstream OutputFile;
  OutputFile.open(file_output, (std::ios::binary | std::ios::out));
  OutputFile<<"altitude [km]"<<'\t'<<"X "<<'\t'
				  <<"Y "<<'\t'
				  <<"Bx"<<'\t'<<"By"<<'\t'
				  <<"Bz"<<'\t'
				  <<"B"<<std::endl;
 /* OutputFile.precision(2);
  OutputFile.setf(std::ios::fixed);*/
  OutputFile<<std::setiosflags(std::ios::scientific);
  OutputFile<<std::setprecision(4);
 
  for (double  altitude=alt0;altitude<alt0 + double(nAlt)*dAlt; altitude+=dAlt){
  	for (double  x=x0; x< x0 + double(nX)*dx; x+=dx){
		for (double  y=y0; y< y0 + double(nY)*dy; y+=dy){
		
			
			double z= altitude-altitude_world_center;
			G4ThreeVector Bfield =
				GetFieldValue(G4ThreeVector(x,y,z));
			Bfield /= nanotesla;
			
			
			OutputFile<<altitude/km<<'\t'<<'\t'
				  <<x/km<<'\t'<<'\t'
				  <<y/km<<'\t'<<'\t'
				  <<Bfield.x()<<'\t'<<Bfield.y()<<'\t'
				  <<Bfield.z()<<'\t'<<Bfield.mag()<<'\t'
				  <<std::endl;	              	
		}
	}
  
  }
  
  OutputFile.close();
   
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the computation: "<<tclock<<" s"<<std::endl;
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::TraceBlineFromDifferentPositions(G4String coorsys, 
				G4double altitude,
				G4double lat0, G4double dLat, G4int nLat,
				G4double long0, G4double dLong, G4int nLong)
{ if (IsFlatGeometry){
  	G4cout<<"This command is not available in the case of flat geometry"<<std::endl;
  	return;
  }
  
  
  PLANETOCOSPrimaryGeneratorAction* thePrimaryAction =
 	dynamic_cast<PLANETOCOSPrimaryGeneratorAction*>(
		const_cast< G4VUserPrimaryGeneratorAction*>(
 		  G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));   
 
 //G4cout<<"test1"<<std::endl;
// G4UImanager * UI = G4UImanager::GetUIpointer();
 //BlineTool* theBlineTool = BlineTool::GetInstance();
 for (double  latitude=lat0; latitude< lat0 + double(nLat)*dLat; 
							latitude+=dLat){
		for (double  longitude=long0; longitude< long0 +double(nLong)*dLong;
							 longitude+=dLong){
			
			/*G4cout<<latitude/degree<<std::endl;
			G4cout<<longitude/degree<<std::endl;*/
			thePrimaryAction->SetPosition(coorsys,
						      altitude,
						      longitude,
						      latitude
						      );	
			MagneticShieldingTool::GetInstance()->Bline();
			
			//UI->ApplyCommand("/PLANETOCOS/DRAW/TraceBline"); 
			
		}
 } 		
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::TraceBlineFromDifferentPositions(
				G4double altitude,
				G4double x0, G4double dx, G4int nx,
				G4double y0, G4double dy, G4int ny)
{ if (!IsFlatGeometry){
  	G4cout<<"This command is not available in the case of spherical geometry"<<std::endl;
  	return;
  }
  PLANETOCOSPrimaryGeneratorAction* thePrimaryAction =
 	dynamic_cast<PLANETOCOSPrimaryGeneratorAction*>(
		const_cast< G4VUserPrimaryGeneratorAction*>(
 		  G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));   
 
 //G4cout<<"test1"<<std::endl;
 //BlineTool* theBlineTool = BlineTool::GetInstance();
 //G4UImanager * UI = G4UImanager::GetUIpointer();
 for (double  x=x0; x< x0 + double(nx)*dx; x+=dx){
		for (double  y=y0; y< y0 + double(ny)*dy; y+=dy){
		
			//G4cout<<x/km<<'\t'<<y/km<<'\t'<<altitude/km<<std::endl;
			thePrimaryAction->SetPosition(G4ThreeVector(x,y,altitude));
						  
			MagneticShieldingTool::GetInstance()->Bline();
			//UI->ApplyCommand("/PLANETOCOS/DRAW/TraceBline");  		      
			
			
		}
 } 		
}

//////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ComputeOrientation()
{ 
 SpaceCoordinatePlanet* theConvertor = 
                                SpaceCoordinatePlanet::GetInstance();
  //G4cout<<"refsystem_world_center"<<refsystem_world_center<<std::endl;
  if (refsystem_world_center == "PLAG"){
   	theConvertor->ComputePLADirectionAndPositionFromPLAG
	                 (0., latitude_world_center,longitude_world_center,
			  0., 0.,WorldCentPLAPosition ,ZenithDirectionInPLA);
    	ZenithDirectionInPLA = -ZenithDirectionInPLA;		  
    	WorldCentPLAPosition+=ZenithDirectionInPLA*altitude_world_center;	  		     
    	return;
  }
    
  G4double cos_lat = std::cos(latitude_world_center);
  G4double  sin_lat = std::sin(latitude_world_center);
  G4double cos_lon =std::cos(longitude_world_center);
  G4double sin_lon =std::sin(longitude_world_center);
  
  WorldCentPLAPosition=(rplanet + altitude_world_center)* 
                                        G4ThreeVector(cos_lat*cos_lon,
                                                      cos_lat*sin_lon,
						      sin_lat);
   
  ZenithDirectionInPLA = WorldCentPLAPosition / WorldCentPLAPosition.mag();
 
						
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetWorldCenterPosition( G4double latitude,
	 			 		       	G4double longitude, 
				 		       	G4String coord_sys)
{ latitude_world_center = latitude;
  longitude_world_center = longitude;
  refsystem_world_center = coord_sys;
  //G4cout<<coord_sys<<"coord_sys"<<std::endl;
  ComputeOrientation();
 // ComputeInterpolMatrixForFlatGeometry(100,100,100,"test_100.out");
 //  ReadInterpolMatrixForFlatGeometry("test_100.out");
  
}

////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ComputeOrientationAndBfieldConst
	                         (G4double altitude, 
	                          G4double latitude, G4double longitude,
	                          G4String reference_system,
			          G4String internal_model,
				  G4String external_model)				 
{ 
  SpaceCoordinatePlanet* theConvertor = 
  	                        SpaceCoordinatePlanet::GetInstance();
 
  G4String old_internal = InternalModelName;
  G4String old_external = ExternalModelName;
  SetInternalField(internal_model);
  SetExternalField(external_model);
  if (reference_system== "PLAG"){ 
   	theConvertor->ComputePLAPositionFromPLAG
	                 (altitude, latitude,longitude,WorldCentPLAPosition);
    	G4ThreeVector position;
    	theConvertor->ComputePLADirectionAndPositionFromPLAG(altitude, 
							      latitude,
							      longitude,
			  				      0., 0., 
							      position,
							     ZenithDirectionInPLA);
    	ZenithDirectionInPLA= -ZenithDirectionInPLA;
	G4bool aBool = 	RotatePosAndField;
	RotatePosAndField = false;		  		     
    	BfieldConst = GetFieldValue(position);
	RotatePosAndField = aBool;
    	// now orient BfieldConst to the world Z axis 
    	if  (RotatePosAndField) {
   		BfieldConst= BfieldConst.rotateZ(-ZenithDirectionInPLA.phi())
                           		.rotateY(-ZenithDirectionInPLA.theta());  
        }      
   	InternalModelName = "CONST";
	return;
    	
  }
    
  G4double cos_lat = std::cos(latitude);
  G4double  sin_lat = std::sin(latitude);
  G4double cos_lon =std::cos(longitude);
  G4double sin_lon =std::sin(longitude);
  WorldCentPLAPosition=(rplanet + altitude) * G4ThreeVector(cos_lat*cos_lon,
                                                       cos_lat*sin_lon,
						       sin_lat);
  G4ThreeVector position = rplanet * WorldCentPLAPosition /(rplanet + altitude);
 						      
  ZenithDirectionInPLA = WorldCentPLAPosition / WorldCentPLAPosition.mag();
  
  G4bool aBool = 	RotatePosAndField;
  RotatePosAndField = false;
  BfieldConst = GetFieldValue(WorldCentPLAPosition);
  RotatePosAndField = aBool; 
    
  // now orient BfieldConst to the world Z axis 
  if  (RotatePosAndField) {    
  	BfieldConst= BfieldConst.rotateZ(-ZenithDirectionInPLA.phi())
                          	.rotateY(-ZenithDirectionInPLA.theta()); 
  }
  InternalModelName = "CONST";					
}


////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector PlanetMagneticField::GetBInterpolFlat(G4ThreeVector pos) const
{G4int i,j,k;
 G4double x,y,z;
 G4int nx = Bx_flat.size()-1;
 G4int ny= Bx_flat[0].size()-1;
 G4int nz= Bx_flat[0][0].size()-1;
 x=pos.x();
 y=pos.y();
 z=pos.z();
 i=int( (x-xmin)/dx);
 if (i <0) i=0;
 else if (i >= nx) i = nx-1;
 float fx = (x- (xmin+i*dx))/dx;

 j=int( (y-ymin)/dy);
 if (j<0) j=0;
 else if (j>= ny) j=ny-1;
 float fy = (y- (ymin+j*dy))/dy;
 
 k=int( (z-zmin)/dz);
 if (k<0) k=0;
 else if (k>=nz) k=nz-1; 
 float fz = (z- (zmin+k*dz)) /dz;
 
/* G4double fac = (1.-fx);
 fac *= (1.-fy);
 fac *=(1.-fz);
 G4double Bx= */

 G4double Bx= Bx_flat[i][j][k]*(1.-fx)*(1.-fy)*(1.-fz) +
 	      Bx_flat[i][j][k+1]*(1.-fx)*(1.-fy)*fz+
	      Bx_flat[i][j+1][k]*(1.-fx)*fy*(1.-fz) +
 	      Bx_flat[i][j+1][k+1]*(1.-fx)*fy*fz+
	      Bx_flat[i+1][j][k]*fx*(1.-fy)*(1.-fz) +
 	      Bx_flat[i+1][j][k+1]*fx*(1.-fy)*fz+
	      Bx_flat[i+1][j+1][k]*fx*fy*(1.-fz) +
 	      Bx_flat[i+1][j+1][k+1]*fx*fy*fz;
 
 G4double By= By_flat[i][j][k]*(1.-fx)*(1.-fy)*(1.-fz) +
 	      By_flat[i][j][k+1]*(1.-fx)*(1.-fy)*fz+
	      By_flat[i][j+1][k]*(1.-fx)*fy*(1.-fz) +
 	      By_flat[i][j+1][k+1]*(1.-fx)*fy*fz+
	      By_flat[i+1][j][k]*fx*(1.-fy)*(1.-fz) +
 	      By_flat[i+1][j][k+1]*fx*(1.-fy)*fz+
	      By_flat[i+1][j+1][k]*fx*fy*(1.-fz) +
 	      By_flat[i+1][j+1][k+1]*fx*fy*fz;	      
	      
 G4double Bz= Bz_flat[i][j][k]*(1.-fx)*(1.-fy)*(1.-fz) +
 	      Bz_flat[i][j][k+1]*(1.-fx)*(1.-fy)*fz+
	      Bz_flat[i][j+1][k]*(1.-fx)*fy*(1.-fz) +
 	      Bz_flat[i][j+1][k+1]*(1.-fx)*fy*fz+
	      Bz_flat[i+1][j][k]*fx*(1.-fy)*(1.-fz) +
 	      Bz_flat[i+1][j][k+1]*fx*(1.-fy)*fz+
	      Bz_flat[i+1][j+1][k]*fx*fy*(1.-fz) +
 	      Bz_flat[i+1][j+1][k+1]*fx*fy*fz;	
	      	      
	      
 
 return G4ThreeVector(Bx,By,Bz)*nanotesla;
 
 
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector PlanetMagneticField::GetPLADipole(G4ThreeVector pos) const
{ G4double xpla,ypla,zpla; 
  xpla=pos.x()/rplanet;
  ypla=pos.y()/rplanet;
  zpla=pos.z()/rplanet;
   
  G4ThreeVector SchiftedPosition = 
  G4ThreeVector(xpla,ypla,zpla)-(DipoleShift/rplanet); 
   
  G4double sth=sin(DipoleTheta);
  G4double cth=cos(DipoleTheta);
  G4double sphi=sin(DipolePhi);
  G4double cphi=cos(DipolePhi);
   
  //coordinates in DipoleSystem
  G4double xd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*cth
                -sth*SchiftedPosition.z();
  G4double yd=SchiftedPosition.y()*cphi-SchiftedPosition.x()*sphi;
  G4double zd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*sth
                +cth*SchiftedPosition.z();
   
   
  //component in dipole system		
  G4double Bxd,Byd,Bzd;
  G4double r2=xd*xd +yd *yd +zd*zd; 
  G4double r5=std::pow(r2,2.5);
  G4double factor=DipoleB0/r5;
  Bxd=-3*factor*xd*zd;
  Byd=-3*factor*yd*zd;
  Bzd=factor*(r2-3*zd*zd);
   
   
  //back to Geo coordinate
  G4double Bx,By,Bz;
  Bx=(cth*cphi*Bxd - sphi*Byd  + sth*cphi*Bzd);
  By=(cth*sphi*Bxd + cphi*Byd  + sth*sphi*Bzd);
  Bz=(-sth*Bxd  + cth*Bzd);
   
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}


////////////////////////////////////////////////////////////////////////////////
//
bool PlanetMagneticField::OutsideSphericalMagnetosphere( G4ThreeVector pos)  const
{ G4double r=pos.mag();
  if (r > RadiusMagnetosphere) return true;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetStepper(G4String aString)
{ if (fStepper) delete fStepper;
   
  if (aString == "ExplicitEuler"){
    	fStepper = new G4ExplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ExplicitEuler is called"<<G4endl;
  }
  else if  (aString == "ImplicitEuler"){
    	fStepper = new G4ImplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ImplicitEuler is called"<<G4endl;
  }
  else if  (aString == "SimpleRunge"){
    	fStepper = new G4SimpleRunge( fEquationOfMotion );
     	G4cout<<"G4SimpleRunge is called"<<G4endl;
  } 
  else if  (aString == "ClassicalRK4"){
    	fStepper = new G4ClassicalRK4( fEquationOfMotion );
     	G4cout<<"G4ClassicalRK4 is called"<<G4endl;
  }  
  else if  (aString == "RKG3_Stepper"){
    	fStepper = new G4RKG3_Stepper( fEquationOfMotion );
     	G4cout<<"G4RKG3_Stepper is called"<<G4endl;
  }
  else if  (aString == "CashKarpRKF45"){
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }   
  else {
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"The selected stepper is not available"<<G4endl;
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }
  
  if (fChordFinder) delete fChordFinder;
  G4double min_step= .000001*re;
  fChordFinder = new G4ChordFinder(this,min_step,fStepper);
  fChordFinder->SetDeltaChord(DeltaChord); 
  
  G4FieldManager* fieldMgr= 
  G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetChordFinder(fChordFinder); 

}
////////////////////////////////////////////////////////////////////
void PlanetMagneticField::SwitchOn()
{G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->SetDetectorField(this);
}
////////////////////////////////////////////////////////////////////
void PlanetMagneticField::SwitchOff()
{G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->SetDetectorField(0);
}

///////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ComputeInterpolMatrixForFlatGeometry(int nx, 
							     int ny, 
							     int nz,
							     G4String file_name)

{ PLANETOCOSGeometryConstruction * theDetector = 
    dynamic_cast<PLANETOCOSGeometryConstruction*>
 	(const_cast<G4VUserDetectorConstruction*> 
	       (G4RunManager::GetRunManager()->GetUserDetectorConstruction())) ;std::ofstream OutputFile;
  if (theDetector->GetGeometryType() =="SPHERICAL"){
  	G4cout<<"This method can only be used for Flat geometry"<<std::endl;
	return;
  }
  
  
  OutputFile.open(file_name, std::ios::binary);
  OutputFile<<nx<<'\t'<<ny<<'\t'<<nz<<std::endl;
  
  
  Bx_flat.clear();
  By_flat.clear();
  Bz_flat.clear();
 
  G4SolidStore* theSolidStore = G4SolidStore::GetInstance();
  unsigned int index= theSolidStore->size();
  for (unsigned int i=0;i<theSolidStore->size();i++){
  	if ((*theSolidStore)[i]->GetName() == "Magnetosphere"){
		index =i;
		i=theSolidStore->size();
	
	} 
  	
  }	
  G4Box* theSolidMagnetosphere = dynamic_cast<G4Box*>((*theSolidStore)[index]);
  
  
  zmin=-theSolidMagnetosphere->GetZHalfLength();
  zmax=theSolidMagnetosphere->GetZHalfLength(); 
  xmin=-theSolidMagnetosphere->GetXHalfLength();
  xmax=theSolidMagnetosphere->GetXHalfLength();
  dx=2.*theSolidMagnetosphere->GetXHalfLength()/double(nx);
  dy=2.*theSolidMagnetosphere->GetYHalfLength()/double(ny);
  dz=(zmax-zmin)/double(nz);
  ymin=-theSolidMagnetosphere->GetYHalfLength();
  ymax=theSolidMagnetosphere->GetYHalfLength();     
  G4double x=xmin-dx;
  G4double y;
  G4double z;
  for (int i=0;i<nx+1;i++){
 	G4cout<<"i "<<i<<std::endl;
 	x+=dx;
	Bx_flat.push_back(std::vector< std::vector<float> >());
	By_flat.push_back(std::vector< std::vector<float> >());
	Bz_flat.push_back(std::vector< std::vector<float> >());
	y=ymin-dy;
	for (int j=0;j<ny+1;j++){
	        y+=dy;
		Bx_flat[i].push_back(std::vector<float>());
	        By_flat[i].push_back(std::vector<float>());
	        Bz_flat[i].push_back(std::vector<float>());
		z=zmin-dz;
		for (int k=0;k<nz+1;k++){
		        
			z+=dz;
			G4ThreeVector Bfield = GetFieldValue(
						G4ThreeVector(x,y,z))/nanotesla;
			Bx_flat[i][j].push_back(Bfield.x());
			By_flat[i][j].push_back(Bfield.y());
			Bz_flat[i][j].push_back(Bfield.z());
			
			OutputFile<<float(Bfield.x())<<'\t'
				  <<float(Bfield.y())<<'\t'
				  <<float(Bfield.z())<<std::endl;;
			
		}
	}
	
 }
}
///////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::ReadInterpolMatrixForFlatGeometry(G4String file_name)

{ PLANETOCOSGeometryConstruction * theDetector = 
    dynamic_cast<PLANETOCOSGeometryConstruction*>
 	(const_cast<G4VUserDetectorConstruction*> 
	       (G4RunManager::GetRunManager()->GetUserDetectorConstruction())) ;std::ofstream OutputFile;
  if (theDetector->GetGeometryType() =="SPHERICAL"){
  	G4cout<<"This method can only be used for Flat geometry"<<std::endl;
	return;
  }
  
  std::ifstream InputFile;
  InputFile.open(file_name, std::ios::in);
  G4int nx,ny,nz;
  InputFile>>nx>>ny>>nz;
  G4cout<<nx<<ny<<nz<<std::endl;
  
  
  Bx_flat.clear();
  By_flat.clear();
  Bz_flat.clear();
  
  G4SolidStore* theSolidStore = G4SolidStore::GetInstance();
  unsigned int index= theSolidStore->size();
  for (unsigned int i=0;i<theSolidStore->size();i++){
  	if ((*theSolidStore)[i]->GetName() == "Magnetosphere"){
		index =i;
		i=theSolidStore->size();
	
	} 
  	
  }	
  G4Box* theSolidMagnetosphere = dynamic_cast<G4Box*>((*theSolidStore)[index]);
  
  zmin=-theSolidMagnetosphere->GetZHalfLength();
  zmax=theSolidMagnetosphere->GetZHalfLength(); 
  xmin=-theSolidMagnetosphere->GetXHalfLength();
  xmax=theSolidMagnetosphere->GetXHalfLength();
  dx=2.*theSolidMagnetosphere->GetXHalfLength()/double(nx);
  dy=2.*theSolidMagnetosphere->GetYHalfLength()/double(ny);
  dz=(zmax-zmin)/double(nz);
  ymin=-theSolidMagnetosphere->GetYHalfLength();
  ymax=theSolidMagnetosphere->GetYHalfLength();  
 
  
  
  G4double x=xmin-dx;
  G4double y;
  G4double z;
  for (int i=0;i<nx+1;i++){
 	G4cout<<"i "<<i<<std::endl;
 	x+=dx;
	Bx_flat.push_back(std::vector< std::vector<float> >());
	By_flat.push_back(std::vector< std::vector<float> >());
	Bz_flat.push_back(std::vector< std::vector<float> >());
	y=ymin-dy;
	for (int j=0;j<ny+1;j++){
	        y+=dy;
		Bx_flat[i].push_back(std::vector<float>());
	        By_flat[i].push_back(std::vector<float>());
	        Bz_flat[i].push_back(std::vector<float>());
		z=zmin-dz;
		for (int k=0;k<nz+1;k++){
		        float a,b,c;
			InputFile>>a>>b>>c;
			Bx_flat[i][j].push_back(a);
			By_flat[i][j].push_back(b);
			Bz_flat[i][j].push_back(c);
			
			
			
		}
	}
	
 }
}
///////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetTiltedDipoleParameterFromGAUSSCoefficients()
{ if (hh_nm.size()>0 ) {
	G4double g10,h11,g11;
 	g10= gg_nm[0];
       	g11= gg_nm[1];
       	h11= hh_nm[1];
	G4ThreeVector dipole_axis_in_PLA = -G4ThreeVector(g11,h11,g10)*nanotesla;
	DipoleB0 = dipole_axis_in_PLA.mag();
	dipole_axis_in_PLA = dipole_axis_in_PLA/DipoleB0 ;
	SpaceCoordinatePlanet::GetInstance()->SetDipoleAxisInPLA(dipole_axis_in_PLA);
	DipoleTheta = dipole_axis_in_PLA.theta() ;
	DipolePhi = dipole_axis_in_PLA.phi() ;
  }
  else {
    	G4cout<<"The spherical harmonic coefficients have not been  defined!"<<std::endl;
	      
  }


}
///////////////////////////////////////////////////////////////////////
//
void PlanetMagneticField::SetEccentredDipoleParameterFromGAUSSCoefficients()
{if (hh_nm.size()>0 ) {
	
	G4double g10,h11,g11,g20,g21,g22,h21,h22;
 	g10= gg_nm[0];
       	g11= gg_nm[1];
       	h11= hh_nm[1];
       	g20= gg_nm[2]/1.5;
       	g21= gg_nm[3]/(std::sqrt(4./3.) * 1.5);
       	h21= hh_nm[3]/(std::sqrt(4./3.) * 1.5);
       	g22= gg_nm[4]/(std::sqrt(1./3.) * 1.5);
       	h22= hh_nm[4]/(std::sqrt(1./3.) * 1.5);
	
	G4ThreeVector dipole_axis_in_PLA = -G4ThreeVector(g11,h11,g10)*nanotesla;
	DipoleB0 = dipole_axis_in_PLA.mag();
	dipole_axis_in_PLA = dipole_axis_in_PLA/DipoleB0 ;
	SpaceCoordinatePlanet::GetInstance()->SetDipoleAxisInPLA(dipole_axis_in_PLA);
	DipoleTheta = dipole_axis_in_PLA.theta() ;
	DipolePhi = dipole_axis_in_PLA.phi() ;
        G4double rd2=DipoleB0*DipoleB0/(nT*nT);
	
	
	//computation of DipoleSchift in GEO according to 
	//Langel, Geomagnetism 1, p 381-390 with a correction (-sign) in equation 228b p386
	//Fraser-Smith, Rev. Geophys., 25, 1,1-16,1987 there the formula seems correct
	//Spenvis bacgroubd url: 
	//   http://www.spenvis.oma.be/spenvis/help/background/magfield/cd.html    
     	G4double sqr3=std::sqrt(3.); 
     	G4double L0,L1,L2,T;
     	L0=2.*g10*g20+sqr3*(g11*g21+h11*h21);
     	L1=-g11*g20+sqr3*(g10*g21+g11*g22+h11*h22);
     	L2=-h11*g20+sqr3*(g10*h21+g11*h22-h11*g22);
     	T=(L0*g10+L1*g11+L2*h11)/(4.*rd2);
     	DipoleShift = G4ThreeVector(re*(L1-g11*T)/(3.*rd2),
                            re*(L2-h11*T)/(3.*rd2),
			    re*(L0-g10*T)/(3.*rd2));
     
      }
  else {
    	G4cout<<"The spherical harmonic coefficients have not been  defined!"<<std::endl;
	      
  }
}
