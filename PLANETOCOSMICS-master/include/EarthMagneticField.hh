#ifndef EarthMAGNETICFIELD_HH
#define EarthMAGNETICFIELD_HH 

#include "globals.hh"
#include"G4ios.hh"

#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"
#include "PlanetMagneticField.hh"

class EarthMagneticFieldMessenger;


class EarthMagneticField : public PlanetMagneticField
{
public:
	 //constructor destructor	       
         EarthMagneticField() ;
	 ~EarthMagneticField() ;
	 
	 
	 //Methods that should be provided for specific planet
	 virtual void SetInternalField(G4String aString);
	 virtual void SetExternalField(G4String aString);
	 virtual void SetMagnetopauseModel(G4String aString);
	 virtual void ComputeBfieldParametersAtReferenceTime();    
         bool OutsideMagnetosphere( G4ThreeVector pos) const;
	 G4ThreeVector GetInternalField (G4ThreeVector) const;
	 G4ThreeVector GetExternalField (G4ThreeVector) const;
	 
	 //Set methods specific to Earth 
	 void SetIopt(G4int val); 
	 //inline void Setnm_igrf(G4int aVal) {nm_igrf=aVal;}
		 
	 
	 //Method for geomagnetic and sw parameters used in the Tsyganenko
	 //models 
	 inline void SetTiltAngle(G4double aVal){TiltAngle = aVal;}
	 inline void SetPdyn(G4double aVal){Pdyn = aVal;}
	 inline void SetDst(G4double aVal){Dst = aVal;}
	 inline void SetImfy(G4double aVal){Imfy = aVal;}
	 inline void SetImfz(G4double aVal){Imfz = aVal;}
	 inline void SetG1(G4double aVal){G1 = aVal;}
	 inline void SetG2(G4double aVal){G2 = aVal;}
	 
	 inline void SetW1(G4double aVal){W1 = aVal;}	 
	 inline void SetW2(G4double aVal){W2 = aVal;}	 
	 inline void SetW3(G4double aVal){W3 = aVal;}	 
	 inline void SetW4(G4double aVal){W4 = aVal;}	 
	 inline void SetW5(G4double aVal){W5 = aVal;}	 
	 inline void SetW6(G4double aVal){W6 = aVal;}	 	 
	 
	 void PrintStormParameterTSY2001();
	 void PrintStormParameterTSY2004();
	 void ReadTSY2001Parameter(G4String nameFile);
	 void ReadTSY2004Parameter(G4String nameFile);
	 void ReadIgrfTable(G4String name_file);
	 void ComputeIgrfCoefAccordingToTime();
	 void SetTiltedDipoleParameterFromIGRF();
	 void SetEccentricDipoleParameterFromIGRF();
	 
	 inline void SetConsiderDipoleShift(G4bool abool)
	                               {ConsiderDipoleShift = abool;}
				       
private:

        // messenger
         EarthMagneticFieldMessenger* theFieldMessenger;

        //Magnetic field model parameters  
	G4double year_igrf;
	G4int iopt;
	std::vector< double >  times_of_data, Pdyn_data, Dst_data, 
			       Imf_gsm_y_data,Imf_gsm_z_data,
	         	       g1_tsy_param_data, g2_tsy_param_data,
			       w1_tsy_param_data, w2_tsy_param_data, w3_tsy_param_data, w4_tsy_param_data, w5_tsy_param_data, w6_tsy_param_data;
        G4int n_tsy_data,nm_igrf;
	G4double TiltAngle, Pdyn,Dst,Imfy,Imfz,G1,G2, W1, W2, W3, W4, W5, W6;
	
	//Igrf coefficient table
	std::vector<float>  igrf_year;
	std::vector< std::vector<float > > h_nm;
	std::vector< std::vector<float > > g_nm;
	std::vector<float>  dh_nm;
	std::vector<float>  dg_nm;
		
	//Igrf coefficient according to reference time 
	std::vector<double> hh_nm1;
	std::vector<double> gg_nm1;
	
	float* hh;
	float* gg;
	float* rec;
	
	
	

//private methods       
private:     
        //IGRF model
	G4ThreeVector  GetIGRFFortran(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRFFortran1(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRF(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRF1(G4ThreeVector pos) const;
	G4ThreeVector  GetIGRFDOUBLE(G4ThreeVector pos) const;  
	
	//Tsyganenko89 model 
	G4ThreeVector  GetTSY89(G4ThreeVector pos) const;
	G4ThreeVector  GetTSYBOB89(G4ThreeVector pos) const;
	
	//Tsyganenko96 model 
	G4ThreeVector  GetTSY96(G4ThreeVector pos) const;
	
	//Tsyganenko2001 model
	G4ThreeVector  GetTSY2001(G4ThreeVector pos) const;
	
	//Tsyganenko2004 model
	G4ThreeVector  GetTSY2004(G4ThreeVector pos) const;
	
	//pointers on the selected Internal and External field model
	G4ThreeVector (EarthMagneticField::* PInternalField)(G4ThreeVector)const;
	G4ThreeVector (EarthMagneticField::* PExternalField)(G4ThreeVector)const;
	
	
	//pointer on magnetopause model
	bool (EarthMagneticField::*POutsideMagnetosphere)
	                                             (G4ThreeVector pos) const;
        
	//IGRF and geomagnetic magnetopause model r>25. re is outside
	//magnetosphere

	bool IGRFOutsideMagnetosphere(G4ThreeVector pos) const; 
	
	
	//tsy89 magnetopause model
	bool TSY89OutsideMagnetosphere(G4ThreeVector pos) const;
	
	//tsy89 magnetopause model
	bool TSY96OutsideMagnetosphere(G4ThreeVector pos) const;
        
	//tsy2001 magnetopause model
	bool TSY2001OutsideMagnetosphere(G4ThreeVector pos) const;
	
	//tsy2004 magnetopause model
	bool TSY2004OutsideMagnetosphere(G4ThreeVector pos) const;
	
	void ComputeTSY2001Parameter(G4double t);	   
	
	void ComputeTSY2004Parameter(G4double t);	   
}; 

#endif
