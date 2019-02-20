#ifndef PlanetMAGNETICFIELD_HH
#define PlanetMAGNETICFIELD_HH 

#include "globals.hh"
#include"G4ios.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"

class G4ChordFinder;
class G4MagIntegratorStepper;
class PlanetEquationOfMotion;

class PlanetMagneticField : public G4MagneticField
{
public:
	 //constructor destructor	       
         PlanetMagneticField(G4String aName);
	 ~PlanetMagneticField();
	     
         // Gives the magnetic field B at a given position defined by yTrack
	 void GetFieldValue(const G4double yTrack[],G4double B[]) const;
	 G4ThreeVector GetFieldValue(const G4ThreeVector position) const; 
	 
         //Methods that should be provided for specific planet
	 virtual void SetInternalField(G4String aString) =0;
	 virtual void SetExternalField(G4String aString) =0;
	 virtual G4ThreeVector GetInternalField(G4ThreeVector) const =0 ;
	 virtual G4ThreeVector GetExternalField(G4ThreeVector) const =0 ;
	 virtual void SetMagnetopauseModel(G4String aString)=0;
	 virtual void ComputeBfieldParametersAtReferenceTime()=0;
	 virtual bool OutsideMagnetosphere(G4ThreeVector pos) const =0;
	 
	 //Sphericalharmonics
	 void ReadGaussCoefficients(G4String name_file);
	 
	 //Set methods 
	 ///////////////////
	 
	 bool SetMotherInternalField(G4String aString);
	 bool SetMotherMagnetopauseModel(G4String aString);
	 void SetTimeOfB(G4double val);	
	 void SetStartDate(G4int year, G4int month,G4int day,
	                   G4int hour, G4int minute,G4int second );
	 void SetStartDate(DateAndTime aDate);
	 		   
	 inline void Setnm_gauss(G4int n){nm_gauss=n;}
	 //Dipole definition
	 void SetDipoleB0(G4double aVal);
	 void SetDipoleAxis(G4double Theta,G4double Phi);
	 inline void SetDipoleShift(G4ThreeVector aVec){DipoleShift=aVec;}
	 
	 //type of geometry		   
	 inline void SetIsFlatGeometry(G4bool aBool){IsFlatGeometry=aBool;}
	 
	 inline void SetRadiusMagnetosphere(G4double aRadius) {RadiusMagnetosphere =aRadius;}
	
	 //Integration paremeters
	 void ResetIntegrationParameters();
	 void SetEpsilon(G4double aVal);
	 void SetDeltaChord(G4double aVal);
	 void SetDeltaIntersection(G4double aVal);
	 void SetStepper(G4String );
	 
	 inline void SetBfieldConst(G4ThreeVector aVec){BfieldConst=aVec;
	 						InternalModelName= "CONST";}
	 inline void SetRotatePosAndField(G4bool abool){RotatePosAndField=abool;}
	 inline void SetAltitudeWorldCenter(G4double aVal)
	                            {altitude_world_center = aVal;} 			       
	 void SetWorldCenterPosition(G4double latitude,
	 			     G4double longitude, 
				     G4String coord_sys);
	 inline void SetHalfLength(G4double aVal ){half_length = aVal;};			     
	  
	 
	 
	 //Get methods
	 /////////////////////
	 
	 inline G4ThreeVector GetDipoleShift(){return DipoleShift;}
	 inline G4double GetDipoleB0(){return DipoleB0;}
	 inline G4double GetDipolePhi(){return DipolePhi;}
	 inline G4double GetDipoleTheta(){return DipoleTheta;}								
	 inline void SetConsiderDipoleShift(G4bool abool){ConsiderDipoleShift = abool;}
	 inline bool GetHasAGlobalField(){return HasAGlobalField;}
	 inline std::vector< G4String > GetListOfInternalFieldModels() {return ListOfInternalFieldModels;}
	 inline std::vector< G4String > GetListOfExternalFieldModels() {return ListOfExternalFieldModels;}			       
	 inline std::vector< G4String > GetListOfMagnetopauseModels() {return ListOfMagnetopauseModels;}
	 inline PlanetEquationOfMotion* GetEquationOfMotion()
	                                {return fEquationOfMotion;}
	
	 inline DateAndTime GetReferenceDate(){return ReferenceDate;}
	 inline G4String GetIntFieldModelName(){return InternalModelName;}
	 inline G4String GetExtFieldModelName(){return ExternalModelName;}
	
	
	 void PrintBfield(G4ThreeVector geo_pos) const; 
	 // compute orientation and BfieldConst at a given altitude and longitude
	 void ComputeOrientation();
	 void ComputeOrientationAndBfieldConst(G4double altitude,
				  G4double latitude, G4double longitude,
	                          G4String reference_system,
			          G4String internal_model,
				  G4String external_model="NOFIELD");			       
				       
	//Compute Bfield at different positions			       
	 void ComputeBfieldAtDifferentPosition(G4String cosys_pos, 
				G4double alt0, G4double dAlt, G4int nAlt,
				G4double lat0, G4double dLat, G4int nLat,
				G4double long0, G4double dLong, G4int nLong,
				G4String file_output) const;
	 
	 void ComputeBfieldAtDifferentPosition(
	   			G4double alt0, G4double dAlt, G4int nAlt,
				G4double x0,G4double  dx0, G4int nX,
				G4double y0,G4double  dy0, G4int nY,
				G4String file_output) const;
	 
	 
	 //Blines
	 void TraceBlineFromDifferentPositions(G4String cosys_pos, 
				G4double altitude,
				G4double lat0, G4double dLat, G4int nLat,
				G4double long0, G4double dLong, G4int nLong);			 			       
	 void TraceBlineFromDifferentPositions( 
				G4double altitude,
				G4double x0, G4double dX, G4int nX,
				G4double y0,G4double  dy0, G4int nY);
	 
	 //SwitchOn, Switch off field
	 void SwitchOn();
	 void SwitchOff();
	 
	 //Interpolation 
	 void ComputeInterpolMatrixForFlatGeometry(int nx, int ny, int nz,
	 						G4String file_name);
	 void ReadInterpolMatrixForFlatGeometry(G4String file_name);
	 void ComputeInterpolMatrixForSphericalGeometry(int nx, int ny, int nz,
	 						G4String file_name);
	 
	 void ReadInterpolMatrixForSphericalGeometry(G4String file_name);									       				
	 void SetTiltedDipoleParameterFromGAUSSCoefficients();
	 void SetEccentredDipoleParameterFromGAUSSCoefficients();
protected:
        
        //Magnetic field model parameters  
	
	std::vector< G4String > ListOfInternalFieldModels;
	std::vector< G4String > ListOfExternalFieldModels;
	std::vector< G4String > ListOfMagnetopauseModels;
        
	G4bool External,Internal,ConsiderDipoleShift; 
	G4bool Magnetopause;
	G4double DipoleB0,DipoleTheta,DipolePhi;
	G4ThreeVector DipoleShift;
	G4double RadiusMagnetosphere;
	G4String InternalModelName, ExternalModelName; 
	
	// planet radius
	G4double rplanet;
	
	//Name
	G4String PlanetName;
	
	
	//Spherical harmonic field
	///////////////////////////
	
	//coefficient for the recurrence equation of legendre assocated function
	std::vector<float> recurrence_coef;
	G4int nm_gauss,nmax;
	
	//spherical coefficient 
	std::vector<float> hh_nm;
	std::vector<float> gg_nm;
	
	//attribute for the integration method 
	///////////////////////////////////////
	
	G4ChordFinder* fChordFinder;
	G4MagIntegratorStepper* fStepper;
	G4double DefaultDeltaChord;
	G4double DeltaChord;
	G4double DefaultDeltaIntersection;
	G4double DefaultEpsilon;
	
	PlanetEquationOfMotion* fEquationOfMotion;
	
	//define the orientation of geomagnetic field and world center for 
	// flat geometry
	//////////////////////////////////////////////////////////////////
	G4ThreeVector ZenithDirectionInPLA;
	G4ThreeVector WorldCentPLAPosition;
	G4double altitude_world_center;
	G4double latitude_world_center;
	G4double longitude_world_center;
	G4String refsystem_world_center;
	G4ThreeVector BfieldConst;
	
	G4bool RotatePosAndField;
	
	//Geometry
	G4bool IsFlatGeometry;
	
	//Global Field
	G4bool HasAGlobalField;
	
	
	//Matrices for  interpolation
	//////////////////////////////
	
	std::vector< std::vector< std::vector<float> > > Bx_flat;
	std::vector< std::vector< std::vector<float> > > By_flat;
	std::vector< std::vector< std::vector<float> > > Bz_flat;
	G4double half_length,xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz;
	
	//Time
	DateAndTime StartDate,ReferenceDate;
	G4double TimeOfB;	

//protected methods
protected:
	//magnetic dipole in PLA coordinate
	G4ThreeVector  GetPLADipole(G4ThreeVector pos) const;
	
	//Constant field
	inline G4ThreeVector GetBConst(G4ThreeVector )const{return BfieldConst;} ;
	
	//Interpolation for flat geometry
	inline G4ThreeVector GetBInterpolFlat(G4ThreeVector pos) const;
	G4ThreeVector  GetGAUSS(G4ThreeVector pos) const;
	
	
	//Pointer of Internal field declared in Mother class
	G4ThreeVector (PlanetMagneticField::* PMotherInternalField)(G4ThreeVector) const;
	
	
	//pointer on magnetopause model declared in Mother class
	bool (PlanetMagneticField::*PMotherOutsideMagnetosphere)
	                                             (G4ThreeVector pos) const;
        
	
	//Spahericalk magnetopuase model with radius defined by 
	bool OutsideSphericalMagnetosphere(G4ThreeVector pos) const; 
	
	
	
		   
} ;

#endif
