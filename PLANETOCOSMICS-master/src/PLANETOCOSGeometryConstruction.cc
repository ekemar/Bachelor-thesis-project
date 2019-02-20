#include "PLANETOCOSGeometryConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4PropagatorInField.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4CSGSolid.hh"
#include "PLANETOCOSGeometryMessenger.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "PLANETOCOSEdepAnalyser.hh"
#include "PLANETOCOSFluxDetectionAnalyser.hh"
#include "PLANETOCOSSteppingAction.hh"

#include "PlanetMagneticField.hh"
#include "PlanetAtmosphereModel.hh"
#include "PlanetSoil.hh"
#include "PlanetManager.hh"
#include "PlanetUnits.hh"
#include "PLANETOCOSSD.hh"
#include "PLANETOCOSSoilSD.hh"
#include "SpaceCoordinatePlanet.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "myfunctions.hh"
#include "G4Orb.hh"
#include "G4FieldManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSGeometryConstruction::PLANETOCOSGeometryConstruction()
{ 
  //Some elements 
  //-----------
  elH = new G4Element( "Hydrogen", "H", 1., 1.00797*g/mole);
  elHe = new G4Element( "Helium", "He", 2., 4.00260*g/mole);
  elC = new G4Element( "Carbon", "C", 6. , 12.01*g/mole );
  elN = new G4Element( "Nitrogen", "N", 7. , 14.00674*g/mole );
  elO = new G4Element( "Oxygen", "O", 8. , 16.*g/mole );
  elAr = new G4Element( "Argon", "Ar", 18., 39.948*g/mole );
  
  new G4Element( "Sodium", "Na", 11., 22.9898*g/mole );
  new G4Element( "Magnesium", "Mg", 12., 24.3051*g/mole );
  new G4Element( "Aluminium", "Al", 13., 26.98154*g/mole );
  new G4Element( "Silicon", "Si", 14., 27.9769*g/mole );
  new G4Element( "Potassium", "K", 19, 39.0983*g/mole);
  new G4Element( "Calcium", "Ca", 20., 40.078*g/mole );
  new G4Element( "Titanium", "Ti", 22., 47.88*g/mole );
  new G4Element( "Iron", "Fe", 26., 55.847*g/mole );
  new G4Element( "Sulfur", "S", 16., 32.064*g/mole );
  new G4Element( "Phosphore", "P", 15., 30.974*g/mole );
  new G4Element( "Chlorine", "Cl", 17., 35.453*g/mole );
  new G4Element( "Chromium", "Cr", 24., 51.996*g/mole );
  new G4Element( "Manganese", "Mn", 25,54.93805*g/mole );

  DefineDefaultGeometry();
  //Caution the messenger should be defined after the DefineDefaultGeometry()
  detectorMessenger = new PLANETOCOSGeometryMessenger(this);
  
  //vacuum definition 
  //-------------------

  Vacuum =new G4Material("Vacuum",1.e-25*g/cm3,
                          1,kStateSolid, 1.e-10*kelvin, 1.e-17*atmosphere );
  Vacuum->AddElement(G4Element::GetElement("Hydrogen"), 1.);
 

  atmoSD = new PLANETOCOSSD("atmoSD",this);
  soilSD = new PLANETOCOSSoilSD("soilSD",this);
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->AddNewDetector(atmoSD);
  SDman->AddNewDetector(soilSD);
  geometry_verbosity =0;
  
  DetectionInAtmosphereInMiddleOfALayer = false;
  detect_below_soil_layers =false; 

}
////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSGeometryConstruction::~PLANETOCOSGeometryConstruction()
{ delete detectorMessenger;
  if (fMagneticField) delete  fMagneticField;
  if (theAtmosphericModel) delete theAtmosphericModel;
   
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* PLANETOCOSGeometryConstruction::Construct()
{ return ReConstructGeometry();
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::SetAtmosphereMaxStepLength(G4double maxStep)
{//G4cout<<maxStep/km<<std::endl; 
 
 theAtmosphereUserLimits->SetMaxAllowedStep(maxStep);
 AtmosphereMaxStepLength=maxStep;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::SetMagnetosphereMaxStepLength(G4double maxStep)
{ //G4cout<<maxStep/km<<std::endl; 
  theMagnetosphereUserLimits->SetMaxAllowedStep(maxStep);
  MagnetosphereMaxStepLength=maxStep;
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::SetMagnetosphereMaxTrackLength(G4double maxTrack)
{ theMagnetosphereUserLimits->SetUserMaxTrackLength(maxTrack);
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::SetMagnetosphereMaxTrackDuration(G4double maxTime)
{ theMagnetosphereUserLimits->SetUserMaxTime(maxTime);
}
////////////////////////////////////////////////////////////////////////////////
//
double PLANETOCOSGeometryConstruction::GetAtmosphereMaxStepLength( )
{
 //G4cout<<maxStep/km<<std::endl; 
 return AtmosphereMaxStepLength;
}
////////////////////////////////////////////////////////////////////////////////
//
double PLANETOCOSGeometryConstruction::GetMagnetosphereMaxStepLength( )
{
 //G4cout<<maxStep/km<<std::endl; 
 return MagnetosphereMaxStepLength;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* PLANETOCOSGeometryConstruction::ReConstructGeometry()


{
  // delete all solids, logical volumes and physical volumes if any 
  
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean(); 
    
  AtmoMaterials.clear();
  logicAtmosphere.clear();
  
  //build the geometry
  //------------------- 
  
  //Initialisation 
  
   G4bool detection_at_soil = false;
   PlanetManager* thePlanetManager = PlanetManager::GetInstance();
   G4double rplanet=thePlanetManager->GetRplanet();
   std::vector< double>   atmo_detector_altitudes;
   std::vector< double>   magneto_detector_altitudes;
   std::vector< double>  *pvector;

   
 
   
  
  //sort the detector_altitudes 
  
    for (unsigned int i=0; i<detector_altitudes.size();i++){
  	G4double altitude = detector_altitudes[i];
	if (altitude <= atmosphere_hmax && ConsiderAtmosphere && theAtmosphericModel) {
		pvector = &atmo_detector_altitudes;
	}
	else {
		pvector = &magneto_detector_altitudes;
		if (altitude <=atmosphere_hmin) {
			detection_at_soil =true;
			pvector=0;
		}	
	} 
	
	if (pvector) { 
		if (pvector->size()>0){
			if (altitude > pvector->front()){
	            		pvector->insert(pvector->begin(),
		                        	altitude);
			}					  
	  		else if (altitude < pvector->back()){
	        			pvector->push_back(altitude);
			}	   
	  		else{  
				G4bool out;
	      	       	       	G4int id =myfunc::locate(*pvector, altitude,out); 
	      	       	       	pvector->insert(pvector->begin()+id+1,altitude);
	      
	     		}
		}
        	else {   
			if (detector_altitudes[i]>0. || pvector ==  &atmo_detector_altitudes)
	  				pvector->push_back(detector_altitudes[i]);
		
      		}
	}	
	
  }
  pvector =0;

  std::vector<int > AtmoIndexDetectionBoundaries;
  std::vector<int > AtmoIndexDetectionLayers;
  bool atmo_detect_at_ground, atmo_detect_at_top;
  
  
  //Atmosphere structure
   
  if  (theAtmosphericModel && ConsiderAtmosphere ) {
  	//Get atmospheric layer structure from atmosphere model
  	//------------------------------------------------------

	if (!DetectionInAtmosphereInMiddleOfALayer){
		theAtmosphericModel->ComputeAtmosphericLayers(atmosphere_hmin,atmosphere_hmax,
                                                depth_percent,min_thickness,
				                max_thickness,
				                atmo_detector_altitudes,
				                detector_depths,
                                                AtmoAltitudes,AtmoDepths,
				                AtmoIndexDetectionBoundaries,  
				                AtmoMaterials);
	
	}
	else {
		theAtmosphericModel->ComputeAtmosphericLayers1(atmosphere_hmin,atmosphere_hmax,
                                                depth_percent,min_thickness,
				                max_thickness,
				                atmo_detector_altitudes,
				                detector_depths,
                                                AtmoAltitudes,AtmoDepths,
						AtmoIndexDetectionLayers,
					  	atmo_detect_at_top,
					   	atmo_detect_at_ground,
				                AtmoMaterials);
	
	}					

	if (geometry_verbosity >0){
  		G4cout<<"Structure of the Atmosphere "<<std::endl;
        	G4cout<<"name"<<'\t'<<"density[g/cm3]"
	      		<<'\t'<<"top altitude[km]"<<'\t'<<"bottom altitude[km]"
	      		<<'\t'<<"top depth[g/cm2]"<<'\t'<<"bottom depth[g/cm2]"
	      		<< std::endl;
		G4cout<<std::setiosflags(std::ios::scientific);
  		G4cout<<std::setprecision(6);      
  		for (unsigned int i=0; i<AtmoMaterials.size() ;i++){
       			G4cout<<AtmoMaterials[i]->GetName()<<'\t';
       			G4cout<<AtmoMaterials[i]->GetDensity()*cm3/g<<'\t';
       			G4cout<<AtmoAltitudes[i]/km<<'\t'<<AtmoAltitudes[i+1]/km<<'\t';
       			G4cout<<AtmoDepths[i]*cm2/g<<'\t'<<AtmoDepths[i+1]*cm2/g<<std::endl;
  		}
		G4cout<<std::setiosflags(std::ios::fixed);
  		G4cout<<std::setprecision(6);
	}		
  }
 	
      			           	 			    
  PLANETOCOSAnalysisManager::GetInstance()->GetEdepAnalyser()->SetAtmosphericLayerAltitudes(&AtmoAltitudes);
  PLANETOCOSAnalysisManager::GetInstance()->GetEdepAnalyser()->SetAtmosphericLayerDepths(&AtmoDepths); 

  //vis attributes
  //---------------- 
  G4VisAttributes * VisAttRed = new G4VisAttributes(G4Colour(1.,0.,0.));
  VisAttRed->SetVisibility(true);		  
  
  G4VisAttributes * VisAttGreen = new G4VisAttributes(G4Colour(0.,1.,0.));
  VisAttGreen->SetVisibility(true);			 
  
  G4VisAttributes * VisAttBlue = new G4VisAttributes(G4Colour(0.,0.,1.));
  VisAttBlue->SetVisibility(true);
  
  G4VisAttributes * VisAttRedGreen = new G4VisAttributes(G4Colour(1.,1.,0.));
  VisAttBlue->SetVisibility(true);
  
  /*G4VisAttributes * VisAttGreenBlue = new G4VisAttributes(G4Colour(0.,1.,1.));
  VisAttBlue->SetVisibility(true);*/
  
  
  
  G4VisAttributes * VisAttInvisible = new G4VisAttributes();
  VisAttInvisible->SetVisibility(false);
  
  //Soil attributes 
  //----
 
  PlanetSoil* theSoil = thePlanetManager->GetSoil();
  std::vector<G4Material*> theSoilMaterials = theSoil->GetSoilMaterials();
  std::vector<G4int> theSoilNbSubLayers = theSoil->GetSoilNbSubLayers();
  G4double SoilThickness =0.;
  std::vector<G4double> theSoilThicknessVector = theSoil->GetSoilThickness();
 
  for (unsigned int i=0; i< theSoilThicknessVector.size();i++){
  	SoilThickness +=theSoilThicknessVector[i];
  }
  
  
  //------------------------------ 
  // Magnetosphere  = world
  //------------------------------ 
  G4double WorldHeight= MagnetosphereH+PlanetH+SoilThickness;
  G4double ground_altitude =atmosphere_hmin;
  G4double top_atmosphere_altitude = ground_altitude;
  
  if  (theAtmosphericModel && ConsiderAtmosphere) { 
     ground_altitude = AtmoAltitudes.back();
     top_atmosphere_altitude = AtmoAltitudes[0];
     WorldHeight += (top_atmosphere_altitude-ground_altitude);
     //G4cout<<"WorldHeight "<<WorldHeight/km<<std::endl; 
  } 
 
     
  if (GeometryType =="SPHERICAL") {
  	fMagneticField->SetAltitudeWorldCenter(-rplanet);
  	fMagneticField->SetRotatePosAndField(false);
  }
  else {
  	fMagneticField->SetAltitudeWorldCenter( ground_altitude+ WorldHeight/2. - PlanetH - SoilThickness);
  	fMagneticField->SetRotatePosAndField(true);
  	fMagneticField->ComputeOrientation();
  }	
  	
  theMagnetosphereUserLimits = new G4UserLimits (.1*rplanet);
  theMagnetosphereUserLimits->SetUserMaxTrackLength(120. *rplanet);
  theMagnetosphereUserLimits->SetUserMaxTime(1000000000. *s);
    
  
  if  (GeometryType =="SPHERICAL"){ 
        G4Orb* solidMagnetosphere = new G4Orb("Magnetosphere", rplanet+MagnetosphereH+top_atmosphere_altitude);
  	/*G4double half_d= rplanet+MagnetosphereH+top_atmosphere_altitude;
	G4Box* solidMagnetosphere =new G4Box("Magnetosphere",half_d, half_d, half_d);*/
	
	logicMagnetosphere= new G4LogicalVolume(solidMagnetosphere, 
	                                  	Vacuum,
                                          	"Magnetosphere", 
			                  	0, 0, 0);
  }
  else  { 
  	G4Box* solidMagnetosphere = new G4Box("Magnetosphere",half_length*1.,half_length*1.,WorldHeight/2.);
	logicMagnetosphere= new G4LogicalVolume(solidMagnetosphere, 
	                                  	Vacuum,
                                          	"Magnetosphere", 
			                  	0, 0, 0);
  } 
  
  logicMagnetosphere->SetVisAttributes(VisAttInvisible);
  logicMagnetosphere->SetUserLimits(theMagnetosphereUserLimits);
   
  G4VPhysicalVolume* physiMagnetosphere= new G4PVPlacement(0,               
                                          		   G4ThreeVector(), 
				          		   "MagnetospherePV",   
                                          		   logicMagnetosphere,     
                                          		   0,              
                                          		   false,                                                	  	   0);             
  if  (GeometryType =="SPHERICAL") {
        ZpositionForAtmosphereBottom = rplanet + ground_altitude;
  }	
  else 
  	ZpositionForAtmosphereBottom =(-WorldHeight/2)+PlanetH+SoilThickness;	
   
  ZpositionAtZeroAltitude=ZpositionForAtmosphereBottom-ground_altitude;
   //G4cout<<" ZpositionForAtmosphereBottom "<<ZpositionForAtmosphereBottom/km<<std::endl;		 
   //G4cout<<" ZpositionAtZeroAltitude "<<ZpositionAtZeroAltitude/km<<std::endl;		 
  
  //altitudes and depth of detector for histograms information 
  detector_depths_for_histo.clear();
  detector_altitudes_for_histo.clear();
  
  
  
  //------------------------
  //Detection layer above atmosphere
  //-------------------------
  
  UpBoundaryIsADetector.clear();
  LowBoundaryIsADetector.clear();
  IsAnAtmoDetectionLayer.clear();
  unsigned int nmagneto_layers = magneto_detector_altitudes.size();
  unsigned int natmo_layers = AtmoMaterials.size();
  unsigned int nlayers_above_soil =nmagneto_layers +natmo_layers;
  G4cout<<nlayers_above_soil<<std::endl;
  unsigned int n_det_below_soil =0;
  if (detect_below_soil_layers){
   	for (unsigned int i=0; i< theSoilNbSubLayers.size();i++){
		n_det_below_soil+=theSoilNbSubLayers[i];
	}
  }
  G4cout<<n_det_below_soil<<std::endl;
  
  if (nlayers_above_soil>0) {
  	UpBoundaryIsADetector.insert(UpBoundaryIsADetector.end(),nlayers_above_soil+n_det_below_soil,-1);
  	LowBoundaryIsADetector.insert(LowBoundaryIsADetector.end(),nlayers_above_soil+n_det_below_soil,-1);
	IsAnAtmoDetectionLayer.insert(IsAnAtmoDetectionLayer.end(),nlayers_above_soil+n_det_below_soil,-1);
  }
 
  
  G4double low_altitude =top_atmosphere_altitude;
  //G4double crossing_factor = 1.01;
  G4double dlayer_for_crossing=0.01*mm;
  for ( int i=int(magneto_detector_altitudes.size())-1;i>=0;i+=-1){
  	G4double up_altitude=magneto_detector_altitudes[i];
	UpBoundaryIsADetector[i]=i;
	LowBoundaryIsADetector[i]=i+1;
	if (i == int( magneto_detector_altitudes.size())-1 )
				LowBoundaryIsADetector[i]=-1;	
	G4VSolid* detection_layer;	
	G4ThreeVector pos;
	if (GeometryType =="SPHERICAL") {
		pos= G4ThreeVector(0.,0.,0.);
		G4Orb* solid_layer1 = new  G4Orb("Detection",rplanet+up_altitude);
	    	G4Orb* solid_layer2 = new  G4Orb("Detection",rplanet+low_altitude);
	    	detection_layer = new G4SubtractionSolid("DetectionLayer",
							  solid_layer1,
			   				  solid_layer2);
	}
	else {  
		pos= G4ThreeVector(0.,0.,(low_altitude+up_altitude)/2.+ ZpositionAtZeroAltitude
								      );
		detection_layer = new G4Box("DetectionLayer",
	                             	     half_length,
			             	     half_length,
	                	             (up_altitude-low_altitude)/2. + dlayer_for_crossing);
	}
	
	std::stringstream astream;
	G4String str_i;
	astream<<i+1;
	astream>>str_i;
	G4LogicalVolume* logic_layer = new G4LogicalVolume(detection_layer ,
						          Vacuum,
                                	                  "DetectionLayer"+str_i,  0, 0, 0);
	logic_layer->SetUserLimits(theMagnetosphereUserLimits);					    
	if (GeometryType !="SPHERICAL")
	      logic_layer->SetVisAttributes(VisAttRedGreen);
	else   logic_layer->SetVisAttributes(VisAttInvisible);
	
	G4VPhysicalVolume* physic_layer;
	physic_layer = new G4PVPlacement(0,             
                                         pos, 
				         "DetectionLayer"+str_i,      
                                         logic_layer,    
                                         physiMagnetosphere,       
                                         false,          
                                         i);
	low_altitude = up_altitude;				 
  }
  
  if (detection_at_soil) LowBoundaryIsADetector.back()=int(magneto_detector_altitudes.size());
  
  for (unsigned int i=0;i<magneto_detector_altitudes.size();i++){
  	detector_depths_for_histo.push_back(0.);
	detector_altitudes_for_histo.push_back(magneto_detector_altitudes[i]);
  }
  
  
  if  (detection_at_soil) {
  	detector_depths_for_histo.push_back(0.); 
	detector_altitudes_for_histo.push_back(ground_altitude);
  }					  
 
 
  
	   
 // Atmosphere construction
 //--------------------------
  theAtmosphereUserLimits=new G4UserLimits (1000.*km);
  if  (theAtmosphericModel && ConsiderAtmosphere){
  
  	
   
 	//boundary detector
 	//-----------------
  	unsigned int nlayers=AtmoMaterials.size();
  	
        if (!DetectionInAtmosphereInMiddleOfALayer){
  		for (int i=0; i < (int) AtmoIndexDetectionBoundaries.size(); i++){
   		
		
		
			unsigned int index = AtmoIndexDetectionBoundaries[i];
        		detector_depths_for_histo.push_back(AtmoDepths[index]);
			detector_altitudes_for_histo.push_back(AtmoAltitudes[index]);
		
		
			unsigned int index1 = AtmoIndexDetectionBoundaries[i]+nmagneto_layers;
			int ii=  i+nmagneto_layers;
			if ( index>0 && index1 <nlayers){
	  			UpBoundaryIsADetector[index1]=ii;
	   			LowBoundaryIsADetector[index1-1]=ii;
			}   
			else if (index == 0) {
				UpBoundaryIsADetector[index1]=ii;
				if (nmagneto_layers >0 )
					LowBoundaryIsADetector[index1-1]=ii;
			}	
			else if (index  == nlayers) 
				LowBoundaryIsADetector[index1-1]=ii;
			        if (detect_below_soil_layers && theSoilMaterials.size()>0) UpBoundaryIsADetector[index1]=ii;;
		}		
  	}
	else {
		for (int i=0; i < (int) AtmoIndexDetectionLayers.size(); i++){
   		
		
		
			unsigned int index = AtmoIndexDetectionLayers[i]; 
			IsAnAtmoDetectionLayer[index+nmagneto_layers]=i+nmagneto_layers;
        		detector_depths_for_histo.push_back((AtmoDepths[index]+AtmoDepths[index+1])/2.);
			detector_altitudes_for_histo.push_back((AtmoAltitudes[index]+AtmoAltitudes[index+1])/2.);
			G4cout<<AtmoAltitudes[index]/km<<'\t'<<AtmoAltitudes[index+1]/km<<std::endl;
		
		
			
		}
		unsigned int index1 = nmagneto_layers;
		int ii=  nmagneto_layers;
		if (atmo_detect_at_top) {
			UpBoundaryIsADetector[index1]=ii;
			if (nmagneto_layers >0 )
				LowBoundaryIsADetector[index1-1]=ii;
		}
		
	
	}
	
	
	
  	// Sensitive Geometry
  	//-------------------
   
  
  	atmoSD->SetpAltitudes(&AtmoAltitudes);
  	atmoSD->SetpDepths(&AtmoDepths); 
	if  (GeometryType =="SPHERICAL") atmoSD->SetZSeaLevel(rplanet);
	else atmoSD->SetZSeaLevel((-WorldHeight/2)+PlanetH+SoilThickness-ground_altitude);
  	atmoSD->SetGeometryType(GeometryType);
    	atmoSD->SetNbMagnetoLayers(nmagneto_layers);
 
  			       
  	if (AtmoMaterials.size() > 0){
  		for (unsigned int  i=0;i<AtmoMaterials.size();i++){
           		G4double Alt1=AtmoAltitudes[i+1];
            		G4double Alt2=AtmoAltitudes[i];
        		G4ThreeVector pos;
			G4VSolid* solid_layer;
			if (GeometryType =="SPHERICAL") {			
				pos= G4ThreeVector(0.,0.,0.);
	 			G4Orb* solid_layer1 = new  G4Orb("Atmosphere",rplanet+Alt2);
	    			G4Orb* solid_layer2 = new  G4Orb("Atmosphere",rplanet+Alt1);
	    			solid_layer = 
	                   		new G4SubtractionSolid("Atmosphere",
								solid_layer1,
			   					solid_layer2);        
			}
			else {
				pos= G4ThreeVector(0.,0.,(Alt1+Alt2)/2.
				                          + ZpositionAtZeroAltitude
							  );
				solid_layer = new  G4Box("Atmosphere",
	                             			half_length,
			            	 		half_length,
	                             			(Alt2-Alt1)/2.+dlayer_for_crossing);
        		}
			std::stringstream astream;
	    		G4String str_i;
	    		astream<<i+1;
	    		astream>>str_i;
			G4String name="AtmosphereLayer"+str_i;
			if ( IsAnAtmoDetectionLayer[i] >=0) name="AtmosphereLayerDet"+str_i;
			
			logicAtmosphere.push_back(new G4LogicalVolume(solid_layer,
	                                                      	AtmoMaterials[i],
                                                              	name,
						              	0, 0, 0));
							
	    		if (GeometryType !="SPHERICAL")
			    logicAtmosphere[i]->SetVisAttributes(VisAttRed);
			else 
			    logicAtmosphere[i]->SetVisAttributes(VisAttInvisible);
	   		logicAtmosphere[i]->SetUserLimits(theAtmosphereUserLimits);
			if ( IsAnAtmoDetectionLayer[i] >=0) {
				logicAtmosphere[i]->SetUserLimits(new G4UserLimits ( (Alt2-Alt1)/5));
			}	
	    		logicAtmosphere[i]->SetSensitiveDetector(atmoSD);
	 
	    		G4VPhysicalVolume* physic_layer;
	    		physic_layer =  new G4PVPlacement(0,               
                                              	  	pos, 
				                  	name,       
                                                  	logicAtmosphere[i],      
                                                  	physiMagnetosphere,              
                                                  	false,         
                                                  	i+nmagneto_layers);  					  
		}  
  	}
   }
  
  
  
  //Soil construction
  //---------------- 
   unsigned int ndet_outside_soil =detector_altitudes_for_histo.size();
   SoilDepthsInLengthUnit.clear();
   SoilDepthsInLengthUnit.push_back(0.);
   SoilDepthsInDepthUnit.clear();
   SoilDepthsInDepthUnit.push_back(0.);
 
   G4double depth=0.;
   G4double cumulative_soil_thickness=0.;
   G4double total_atmosphere_depth= 0.;
   
   if  (theAtmosphericModel && ConsiderAtmosphere) { 
      total_atmosphere_depth= AtmoDepths.back();
   } 	 
   G4double soil_layer_topZ =ZpositionForAtmosphereBottom;
   G4int ind_soil_layer=-1;
   
  // Sensitive Detector for soil
  //-------------------
   
   
   
   PLANETOCOSAnalysisManager::GetInstance()->GetSoilEdepAnalyser()->SetSoilLayerDepthsInLengthUnit(&SoilDepthsInLengthUnit);
   PLANETOCOSAnalysisManager::GetInstance()->GetSoilEdepAnalyser()->SetSoilLayerDepthsInDepthUnit(&SoilDepthsInDepthUnit); 
    
   
   soilSD->SetZGroundLevel(ZpositionForAtmosphereBottom);
   soilSD->SetGeometryType(GeometryType);
   soilSD->SetNbLayersAboveSoil(nmagneto_layers+natmo_layers);
   
   
   
   if (theSoilThicknessVector.size()>0){
   	G4FieldManager* soilFManager = new G4FieldManager(0,0,true);
	soilSD->SetpDepthsInDepthUnit(&SoilDepthsInDepthUnit);
   	soilSD->SetpDepthsInLengthUnit(&SoilDepthsInLengthUnit);
	
	for (unsigned int i=0; i< theSoilMaterials.size();i++){
		G4double thickness = theSoilThicknessVector[i]/theSoilNbSubLayers[i];
		for (unsigned int j=0; j<theSoilNbSubLayers[i];j++){
			cumulative_soil_thickness+=thickness;
			G4ThreeVector pos;
			G4VSolid* solid_soil_layer;
			
			ind_soil_layer+=1;
			std::stringstream astream;
             		G4String str_i,str_j;
  			astream<<i+1<<'\t'<<j+1;
  			astream>>str_i>>str_j;
			
			
			
  			G4String name="Soil"+str_i+"_"+str_j;
		
		
			if (GeometryType =="SPHERICAL") {
				pos = G4ThreeVector(0.,0.,0.);
				G4Orb* solid_layer1 = new  G4Orb("Soil",soil_layer_topZ);
	    			G4Orb* solid_layer2 = new  G4Orb("",soil_layer_topZ-thickness);
	    			solid_soil_layer = 
	                   			new G4SubtractionSolid("Soil",
									solid_layer1,
			   						solid_layer2);  
			}
			else {  pos = G4ThreeVector(0.,0.,soil_layer_topZ-thickness/2.);
	    			solid_soil_layer = new  G4Box("Soil",
	                             	 			half_length,
			            	 			half_length,
	                             	 			thickness/2.);
			}					
			G4LogicalVolume* logic_soil_layer;
			logic_soil_layer = new G4LogicalVolume(solid_soil_layer,
	                                               	theSoilMaterials[i],
                                                       	theSoilMaterials[i]->GetName(),  0, 0, 0);
			logic_soil_layer->SetVisAttributes(VisAttInvisible);
			logic_soil_layer->SetUserLimits(theAtmosphereUserLimits);
			logic_soil_layer->SetFieldManager(soilFManager,false);
			logic_soil_layer->SetSensitiveDetector(soilSD);	
			
			if (detect_below_soil_layers) name=G4String("Layer")+name;		     
	        	G4VPhysicalVolume* physic_layer;
	    			physic_layer =  new G4PVPlacement(0,               
                                              	  		pos, 
				                  		name,       
                                                  		logic_soil_layer,      
                                                  		physiMagnetosphere,              
                                                  		false,         
                                                  		ind_soil_layer+nmagneto_layers+natmo_layers);
			soil_layer_topZ+= -thickness;
			G4double alt=ground_altitude-cumulative_soil_thickness;
			depth += theSoilMaterials[i]->GetDensity()*thickness;
			SoilDepthsInLengthUnit.push_back(cumulative_soil_thickness);
			SoilDepthsInDepthUnit.push_back(depth);
			
			if (detect_below_soil_layers) {
				
				detector_depths_for_histo.push_back(depth+total_atmosphere_depth);
				detector_altitudes_for_histo.push_back(alt);
			
				unsigned int  index=nmagneto_layers+natmo_layers+ind_soil_layer;
				unsigned int ndet_below=ndet_outside_soil+ind_soil_layer;
				G4cout<<"TEST ndet_below"<<ind_soil_layer<<'\t'<<ndet_below<<'\t'<<detection_at_soil<<std::endl;
				G4cout<<LowBoundaryIsADetector[0]<<std::endl;
				LowBoundaryIsADetector[index]=ndet_below;
				G4cout<<"TEST1"<<std::endl;
				if (  ind_soil_layer> 0 || detection_at_soil){
	  				UpBoundaryIsADetector[index]=ndet_below-1;
	   		
				}   
				else {
					UpBoundaryIsADetector[index]=-1;
				}
				G4cout<<"TEST2"<<std::endl;
			
			
			
			}
		}
	}	
  
   } 
   
   
   int ndet_layers = int(detector_altitudes_for_histo.size());
   G4cout<<"n_layers"<<ndet_layers<<std::endl;
 

   //output of detector altitudes added
   //AM
   for (int alt=0; alt<detector_altitudes_for_histo.size(); ++alt) {
   		G4cout<<"Detector "<<alt+1<<" altitude: "<<detector_altitudes_for_histo[alt]/km<<" km"<<std::endl;
   		}
   	








  if (ndet_layers  >0){
   	PLANETOCOSFluxDetectionAnalyser* theFluxDetectionAnalyser =
			PLANETOCOSAnalysisManager::GetInstance()->GetFluxDetectionAnalyser();	
  	theFluxDetectionAnalyser->InitialiseHistograms(ndet_layers);
	theFluxDetectionAnalyser->SetDetectorAltitudes(&detector_altitudes_for_histo);		   
   	theFluxDetectionAnalyser->SetDetectorDepths(&detector_depths_for_histo);		   
   	PLANETOCOSAnalysisManager::GetInstance()->SetDetectorAltitudes(&detector_altitudes_for_histo);
	PLANETOCOSAnalysisManager::GetInstance()->SetDetectorDepths(&detector_depths_for_histo);
       
        PLANETOCOSSteppingAction* theStepAction =
   			dynamic_cast< PLANETOCOSSteppingAction* >
   			  	(const_cast< G4UserSteppingAction* >(
				    			G4RunManager::GetRunManager()
									->GetUserSteppingAction()));
	theStepAction->SetpAltitudes(&detector_altitudes_for_histo);								
	
   }
   else {
        PLANETOCOSAnalysisManager::GetInstance()->SetDetectorAltitudes(0);
	PLANETOCOSAnalysisManager::GetInstance()->SetDetectorDepths(0);
   }

   	    
  //----------------------
  //Planet
  //-----------------------
  G4String PlanetName = PlanetManager::GetInstance()->GetPlanetName(); 
  G4ThreeVector pos;
  G4VSolid* solidPlanet;
  if (GeometryType =="SPHERICAL") {
	pos = G4ThreeVector(0.,0.,0.);
	solidPlanet=new G4Orb(PlanetName,soil_layer_topZ);
  }
  else {pos = G4ThreeVector(0.,0.,(-PlanetH/2.)+soil_layer_topZ);
   	solidPlanet= new G4Box(PlanetName,half_length,half_length, PlanetH/2.);
  }
  G4LogicalVolume* logicPlanet = new G4LogicalVolume(solidPlanet,
  						    Vacuum,
						    PlanetName,
						    0,0,0);
  logicPlanet->SetVisAttributes(VisAttBlue);
  logicPlanet->SetUserLimits(theAtmosphereUserLimits);
  G4VPhysicalVolume* physiPlanet;
  physiPlanet =new G4PVPlacement(0,pos, PlanetName,        
                                logicPlanet,       
                                physiMagnetosphere,        
                                false,         
                                0);
  
  //Set type of geoemtry to planet manager
  //------------------------------------
  if (GeometryType =="SPHERICAL"){
  	thePlanetManager->GetMagneticField()->SetIsFlatGeometry(false);
  	SetMagnetosphereMaxStepLength(rplanet);
	SetAtmosphereMaxStepLength(5.*km);
  
  }		
  else {
  	thePlanetManager->GetMagneticField()->SetIsFlatGeometry(true);
	SetMagnetosphereMaxStepLength(5.*km);
	SetAtmosphereMaxStepLength(5.*km);
  }

  return physiMagnetosphere;

  
} 
////////////////////////////////////////////////////////////////////////////////
// 
void PLANETOCOSGeometryConstruction::UpdateGeometry()
{ G4RunManager::GetRunManager()->SetGeometryToBeOptimized(true);
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4RunManager::GetRunManager()->DefineWorldVolume(ReConstructGeometry());
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::DefineDefaultGeometry()
{ 
  PlanetManager* thePlanetManager = PlanetManager::GetInstance();
  G4String planet_name = thePlanetManager->GetPlanetName();
  
  //default value of geometrical parameter
  //--------------
  
  
  SetDepthPercent(5.);
  GeometryType ="SPHERICAL";
  MagnetosphereH = 60 *thePlanetManager->GetRplanet();
  atmosphere_hmax = 100.*km;
  atmosphere_hmin = 0.*km;
  detector_altitudes.clear();
  detector_depths.clear();
  max_thickness =5.*km;
  min_thickness =.0001*km;
  half_length =500*km;
  PlanetH = 10*km;

  //Atmosphere
  //----------
  theAtmosphericModel =thePlanetManager->GetAtmosphereModel();
  theAtmosphericModel->SetGeometry(this);
  if (!theAtmosphericModel) {
  	ConsiderAtmosphere = false;
  }
  else {
  	theAtmosphericModel->SetDefault();
	ConsiderAtmosphere = true;
	if (planet_name == "Mercury") ConsiderAtmosphere = false;
  }

  //Soil
  //----
  PlanetSoil* theSoil = thePlanetManager->GetSoil();
  theSoil->SetDefault();	
  
  //Switch off magnetic field
  //-------------------------
  fMagneticField = thePlanetManager->GetMagneticField();
  if (fMagneticField)  fMagneticField->SwitchOff(); 

  
  
  
}

////////////////////////////////////////////////////////////////////////////////
//
G4int PLANETOCOSGeometryConstruction::GetNbOfLayers()
{ return int (logicAtmosphere.size());
}
////////////////////////////////////////////////////////////////////////////////
//
G4LogicalVolume* PLANETOCOSGeometryConstruction::GetLogicalLayer(G4int i)        
{ if (i< int (logicAtmosphere.size()) ) return  logicAtmosphere[i];
  else return 0;	  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::CreateElement(G4String el_name,G4String el_symbol,G4double z,G4double a)
{ if (!G4Element::GetElement(el_name) && !G4Element::GetElement(el_symbol)) {
  	new G4Element(el_name, el_symbol, z, a);
		
  }
  else if (G4Element::GetElement(el_name)){
  	G4cout<<"An element with name "<<el_name<<" is already defined."<<std::endl;
  }
  else{
  	G4cout<<"An element with symbol "<<el_symbol<<" is already defined."<<std::endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSGeometryConstruction::ListElements()
{ const G4ElementTable* theElementTable = G4Element::GetElementTable();
  for (unsigned int i =0;i< theElementTable->size(); i++){
  	G4cout<<(*theElementTable)[i]->GetName()<<'\t';
	G4cout<<(*theElementTable)[i]->GetSymbol()<<'\t';
	G4cout<<(*theElementTable)[i]->GetZ()<<'\t';
	G4cout<<(*theElementTable)[i]->GetA()<<std::endl;
  
  }
	
}
