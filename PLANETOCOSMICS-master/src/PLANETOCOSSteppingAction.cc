
#include "PLANETOCOSSteppingAction.hh"
#include "PLANETOCOSEventAction.hh"
#include "G4SteppingManager.hh"
#include "PlanetMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "PLANETOCOSSteppingActionMessenger.hh"
#include "G4SDManager.hh"
#include "PLANETOCOSSD.hh"
#include "PLANETOCOSTrackInformation.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "PLANETOCOSGeometryConstruction.hh"
#include "G4RunManager.hh"
#include "PlanetManager.hh"
#include"DurationManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
PLANETOCOSSteppingAction::PLANETOCOSSteppingAction()
{ UntrackedParticleDic=new PartDictionary();
  theMessenger = new PLANETOCOSSteppingActionMessenger(this);
  DetectFlux=false;
  last_trackID=100000000;
  StopAtMagnetopause =false;
  MagnetopauseOutFactor = 1.;
  LowestAltitudeNeeded =1e12*km;
  MaxNumberOfTurnAroundThePlanet=20000000000000000000000000000000000000000000000.;
  
  StopAltitudeForUpward = 1.e50*km;
  StopUpwardPrimary =false;
  
  StopDownwardPrimary =false;
  StopAltitudeForDownward = -1.e50*km;
  LastStepWasOnBoundary = false;
  
}
////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSteppingAction::UserSteppingAction(const G4Step* aStep)
{ DurationManager* theDurationManager = DurationManager::GetInstance();
  if (!theDurationManager->CheckDuration()) {
  	PLANETOCOSAnalysisManager::GetInstance()->SwitchOffRegisterResults();
  	G4cout<< "The event will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
  }
  const G4VPhysicalVolume* currentVolume = aStep->GetPostStepPoint()
                                                       ->GetPhysicalVolume();
  const G4VTouchable* preStepTouchable = aStep->GetPreStepPoint()->GetTouchable();
  
  G4Track* aTrack=aStep->GetTrack();
  G4String aParticleName=aTrack->GetDefinition()->GetParticleName();
  nStep =aTrack->GetCurrentStepNumber();
  if (nStep ==1) LastStepWasOnBoundary =false;
  G4int trackID =aTrack->GetTrackID();
  G4ThreeVector position   = aTrack->GetPosition();
  G4ThreeVector pre_position   = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector direction1   = position - pre_position;
  G4double energy =aTrack->GetKineticEnergy();
  PLANETOCOSGeometryConstruction* theGeometry =
   			dynamic_cast< PLANETOCOSGeometryConstruction* >
   			  	(const_cast< G4VUserDetectorConstruction* >(
				    G4RunManager::GetRunManager()
						->GetUserDetectorConstruction()));
  
  G4String geometry_type =  theGeometry->GetGeometryType();
  
  G4double altitude, pre_altitude;
  if ( geometry_type == "SPHERICAL"){
  	G4double rplanet = PlanetManager::GetInstance()->GetRplanet();
  	pre_altitude = pre_position.mag()-rplanet;
        altitude = position.mag()-rplanet;
  }
  else {
  	G4double ZposZeroAlt= theGeometry->GetZpositionAtZeroAltitude();
  	pre_altitude = pre_position.z()-ZposZeroAlt;
	altitude = position.z()-ZposZeroAlt;
  
  }



//  G4cout<<"test1"<<std::endl;
#ifdef ANALYSIS_SOIL
  G4double ZorR = position.z();
  if (geometry_type=="SPHERICAL") ZorR =position.mag();
  if (ZorR <LowestZorRReached) LowestZorRReached = ZorR;
#endif
  G4double half_length = theGeometry->GetHalfLength();
  PLANETOCOSTrackInformation* aTrackInfo 
  			= dynamic_cast<PLANETOCOSTrackInformation*>
				            (aTrack->GetUserInformation());
  
//  G4cout<<"test2"<<std::endl;
  if (geometry_type=="SPHERICAL" && currentVolume) {
  	
	if (pre_position.z()*position.z()<0) aTrackInfo->AddOneCrossingOfEquator(); 
   	G4ThreeVector vec1,vec2;
   	vec1=G4ThreeVector(pre_position.x(),pre_position.y(),0 );
   	vec2=G4ThreeVector(position.x(),position.y(),0 );
   	G4double dTurn= std::abs(vec1.angle(vec2))/2./3.1415926;
   	G4ThreeVector vec3 = vec1.cross(vec2);
   	if (vec3.z() <0) dTurn=-dTurn;
   	aTrackInfo->AddNbOfPlanetTurn(dTurn);
	
	G4String VName = currentVolume->GetName();
	if ( VName == "MagnetospherePV" || VName.contains("DetectionLayer"))
		     aTrackInfo->SetWasAlreadyInMagnetosphere(true);
					
	LastLowestAltitude = aTrackInfo->GetLowestAltNeededForThisTrack();
	//G4cout<<LastLowestAltitude<<std::endl;
	if (altitude < LastLowestAltitude){
		LastLowestAltitude = altitude;
		aTrackInfo->SetLowestAltNeededForThisTrack(LastLowestAltitude);
		
	}
  }
//  G4cout<<"test3"<<std::endl; 



//detection of particles 
//-------------------------
G4bool step_at_boundary = 
  		(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
LastStepWasOnBoundary =
 		(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary);

//detection at detection boundary
//------------------------------

if (  step_at_boundary &&
      (geometry_type == "SPHERICAL" ||
      (std::abs(position.x())< half_length && std::abs(position.y())<half_length))){
 
      
	
     	const G4VTouchable* postStepTouchable = aStep->GetPostStepPoint()->GetTouchable();
   	
	G4ThreeVector direction   = aTrack->GetMomentumDirection();
        if ( geometry_type == "SPHERICAL"){
	     			direction.rotateZ (-position.phi()).rotateY(-position.theta()); 
				direction1.rotateZ (-position.phi()).rotateY(-position.theta());
	}
  	
	G4double cos_th = std::abs(direction.z());
// do not detect and kill e- produced with small cos_theta on boundary or blocked on boundary
	
	G4bool detect = true;
	/*if ( aParticleName.contains("e-") ){
		if (nStep == 1 && cos_th <= 0.05) {
			detect = false;
			aTrack->SetTrackStatus(fStopAndKill);
			return;	
		}	
		else if (LastStepWasOnBoundary && cos_th <= 0.05){
			detect= false;
			aTrack->SetTrackStatus(fStopAndKill);
			return;	
			
		}
	
	}
	*/

	
	
	//G4cout<<"Detect flux"<<DetectFlux<<std::endl;
	//Flux detection if any 
  	//---------------------
	if (DetectFlux && postStepTouchable->GetVolume()){
	
		G4int UpDetBound1 = -1;
       		G4int LowDetBound1 = -1; 
      		G4int layer1 = -1;
       		if (preStepTouchable->GetVolume()->GetName().contains("Layer")){ 
              		layer1  = preStepTouchable->GetVolume()->GetCopyNo();
			UpDetBound1 = theGeometry->GetUpBoundaryIsADetector(layer1);
	       		LowDetBound1 = theGeometry->GetLowBoundaryIsADetector(layer1);
			//G4cout<<layer1<<'\t'<<UpDetBound1<<'\t'<<LowDetBound1<<std::endl;
			
	       	}
		/*G4cout<<preStepTouchable->GetVolume()->GetName()<<std::endl;
                G4cout<<"layer1 "<<layer1<<std::endl;
		G4cout<<"UpDetBound1 "<<UpDetBound1<<std::endl;
		G4cout<<"LowDetBound1 "<<LowDetBound1<<std::endl;*/
                G4int layer2 = -1;
		G4int UpDetBound2 = -1;
       		G4int LowDetBound2 = -1;
       	
		if (postStepTouchable->GetVolume()->GetName().contains("Layer")){ 
              		layer2  = postStepTouchable->GetVolume()->GetCopyNo();
			UpDetBound2 = theGeometry->GetUpBoundaryIsADetector(layer2);
	       		LowDetBound2 = theGeometry->GetLowBoundaryIsADetector(layer2);
		} 
		/*G4cout<<postStepTouchable->GetVolume()->GetName()<<std::endl;
                G4cout<<"layer2 "<<layer2<<std::endl;
		G4cout<<"UpDetBound2 "<<UpDetBound2<<std::endl;
		G4cout<<"LowDetBound2 "<<LowDetBound2<<std::endl;*/
        
		G4int DetectorBoundary;
        	DetectorBoundary = -1; 
	 	G4ThreeVector normal(0.,0.,0.);
		
		
		if (layer1 > -1 && layer2 >-1) {// in detection region
        		if (layer1 < layer2) { //moving downward
				DetectorBoundary = LowDetBound1;
				normal = G4ThreeVector(0.,0.,-1.);
			}	
	 		else { //moving upward
				DetectorBoundary = LowDetBound2;
				normal = G4ThreeVector(0.,0.,1.);
				if (geometry_type == "SPHERICAL"){
					if (LastLowestAltitude <LowestAltitudeNeeded)
						LowestAltitudeNeeded = LastLowestAltitude;
					if (!aTrackInfo->GetHasBeenAlreadyUpwardDetected()){
						aTrackInfo->SetHasBeenAlreadyUpwardDetected(true);
						aTrackInfo->SetDetectionEnergy(energy);
					}
				}	
			}
			//G4cout<<DetectorBoundary <<std::endl;	
		}
       		else if (layer1 >-1 ) {// leaving the detection region
        		if (postStepTouchable->GetVolume()->GetName() == "MagnetospherePV") { //from the top
			 	DetectorBoundary = UpDetBound1;
				normal = G4ThreeVector(0.,0.,1.);
				if (geometry_type == "SPHERICAL"){
					if (LastLowestAltitude <LowestAltitudeNeeded)
						LowestAltitudeNeeded = LastLowestAltitude;
					if (!aTrackInfo->GetHasBeenAlreadyUpwardDetected()){
						aTrackInfo->SetHasBeenAlreadyUpwardDetected(true);
						aTrackInfo->SetDetectionEnergy(energy);
					}
				}
			}	 
	  		else   {//from the bottom 
				DetectorBoundary = LowDetBound1;
				normal = G4ThreeVector(0.,0.,-1.);
			}	
		}
       		else if (layer2 > -1) {// entering the detection region
       			if (preStepTouchable->GetVolume()->GetName() == "MagnetospherePV") {//from the top
				DetectorBoundary = UpDetBound2;
				normal = G4ThreeVector(0.,0.,-1.);
			}	 
	  		else  { //from the bottom          
		 		DetectorBoundary = LowDetBound2;
				normal = G4ThreeVector(0.,0.,1.);
				if (geometry_type == "SPHERICAL"){
					if (LastLowestAltitude <LowestAltitudeNeeded)
						LowestAltitudeNeeded = LastLowestAltitude;
					if (!aTrackInfo->GetHasBeenAlreadyUpwardDetected()){
						aTrackInfo->SetHasBeenAlreadyUpwardDetected(true);
						aTrackInfo->SetDetectionEnergy(energy);
					}
				}
			}	 
		}
		
		//G4cout<<DetectorBoundary <<std::endl;	
       		// if DetectorBoundary differs from -1 a detector boundary has been
       		// reached and a hit flux will be stored
       		if (DetectorBoundary > -1 && aStep->GetStepLength() != 0.){
         		G4int PDGCode = aTrack->GetDefinition()->GetPDGEncoding();
	  		if (direction.dot(normal) >0){
				direction1=-direction1; 
				direction=-direction;   
	  			G4double theta = direction.theta();
				G4double theta1 = direction1.theta();
				/*if ( aParticleName.contains("e-") ){
					G4cout<<"Test "<<theta/degree<<'\t'<<theta1/degree<<std::endl;
				}
				*/	
				 
	  			G4double azimuth =  180.*degree -direction.phi();
				   
	  			G4double weight  = aTrack->GetWeight();
				G4double lat_or_y, long_or_x;
				if (geometry_type == "SPHERICAL"){
					lat_or_y = 90.*degree -position.theta();
					long_or_x = position.phi();
					if  (long_or_x>180.*degree) long_or_x = 
							long_or_x-360*degree; 
				}
				else {
					lat_or_y =position.y();
					long_or_x =position.x();
					
				}	
				
				PLANETOCOSFluxHit* aHit =new  PLANETOCOSFluxHit(weight, 
								      PDGCode,energy,
                                      				      azimuth, theta, 
								      lat_or_y,long_or_x,
								      DetectorBoundary);
					  			       
	  			PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                       			(G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
          			sd->GetDetectorFluxCollection()->insert(aHit);
				direction=-direction;
#ifdef TEST_ELEC_AT_BOUNDARY

				if ( aParticleName.contains("e-") ){
					//if (aStep->GetPostStepPoint()){
						//G4String ProcessName =aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
						if (aTrack->GetCreatorProcess()) {
						        G4cout<<"Some problems "<<aTrackInfo<<std::endl;
							G4int prevDetBound=aTrackInfo->GetNbLastDetector();
							G4double costhAtPrevDet=aTrackInfo->GetCosthAtLastDetector();
							G4int nstep_at_prev_det=aTrackInfo->GetNStepAtLastDetector();
							G4int nstep = aTrack->GetCurrentStepNumber();
							G4String ProcessName =aTrack->GetCreatorProcess()->GetProcessName();
							G4cout<<ProcessName<<std::endl;
							if (ProcessName == "conv") ProcessName="conve-";
							/*if (std::abs(std::cos(theta)) >0.2){
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),ProcessName,nstep);	
							}
							else if (prevDetBound != DetectorBoundary ) {
							
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),ProcessName,nstep);		
								
							}*/
							PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),"conve-",nstep);
							if (aStep->GetPreStepPoint()){
								G4ThreeVector direction1   = aStep->GetPreStepPoint()->GetMomentumDirection();
        							if ( geometry_type == "SPHERICAL")
	     								direction1.rotateZ (-position.phi()).rotateY(-position.theta());
  								direction1=-direction1;
								G4double theta1 = direction1.theta();
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta1),"conve-",nstep);
								
								
								
							}
						
							aTrackInfo->SetNbLastDetector(DetectorBoundary);
							aTrackInfo->SetNStepAtLastDetector(nstep);
							aTrackInfo->SetCosthAtLastDetector(std::cos(theta));
							
						}
						else {
							
						        G4cout<<"Some problems "<<aTrackInfo<<std::endl;
							G4int prevDetBound=aTrackInfo->GetNbLastDetector();
							G4double costhAtPrevDet=aTrackInfo->GetCosthAtLastDetector();
							G4int nstep_at_prev_det=aTrackInfo->GetNStepAtLastDetector();
							G4int nstep = aTrack->GetCurrentStepNumber();
							G4String ProcessName ="compt";
							G4cout<<ProcessName<<std::endl;
							if (ProcessName == "conv") ProcessName="conve-";
							/*if (std::abs(std::cos(theta)) >0.2){
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),ProcessName,nstep);	
							}
							else if (prevDetBound != DetectorBoundary ) {
							
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),ProcessName,nstep);		
								
							}*/
							
							PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),"conve-",nstep);
							if (aStep->GetPreStepPoint()){
								G4ThreeVector direction1   = aStep->GetPreStepPoint()->GetMomentumDirection();
        							if ( geometry_type == "SPHERICAL")
	     								direction1.rotateZ (-position.phi()).rotateY(-position.theta());
  								direction1=-direction1;
								G4double theta1 = direction1.theta();
								PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta1),"conve-",nstep);
								
								
								
							}
							aTrackInfo->SetNbLastDetector(DetectorBoundary);
							aTrackInfo->SetNStepAtLastDetector(nstep);
							aTrackInfo->SetCosthAtLastDetector(std::cos(theta));
						
						
						}	
					//}
	
				}
				if ( aParticleName.contains("e-") ){
					//if (aStep->GetPostStepPoint()){
						//G4String ProcessName =aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
						if (aTrack->GetCreatorProcess()) {
							G4int n_step = aTrack->GetCurrentStepNumber();
							G4String ProcessName =aTrack->GetCreatorProcess()->GetProcessName();
							if (ProcessName == "conv") ProcessName="conve+";
							//PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),"eIoni",n_step);
						}	
					//}
	
				}
					
		
#endif 
			}
			
			else {//test to remove the e- flux peak close to the horizontal
				  
	  			G4double theta = direction.theta();
				G4double theta1 = direction1.theta();
				//G4cout<<"Test1 "<<theta/degree<<'\t'<<theta1/degree<<std::endl;
				
				/*G4cout<<"Kill the particle "<<std::endl;
				
				//aTrack->SetTrackStatus(fStopAndKill);
				
				G4cout<<direction<<std::endl;
				G4cout<<normal<<std::endl;*/
			
			}	
	  
	 	}
	
	
	}
	

  }
  
  
  
//detection in  detection boundary
//------------------------------
  if (preStepTouchable->GetVolume() && DetectFlux){
   	if (preStepTouchable->GetVolume()->GetName().contains("AtmosphereLayerDet")){
		G4int layer1 = -1;
		layer1  = preStepTouchable->GetVolume()->GetCopyNo();
		G4int ndet = theGeometry->GetIsAnAtmosphericDetectionLayer(layer1);
		G4double alt_detection=  (*pAltitudes)[ndet];
		
		//G4cout<<ndet<<'\t'<<alt_detection/km<<'\t'<<pre_altitude/km<<'\t'<<altitude/km<<std::endl;
		if ((pre_altitude-alt_detection)*(altitude-alt_detection) <0){
			G4int PDGCode = aTrack->GetDefinition()->GetPDGEncoding();
			G4ThreeVector direction   = aTrack->GetMomentumDirection();
        		if ( geometry_type == "SPHERICAL")
	     			direction.rotateZ (-position.phi()).rotateY(-position.theta());
	  		direction=-direction;   
	  		G4double theta = direction.theta();
			G4double azimuth =  180.*degree -direction.phi();
			G4double weight  = aTrack->GetWeight();
			G4double lat_or_y, long_or_x;
			if (geometry_type == "SPHERICAL"){
				lat_or_y = 90.*degree -position.theta();
				long_or_x = position.phi();
				if  (long_or_x>180.*degree) long_or_x = 
							long_or_x-360*degree; 
			}
			else {
				lat_or_y =position.y();
				long_or_x =position.x();
					
			}
			//G4cout<<alt_detection/km<<'\t'<<pre_altitude/km<<'\t'<<altitude/km<<std::endl;
			
			PLANETOCOSFluxHit* aHit =new  PLANETOCOSFluxHit(weight, 
								      PDGCode,energy,
                                      				      azimuth, theta, 
								      lat_or_y,long_or_x,
								      ndet);
			PLANETOCOSSD* sd  = dynamic_cast<PLANETOCOSSD*>
                       			(G4SDManager::GetSDMpointer()->FindSensitiveDetector("atmoSD"));
          		sd->GetDetectorFluxCollection()->insert(aHit);
#ifdef TEST_ELEC_AT_BOUNDARY
			G4int nstep = aTrack->GetCurrentStepNumber();
			if ( aParticleName.contains("e-") ){
				
				PLANETOCOSAnalysisManager::GetInstance()->RegisterElecProducedAtBoundary(std::cos(theta),"eIoni",nstep);
				if (std::abs(std::cos(theta))<0.2){
					//G4cout<<std::cos(theta)<<'\t'<<nstep<<'\t'<<ndet<<std::endl; 
				
				}
			}
#endif				
		}						      
		
	
	}
   }
  

  //check if downward/upward particle should be stopped because crossing a given altitude 
  
  if (pre_altitude >= StopAltitudeForDownward &&   altitude<= StopAltitudeForDownward ){
  	if (aTrack->GetTrackID()!=1 || StopDownwardPrimary)
			aTrack->SetTrackStatus(fStopAndKill);
			//G4cout<<"STOP1"<<std::endl;
			return;	
  }
  if (altitude >= StopAltitudeForUpward &&   pre_altitude<= StopAltitudeForUpward ){
  	if (aTrack->GetTrackID()!=1 || StopUpwardPrimary)
			aTrack->SetTrackStatus(fStopAndKill);
			//G4cout<<"STOP2"<<std::endl;
			return;	
  }
  


  LastStepWasOnBoundary =step_at_boundary;




  //Some Stop particle condition 
  //-------------------------------	
  if (currentVolume ){
  	const G4String name = currentVolume->GetName();
    	//G4cout<<"Stepping2"<<std::endl; 
	//Stop at the planet core below soil
	//----------------------------------
	if (name == PlanetManager::GetInstance()->GetPlanetName())  {
#ifdef DETECT_AFTER_SOIL
		PLANETOCOSAnalysisManager::GetInstance()->DetectParticleAfterSoil(aStep);
#endif
    		//G4cout<<"STOP3"<<std::endl;
		aTrack->SetTrackStatus(fStopAndKill);
    		return;
    	}	
    	//G4cout<<"Stepping22"<<std::endl; 
	
	G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      
        //check if the particle is outside the magnetosphere
	//---------------------------------------------------
      	PlanetMagneticField* theField = PlanetManager::GetInstance()->GetMagneticField();
      	if (StopAtMagnetopause && theField){
      		if (!PrimaryAlreadyInMagnetosphere){
      			if (aTrack->GetTrackID()==1 && !theField->OutsideMagnetosphere(position)){
      				PrimaryAlreadyInMagnetosphere = true;	
					//G4cout<<"Primary is from now in the magnetosphere"<<std::endl;
			}	
      		}	
      		else {
      			if (theField->OutsideMagnetosphere(position/MagnetopauseOutFactor)){
           			aTrack->SetTrackStatus(fStopAndKill);
				//G4cout<<"STOP4"<<std::endl;
	    			return;
			}
      		}
      	} 				
			 
       //G4cout<<"Stepping222"<<std::endl; 
       //check the user limits
       //--------------------
      	G4UserLimits* theUserLimits = 0;
	if (currentVolume->GetLogicalVolume())
                  	theUserLimits = currentVolume->GetLogicalVolume()->GetUserLimits();			 
      	
	if (theUserLimits) {
	
		//stop the particle if the track length is greater than the maximum
      	        // allowed track length
      		
		G4double max_track_length =
                      	theUserLimits->GetUserMaxTrackLength(*aTrack);
      
      		if (aTrack->GetTrackLength() >= max_track_length){
   			aTrack->SetTrackStatus(fStopAndKill);
			//G4cout<<"STOP5"<<std::endl;
	    		return;
		} 
			 
			 
     		//stop the particle if the local time of the particle is greater than the maximum allowed
     		// local time of a trajectory
      
      		G4double max_local_time =
                     		theUserLimits->GetUserMaxTime(*aTrack);
      		if (aTrack->GetLocalTime() >= max_local_time){
     			aTrack->SetTrackStatus(fStopAndKill);
	    		G4cout<<"Stop  at local time ="<<aTrack->GetLocalTime()/s<< "second"<<std::endl;
	   		 return;
		}   			 					    			 
     	} 
  }
  
 //Stop if the particle reach thee maximum allowed number of planet turn 
  if (geometry_type=="SPHERICAL") {
  	//G4cout<<aTrackInfo->GetNbOfPlanetTurn()<<std::endl;
 	if (std::abs(aTrackInfo->GetNbOfPlanetTurn())>= MaxNumberOfTurnAroundThePlanet){
		//G4cout<<"STOP6"<<std::endl;
        	aTrack->SetTrackStatus(fStopAndKill);
 	}
  }	
 
  
 // stop the particles with energies below the particle cutoff energy 
 //------------------------------------------------------------------
  
  //G4cout<<"Stepping3"<<std::endl;
  
  
  
  if (aParticleName.contains("[")) aParticleName="GenericIon";
  if ( UntrackedParticleDic->find(aParticleName) != UntrackedParticleDic->end()){
     	G4double emin=(*UntrackedParticleDic)[aParticleName];
      	G4double ekin=aTrack->GetKineticEnergy();
      	if (ekin <= emin) {
	      // pass the track to the sensitive detector
       		G4SDManager* SDman   = G4SDManager::GetSDMpointer();
       		PLANETOCOSSD* sd = dynamic_cast<PLANETOCOSSD*>
                               		(SDman->FindSensitiveDetector("atmoSD") );
     		       
       		//G4cout<<"Test step1"<<std::endl;
		if ( !aParticleName.contains("nu_"))
	                     	sd->RegisterEkinOfKilledParticle( aTrack);
		//G4cout<<"Test step2"<<std::endl;  
       		aTrack->SetTrackStatus(fStopAndKill);
	
	}
	
   }

   
 // stop particle blocked at a boundary (infinitive loop)
 //---------------------------------------------------------------------------------------------
 
  if (trackID == last_trackID) {
    	if ( (nStep == last_stepNb+1) && (aStep->GetStepLength()<.001*mm || 
	      (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary && 
	           aStep->GetStepLength()<10.*mm  ))) { 
      		nbad++;
       		if (nbad >20){ 
         		aTrack->SetTrackStatus(fStopAndKill);
			if (energy > 0.1*MeV) {
				G4cout<<"Particle boundary problem"<<std::endl;
				G4cout<<aParticleName<<"of "<<energy/MeV<<std::endl;
				PLANETOCOSAnalysisManager::GetInstance()->SwitchOffRegisterResults();
			}	
	  		nbad=0; 
	 	}
      	}
     	else nbad=0;
  }
  else {
    	nbad =0;
  }
  
  last_trackID= trackID;
  last_stepNb= nStep; 
}
///////////////////////////////////////////////////////////////////////////////////////////
//	
void PLANETOCOSSteppingAction::
           AddUntrackedParticle(G4String aParticleName, G4double aEkin)
{ if (UntrackedParticleDic->find(aParticleName) != UntrackedParticleDic->end()){
  	G4cout <<" the  particle is already considered as untracked particle change the energy"  <<G4endl;
  	UntrackedParticleDic->erase(UntrackedParticleDic->find(aParticleName));
  	(*UntrackedParticleDic)[aParticleName]=aEkin; 
  }  
  else {
  	G4ParticleDefinition* aParticleDefinition =
              		G4ParticleTable::GetParticleTable()->FindParticle(aParticleName);
   	if (aParticleDefinition)
     		(*UntrackedParticleDic)[aParticleName]=aEkin;
   	else 
     		G4cout <<" the  particle does not exist in the particle table "  <<G4endl;	    
  }	       
}
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSteppingAction::RemoveUntrackedParticle(G4String aParticleName)
{  if (UntrackedParticleDic->find(aParticleName) != UntrackedParticleDic->end())
   	UntrackedParticleDic->erase(UntrackedParticleDic->find(aParticleName));
   else
     	G4cout <<" the  particle was not in the list of untracked particles "  <<G4endl; 
}
///////////////////////////////////////////////////////////////////////////////////////////
//
void PLANETOCOSSteppingAction::ListUntrackedParticles()
{   
}
