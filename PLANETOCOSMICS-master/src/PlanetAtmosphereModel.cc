#include "PlanetAtmosphereModel.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4UnitsTable.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"

#include "G4PropagatorInField.hh"

#include "myfunctions.hh"
//#include "IAtmosphere.hh"
#include "PlanetUnits.hh"
#include "PlanetManager.hh"
#include "G4RunManager.hh"
#include "PLANETOCOSGeometryConstruction.hh"



PlanetAtmosphereModel::PlanetAtmosphereModel(G4String )
{
  ListOfAtmosphericModels.clear();
  ListOfAtmosphericModels.push_back("TABLE");
  WithMaterial =true;
  
  
   
  
   
}
////////////////////////////////////////////////////////////////////////////////
//
PlanetAtmosphereModel::~PlanetAtmosphereModel()
{ ;
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetAtmosphereModel::ComputeAtmosphericLayers(G4double alt_min,G4double alt_max,
                                   G4double perc_depth, G4double min_tickness,
				   G4double max_tickness,
				   const std::vector<double> flux_detector_altitudes, 
                                   const std::vector<double> flux_detector_depth,
				   std::vector< double >& Altitudes,
				   std::vector< double >& Depth,
				   std::vector< int >& IndexDetectionBoundary,
				   std::vector< G4Material* >& Materials)
{ if (AtmosphericModel != "TABLE")
           ComputeAtmosphereTableFromModel(alt_min,alt_max);
  
  //Problem in initialisation of MSISE90 model to be checked
  if (AtmosphericModel == "MSISE90")
	    ComputeAtmosphereTableFromModel(alt_min,alt_max);
  
  ComputeAtmosphericLayersFromTable(alt_min, alt_max, perc_depth, min_tickness,
				          max_tickness, flux_detector_altitudes, 
                                          flux_detector_depth, Altitudes,Depth,
				          IndexDetectionBoundary,Materials);
 			         
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetAtmosphereModel::ComputeAtmosphericLayers1(G4double alt_min,G4double alt_max,
                                           G4double perc_depth, G4double min_tickness,
				           G4double max_tickness,
				           const std::vector<double> flux_detector_altitudes, 
                                           const std::vector<double> flux_detector_depth,
				           std::vector< double >& Altitudes,
				           std::vector< double >& Depth,
				           std::vector< int >& IndexDetectionLayers,
					   bool detect_at_top,
					   bool detect_at_ground,
				           std::vector< G4Material* >& Materials)
{ if (AtmosphericModel != "TABLE")
           ComputeAtmosphereTableFromModel(alt_min,alt_max);
  
  //Problem in initialisation of MSISE90 model to be checked
  if (AtmosphericModel == "MSISE90")
	    ComputeAtmosphereTableFromModel(alt_min,alt_max);
  
  ComputeAtmosphericLayersFromTable1(alt_min, alt_max, perc_depth, min_tickness,
				          max_tickness, flux_detector_altitudes, 
                                          flux_detector_depth, Altitudes,Depth,
					  IndexDetectionLayers,
					  detect_at_top,
					  detect_at_ground,
				          Materials);
 			         
}						   
////////////////////////////////////////////////////////////////////////////////
//
void PlanetAtmosphereModel::ReadAtmosphereComposition(G4String nameFile)
{  
  
 
  std::fstream File_Input(nameFile,  std::ios::in);
  G4double avogadro;
  avogadro =6.022e+23;
  
  // G4double exp_val =std::exp(1.);
  std::vector< std::vector< double >* > p_data_to_read;
  std::vector< G4String>  atome_molecule_names;
  std::vector< double >   atome_molecule_mass;
  std::vector< std::vector <unsigned int> >   atome_molecule_elemental_composition;
  std::vector<double >    data_units;
  std::vector< std::vector <double>* >  atmosphere_composition;

  //clear vectors
  //------------------
  elements_in_atmosphere.clear();
  atmosphere_mass_density_elements.clear();
  atmosphere_altitude.clear();
  atmosphere_depth.clear();
  atmosphere_mass_density.clear();
  atmosphere_number_density.clear();
  atmosphere_temperature.clear();
  atmosphere_pressure.clear();

  O2_given=false;
  N2_given=false;
 
  //reset the units
  //------------------
  altitude_unit =km;
  mass_density_unit =g/cm3;
  number_density_unit = 1./cm3;
  temperature_unit =kelvin;
  pressure_unit = atmosphere;
  number_of_particle_composition=true;
  
  
  bool SetGroundAltitude =false;
  bool SetTopAltitude =false;
  G4double ground_altitude;
  G4double top_altitude;
  
  //read the first data line
  //if file start by \comments
  //the first line start after \data or \definition
 
  char ch;
  G4String first_word;
  std::stringstream* aline = new std::stringstream(); 
  
  ch='\n';
  // avoid blank lines
  while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
  

  while (ch != '\n'){
 	*aline << ch;
  	File_Input.get(ch);
  }
  *aline>>first_word;
  
  // read the section of comments if exist
  //-----------------------------------------
 
  if (first_word=="\\comments"){
   	while (first_word != "\\data" && 
	       first_word != "\\definition" 
	       && !File_Input.eof()) {
     		
		delete aline;
      		aline = new std::stringstream();
      		File_Input.get(ch);
      		*aline<<ch;
      		while (ch != '\n' && !File_Input.eof()){
        		File_Input.get(ch);
	 		*aline<<ch;
        	}
      		*aline>>first_word;
     	}
      
     
   	if (File_Input.eof()) { 
     		std::cout<<"No atmospheric data have been read"<<std::endl;
       		return;
     	}
  }
 
 //read the section of definition if exist
 //----------------------------------------
  
 if (first_word=="\\definition"){
  	while (first_word != "\\data" && !File_Input.eof()){
     		delete aline;
      		aline = new std::stringstream();
      		File_Input.get(ch);
      		while (ch =='\n' && !File_Input.eof()) 
			File_Input.get(ch);
      		*aline<<ch;
      		while (ch != '\n' && !File_Input.eof()){
        		File_Input.get(ch);
	 		*aline<<ch;
        	}
      		
		*aline>>first_word;
      		if (first_word.find("\\mass_density_unit") == 0 ||
          	    first_word.find("\\number_density_unit") == 0 ||
	  	    first_word.find("\\altitude_unit") == 0 ||
	      	    first_word.find("\\pressure_unit") == 0 ||
	      	    first_word.find("\\temperature_unit") == 0 ||
	  	    first_word.find("\\type_of_composition") == 0 || 
		    first_word.find("\\ground_altitude") == 0 ||
		    first_word.find("\\top_altitude") == 0) {
	 
	   		str_size last= first_word.find("{");
	    		G4String tag_name = first_word;
	    		tag_name = tag_name.remove(last);
	    		tag_name = tag_name.remove(0,1);
            		G4String arg_val = first_word.remove(0,last+1);
	    		if (arg_val.find("}")  == arg_val.size()-1) { 
	     			arg_val = arg_val.remove(arg_val.size()-1);
	      			/*G4cout<<"you have selected the following "<<tag_name<<" : ";
	      			G4cout<<arg_val<<std::endl;*/
	      			if (tag_name.contains("unit")){
	        			G4double aVal =G4UnitDefinition::GetValueOf(arg_val);
	         			if (aVal != 0) {
                    				//G4cout<<"value of unit: "<<aVal<<std::endl;
	             				if (tag_name =="mass_density_unit") mass_density_unit=aVal;
		     				else if (tag_name =="number_density_unit") number_density_unit=aVal;
		     				else if (tag_name =="altitude_unit") altitude_unit=aVal;
		     				else if (tag_name =="pressure_unit") pressure_unit=aVal;
		    
	            			}
	         			else G4cout<<"this unit is not defined in Geant4"<<std::endl;  
	        		}
	      			else if (tag_name == "type_of_composition"){
	       				if (arg_val == "number_of_particles")
	                          		number_of_particle_composition= true;
					else if (arg_val == "mass_composition")
		 		  		number_of_particle_composition= false;
					else G4cout<<"type_of_composition is not valid"<<std::endl;		   
	       			}
				else if (tag_name == "ground_altitude"){
	       				std::stringstream astream;
					astream<<arg_val;
					astream>>ground_altitude;
					SetGroundAltitude =true;		   
	       			}
				else if (tag_name == "top_altitude"){
					std::stringstream astream;
					astream<<arg_val;
					astream>>top_altitude;
					SetTopAltitude =true;
	       						   
	       			}
					  
             		}
            	}
	    
    		if (File_Input.eof()){ 
     			std::cout<<"No atmospheric data have been read"<<std::endl;
       			return;
     		}  
    	}
  }    
 
  if  (first_word=="\\data") { 
  	
	//this is the first line to be read
    	delete aline;
    	aline = new std::stringstream();
    	ch='\n';
    	// avoid blank lines
    	while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
   
    	while (ch != '\n' && !File_Input.eof()){
		*aline<<ch;
      		File_Input.get(ch);
      	}      
    	*aline>>first_word; 
  }
 

 
 //definition of the data to be read
  
  std::vector< G4String > names_of_data;
  names_of_data.push_back(first_word); 
  G4String word;
  G4String last_word="";
  while (!aline->eof()) {
  	*aline>>word; 
      	if (word != last_word) names_of_data.push_back(word);
      	last_word =word;
  }
        
  for (unsigned int i=0;i<names_of_data.size();i++) {
    	//G4cout<<names_of_data[i]<<std::endl;
	G4String lower_name = names_of_data[i];
	lower_name.toLower();
	if (names_of_data[i] =="O2") O2_given=true;
	if (names_of_data[i] =="N2") N2_given=true; 
	if (lower_name == "altitude" ){
		p_data_to_read.push_back(&atmosphere_altitude);
		data_units.push_back(altitude_unit);
	}
	else if (lower_name == "temperature" ){
		p_data_to_read.push_back(&atmosphere_temperature);
		data_units.push_back(temperature_unit);
	}
			  
	else if (lower_name == "pressure"){
	  	p_data_to_read.push_back(&atmosphere_pressure);
		data_units.push_back(pressure_unit);
	}
			   
	else if (lower_name == "density" ){
	    	p_data_to_read.push_back(&atmosphere_mass_density);
		data_units.push_back(mass_density_unit);
	}
	else if (G4Element::GetElement(names_of_data[i])){
		G4Element* theElement = G4Element::GetElement(names_of_data[i]);
	      	atome_molecule_names.push_back(names_of_data[i]);
	      	atome_molecule_mass.push_back(theElement->GetA());		
	      	p_data_to_read.push_back(new std::vector<double>() );		
	      	atmosphere_composition.push_back(p_data_to_read[p_data_to_read.size()-1]);
	     
	      	unsigned int index=10000; 
	      	for (unsigned int i=0;i<elements_in_atmosphere.size();i++){
	         	if (theElement == elements_in_atmosphere[i]){
		   		index=i;
		    		atome_molecule_elemental_composition
		          		.push_back(std::vector<unsigned int>(1,i));
		    	i=elements_in_atmosphere.size();		     
		   	}
		}
	      	if (index ==  10000){
	         	//new element
		  	elements_in_atmosphere.push_back(theElement);
		  	atome_molecule_elemental_composition.push_back
		               	(std::vector<unsigned int>(1,elements_in_atmosphere.size()-1));
		}
	      	if (number_of_particle_composition)
	            	        data_units.push_back(number_density_unit);
	      	else data_units.push_back(mass_density_unit);		        
	    	
	     	}
	else{
		std::vector <G4Element* > aVector 
	                        =ComputeCompositionFromChemicalFormula(names_of_data[i]);
		if (aVector.size() ==0) {
	      		G4cout<< names_of_data[i] <<" was not recognised"<<std::endl;
	       		G4cout<< " No data have been read"<<std::endl;
	       		return;
	      	}
	    	else {
	    		atome_molecule_names.push_back(names_of_data[i]);
	       		atome_molecule_elemental_composition.push_back(std::vector <unsigned int>());
	       		unsigned int last = atome_molecule_elemental_composition.size()-1;
	       		p_data_to_read.push_back(new std::vector<double>() );
	       		atmosphere_composition.
	               		push_back(p_data_to_read[p_data_to_read.size()-1]);
	       		G4double mass =0.; 
	       		for (unsigned i=0;i < aVector.size();i++) {
				mass+=aVector[i]->GetA();
		  		unsigned int index=elements_in_atmosphere.size();
	          		for (unsigned int j=0;j<elements_in_atmosphere.size();j++){
	            			if (aVector[i] == elements_in_atmosphere[j]){
		      				index=j;
		       				atome_molecule_elemental_composition[last].push_back(j);
		       				j=elements_in_atmosphere.size();		     
		      			}
		    		}
		  		if (index ==  elements_in_atmosphere.size()){
	            			//new element
		     			elements_in_atmosphere.push_back(aVector[i]);
		     			atome_molecule_elemental_composition[last]
		                		 .push_back(elements_in_atmosphere.size()-1);
		    		} 
			}
	       		atome_molecule_mass.push_back(mass);
	       		if (number_of_particle_composition) 
				data_units.push_back(number_density_unit);
	       		else data_units.push_back(mass_density_unit); 
		}
	}
  } 
 
 
  do {
  	delete aline;
      	aline = new std::stringstream();
      	ch='\n';
      	// avoid blank lines
      	while (ch =='\n' && !File_Input.eof()) File_Input.get(ch);
   	while (ch != '\n' && !File_Input.eof()){
      		*aline<<ch;
       		File_Input.get(ch);
      	}      
      	if (!File_Input.eof())
        for (unsigned int i=0;i<names_of_data.size();i++){
        	G4double value;
         	*aline>> value;
        	// G4cout<<value<<" "<<i<<std::endl;
		// G4cout<<data_units[i]<<std::endl;
		// G4cout<<value*data_units[i]*cm3/mg<<" "<<i<<std::endl;
        	// G4cout<<p_data_to_read.size()<<std::endl;
        	// G4cout<<p_data_to_read[i]->size()<<std::endl;
	 	(*p_data_to_read[i]).push_back(value*data_units[i]);
  	}  
 }
 while (!File_Input.eof());

  
 // compute mass density if it was not provided in file 
 //-----------------------------------------------------------
 
  if (atmosphere_mass_density.size() <= 0) { // density has to be computed
 	for (unsigned int i=0; i<atmosphere_altitude.size() ; i++) {
     		G4double density = 0.;
      		for (unsigned int j=0; j<atome_molecule_mass.size() ; j++){
       			if (number_of_particle_composition)
          			density+= (*atmosphere_composition[j])[i]*atome_molecule_mass[j]/avogadro;
			else density+=(*atmosphere_composition[j])[i];
	  
       		}
      		atmosphere_mass_density.push_back(density);
    	} 
  }
   
  // compute number density 
  //-------------------------
  
  std::vector< G4double > factor = 
            	std::vector< G4double >(atome_molecule_mass.size(),1.);
  for (unsigned int j=0; j<atome_molecule_mass.size() ; j++){
 	if (!O2_given && (atome_molecule_names[j] == "Oxygen" || 
					atome_molecule_names[j] == "O"))
                                                              	factor[j]=.5;
      	else if (!N2_given && (atome_molecule_names[j] == "Nitrogen" || 
								atome_molecule_names[j] == "N"))
                                                                factor[j]=.5;									      
  } 	    

  for (unsigned int i=0; i<atmosphere_altitude.size() ; i++){
  	G4double density = 0.;
      	G4double ndensity=0.;
      	for (unsigned int j=0; j<atome_molecule_mass.size() ; j++){
       		if (number_of_particle_composition){
          		density+= (*atmosphere_composition[j])[i]*atome_molecule_mass[j]/avogadro;
	   		ndensity+=(*atmosphere_composition[j])[i]*factor[j];
	  	} 
		else {
	  		density+=(*atmosphere_composition[j])[i];
	   		ndensity+= (*atmosphere_composition[j])[i]*factor[j]
	                                       *avogadro/atome_molecule_mass[j];
		}				       
	}
       	atmosphere_number_density.push_back(ndensity*atmosphere_mass_density[i]/density);
  } 
   
   
   
   
  //compute mass density for the different elements contained in the composition
  //----------------------------
   
  for (unsigned int i=0; i<elements_in_atmosphere.size() ; i++)
      		atmosphere_mass_density_elements.push_back(std::vector<double> (atmosphere_altitude.size(),0.));
      
  if (number_of_particle_composition){
  	for (unsigned int i=0; i<atmosphere_altitude.size() ; i++){
       		G4double total = 0.;
         	for (unsigned int j=0; j<atome_molecule_mass.size() ; j++)
                   		total+=(*atmosphere_composition[j])[i]*atome_molecule_mass[j];
	 	   
         	for (unsigned int j=0; j<atome_molecule_mass.size() ; j++){
	   		for (unsigned int k=0; 
			     k<atome_molecule_elemental_composition[j].size();
			     k++){
	      			unsigned int index = atome_molecule_elemental_composition[j][k];
	        		atmosphere_mass_density_elements[index][i]+=
		   				(*atmosphere_composition[j])[i]*elements_in_atmosphere[index]->GetA()
                                                          *atmosphere_mass_density[i]/total; 
	      		}
	      	}
	      
       }
  }
  else {
 	for (unsigned int i=0; i<atmosphere_altitude.size() ; i++) {
       		G4double total = 0.;
         	for (unsigned int j=0; j<atome_molecule_mass.size() ; j++)
                   				total+=(*atmosphere_composition[j])[i];
		for (unsigned int j=0; j<atome_molecule_mass.size() ; j++){
	   		for (unsigned int k=0; k<atome_molecule_elemental_composition[j].size() ; k++){
	      			unsigned int index = atome_molecule_elemental_composition[j][k];
	        		atmosphere_mass_density_elements[index][i]+=
		   		(*atmosphere_composition[j])[i]*elements_in_atmosphere[index]->GetA()
		                    					*atmosphere_mass_density[i]/total/atome_molecule_mass[j]; 	   
	
	      		}
	      	}
      
   	}
  }      
       		 
   
   
 //compute atmospheric depth
 //-------------------------
 atmosphere_depth.clear();
 atmosphere_depth.insert(atmosphere_depth.end(),atmosphere_mass_density.size(),0.);
  
 for (unsigned int i=1; i<atmosphere_altitude.size() ; i++){ 
     	if (atmosphere_altitude[0]>atmosphere_altitude[1]){
       		G4double dens1=atmosphere_mass_density[i];
        	G4double dens0=atmosphere_mass_density[i-1];
		G4double h=atmosphere_altitude[i-1]-atmosphere_altitude[i];
		G4double depth_layer =0;
		if (std::abs(2.*(dens1-dens0)/(dens0+dens1))  <1e-12 ) depth_layer=(dens1+dens0)*h/2.;  
		else depth_layer=(dens1-dens0)*h/std::log(dens1/dens0);
		atmosphere_depth[i]=atmosphere_depth[i-1]+depth_layer;
		/*G4cout<< atmosphere_depth[i]*cm2/g <<'\t'
	      		<< atmosphere_depth[i]*9.81*m/s/s/(100.*pascal)<<'\t'
	      		<< atmosphere_pressure[i]/(100.*pascal)<<'\t'
	      		<< atmosphere_altitude[i]/km <<std::endl;*/
    	}
      	else{
      		G4int i1 = atmosphere_altitude.size()-i-1;
       		G4int i0 = i1+1;
       		G4double dens1=atmosphere_mass_density[i1];
       		G4double dens0=atmosphere_mass_density[i0];
       		G4double h=atmosphere_altitude[i0]-atmosphere_altitude[i1];
       		G4double depth_layer =0;
       		if (std::abs(2.*(dens1-dens0)/(dens0+dens1))  <1e-12 ) depth_layer=(dens1+dens0)*h/2.;  
       		else depth_layer=(dens1-dens0)*h/std::log(dens1/dens0);
       		atmosphere_depth[i1]=atmosphere_depth[i0]+depth_layer;
       		/*G4cout<< atmosphere_depth[i1]*cm2/g <<std::endl;*/
       
      	}        
  }
     
  //k_boltzman
  G4double k_boltzmann;
  k_boltzmann =  1.381e-23*m*newton/kelvin;
  
  // pressure and temperature
  if (atmosphere_pressure.size() == 0 ||  atmosphere_temperature.size() ==0){
  	if (atmosphere_pressure.size() == 0 ){
       		//if the pressure is not provided
         	//the pressure is given by the atmospheric depth * accelaration at earth
         	//                p=depth* 9.81 m/s2
		G4double gravity_at_surface=
			PlanetManager::GetInstance()->GetGplanet()*9.81*m/s/s;
        	for (unsigned int i=0; i<atmosphere_mass_density.size() ;i++)
               		atmosphere_pressure.push_back(1e-3 *pascal 
				  + atmosphere_depth[i]*gravity_at_surface);
       	       
	}
  	if (atmosphere_temperature.size() == 0 ){
  		//if the temperature is not provided
         	//the temperature is deduced from
         	//                p=nkT
	 	G4double k_boltzmann =  1.381e-23*m*newton/kelvin;
		// G4cout<<"k_boltzmann "<<k_boltzmann<<std::endl;
         	for (unsigned int i=0; i<atmosphere_mass_density.size() ;i++)
           		atmosphere_temperature.push_back(atmosphere_pressure[i]/k_boltzmann/atmosphere_number_density[i]);  
	}
  }
     

 //remove the data that we do not need anymore
  p_data_to_read.clear();
  for (unsigned int i = 0; i<atmosphere_composition.size();i++){ 
  	atmosphere_composition[i]->clear();
	delete atmosphere_composition[i];
  } 
  atmosphere_composition.clear();
  atome_molecule_elemental_composition.clear();
  atome_molecule_mass.clear();
  atome_molecule_names.clear();
 
  WriteAtmosphereComposition("Atmo.txt");
  
  
  
  //Set top and ground altitzude in geometry if needed
  
  if (SetGroundAltitude){
  	theGeometry->SetAtmosphereHmin(ground_altitude*altitude_unit);
  }
  if (SetTopAltitude){
  	theGeometry->SetAtmosphereHmax(top_altitude*altitude_unit);
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetAtmosphereModel::WriteAtmosphereComposition(G4String nameFile)
{ std::fstream File_Output(nameFile,  std::ios::out);
 
  File_Output<<"altitude[km]"<<'\t'
            <<"temperature[k]"<<'\t'
            <<"pressure[hPa]"<<'\t'
	    <<"depth[mg/cm2]"<<'\t'
            <<"density[mg/cm3]"<<'\t'
	    <<"ndensity[#/cm3]";
  	    
	    
  for (unsigned int j=0; j<elements_in_atmosphere.size() ; j++)
	       File_Output<<'\t'<<elements_in_atmosphere[j]->GetName()<<"[%mass]";	    
 
  File_Output<<std::endl;    
     
  for (unsigned int i=0; i<atmosphere_altitude.size() ;i++){
   	File_Output<<atmosphere_altitude[i]/km<<'\t'
                <<atmosphere_temperature[i]/kelvin<<'\t'
		<<atmosphere_pressure[i]/100./pascal<<'\t'
		<<atmosphere_depth[i]*cm2/mg<<'\t'
		<<atmosphere_mass_density[i]*cm3/mg<<'\t'
		<<atmosphere_number_density[i]*cm3<<'\t';
     	for (unsigned int j=0; j<elements_in_atmosphere.size() ; j++)
	     	File_Output<<'\t'<<atmosphere_mass_density_elements[j][i]*100.
			                                           /atmosphere_mass_density[i];
	      
     	File_Output<<std::endl;
  }	  
  
  File_Output.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetAtmosphereModel::SetAtmosphericModel(G4String aName)
{ for (unsigned i=0; i<ListOfAtmosphericModels.size();i++){
 	if (aName == ListOfAtmosphericModels[i]) {
		AtmosphericModel=aName;
		return;
	}
  }

  G4cout<<"The atmospheric model that you have selected is not available"<<std::endl;
 
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector <G4Element* >  PlanetAtmosphereModel::
                                ComputeCompositionFromChemicalFormula(G4String aFormula)
{ G4String lower_formula =aFormula;
  lower_formula.toLower();
  std::vector< G4String > substrings;
  std::vector <G4Element* > aVector;
  
  unsigned int ii=0;
  substrings.push_back(aFormula[ii]);
  
  for (unsigned int i=1;i<aFormula.length();i++){
  	if (lower_formula[i] != aFormula[i]) substrings.push_back("");
	substrings[substrings.size()-1] =substrings[substrings.size()-1] + aFormula[i];
  }
  
  for (unsigned int i=0;i<substrings.size();i++){
  	//G4cout<<substrings[i]<<std::endl;
	G4String lower =substrings[i];
	lower.toLower();
	G4String upper =substrings[i];
	upper.toUpper();
	G4String element_name="";
	G4int nb_of_elements;
	nb_of_elements=1;
	
	unsigned int n=substrings[i].length();
	for (unsigned int j=0;j<substrings[i].length();j++) {
	    	if (substrings[i][j] == lower[j] &&  
		             substrings[i][j] == upper[j]){
	       		n=j;
			j=substrings[i].length();
	   	} 
             	else element_name=element_name+substrings[i][j];
	     	if (n==substrings[i].length()) nb_of_elements=1;
	     	else {
	       		G4String nb_of_elements_string = substrings[i].remove(0,n); 
		 	std::stringstream* convert_stream =new std::stringstream();
		 	*convert_stream<<nb_of_elements_string;
		 	*convert_stream>>nb_of_elements;
		 	delete convert_stream;
		 	convert_stream= new std::stringstream();
		 	*convert_stream<<nb_of_elements;
		 	G4String test_string;
		 	*convert_stream>>test_string;
		 	if (test_string != nb_of_elements_string){
				G4cout<<"Invalid symbol"<<nb_of_elements<<std::endl;
		    	aVector.clear();
		    	substrings.clear();
		    	return aVector;
		   	}
		}   
	}
	      
	// find if the element exist
	const G4ElementTable* theElementTable = G4Element::GetElementTable();
	unsigned int j=0;
	unsigned int index=(*theElementTable).size();
	while (j< (*theElementTable).size()){ 
		if (element_name == (*theElementTable)[j]->GetSymbol()) {
			aVector.insert(aVector.end(),nb_of_elements,(*theElementTable)[j]);
		     	index=j;
		     	j=(*theElementTable).size(); 
		}
		j++; 
	}  
	
	if (index == (*theElementTable).size()){
		G4cout<<"the element "<<element_name<<" was not found in the element table"
	                                            <<std::endl;
		aVector.clear();
		substrings.clear();
		return aVector;				      
	}	
  } 
  substrings.clear();
  return aVector;	
} 
////////////////////////////////////////////////////////////////////////////////
//   
void PlanetAtmosphereModel
          ::ComputeAtmosphericLayersFromTable( G4double alt_min,G4double alt_max,
                                   G4double perc_depth, G4double min_thickness,
				   G4double max_thickness,
				   const std::vector<double> flux_detector_altitudes, 
                                   const std::vector<double> flux_detector_depths,
				   std::vector< double >& Altitudes,
				   std::vector< double >& Depths,
				   std::vector< int >& IndexDetectionBoundaries,
				   std::vector< G4Material* >& Materials)
{ Altitudes.clear();
  Depths.clear();
  IndexDetectionBoundaries.clear();
  Materials.clear();
  if (atmosphere_altitude.size() == 0){
  	G4cout<<"No atmospheric data were found"<<std::endl;
    	return;
  }
  
  if ( alt_max<alt_min){
  	G4cout<<"The following conditions was not fulfilled"<<std::endl;
    	G4cout<<"alt_min < alt_max"<<std::endl;
    	return;
  }
 
   
 
 //variable determining the increasing sense of depth
  unsigned int f1=0;
 
  if (atmosphere_altitude[0] <atmosphere_altitude[atmosphere_altitude.size()-1]){
 	f1=atmosphere_altitude.size()-1;
  }

  
  G4double altitude_top =std::min(atmosphere_altitude[f1],alt_max);
 
 
  //find depth_min corresponding at Altitudes[0] and alt_min
  
  G4double depth_min = InterpolateDepthAtAltitude(altitude_top);
  G4double depth_max = InterpolateDepthAtAltitude(alt_min);
       
       
 //compute Altitudes 		
 // the same percent of total depth is containd in one layer except if the 
 // layer tickness is greater than a ximum value or lower than a minimum value		  
 
  G4double d_depth= perc_depth*(depth_max-depth_min)/100.;
  G4double depth_up=depth_min;
  G4double altitude_up =altitude_top;
  
  while( altitude_up >alt_min+0.0000001*mm){
  	Altitudes.push_back(altitude_up);
    	Depths.push_back(depth_up-depth_min);
    
    	G4double depth=depth_up + d_depth;
    	//G4cout<<"depth"<<depth*m2/g<<std::endl;
    	G4double altitude=InterpolateAltitudeAtDepth(atmosphere_altitude,
                                                     atmosphere_depth,
    				                     atmosphere_mass_density,
						     depth);
    	//check if altitude respect the altitude criterium
      
       	// big layers are divided into sublayers with same thickness
      	if ((altitude_up-altitude) > max_thickness) {
		G4int n_sub_layers = int ((altitude_up-altitude) /max_thickness) +1;
       		G4double dalt=(altitude_up-altitude)/n_sub_layers;
       		for (int j=1; j< n_sub_layers;j++){
			G4double altitude1 =altitude_up - dalt*j;
	  		G4double depth1= InterpolateDepthAtAltitude(altitude1);
	  		Altitudes.push_back(altitude1);
	  		Depths.push_back(depth1-depth_min);	 
	  
	     	}
      	}
     	
	//for small layer, a layer that contains 2,3.,4.5,... times the d_depth is constructed   
     	else if ((altitude_up-altitude) < min_thickness) {
		G4double depth1= InterpolateDepthAtAltitude(altitude_up - min_thickness); 
       		depth=depth + d_depth * ( int ((depth1-depth)/d_depth) +1);
       		altitude=InterpolateAltitudeAtDepth(atmosphere_altitude,
                                                atmosphere_depth,
    				                atmosphere_mass_density,depth);
      	}
     
     	depth_up=depth;
     	altitude_up=altitude;
     	
      
  }
  if (  std::abs(Altitudes.back()-alt_min) > 0.001* mm){
  		Altitudes.push_back(alt_min);    
      		Depths.push_back(depth_max); 
  }
  else Altitudes.back()=alt_min;
        
      
  
 // make aflux detection layer altitude vector decreasing and
 //increasing monotically with altitude and depth respectively  
 
  std::vector< double>   detector_altitudes; 
  
  for (unsigned int i=0; i<flux_detector_altitudes.size();i++){
  	G4double altitude = flux_detector_altitudes[i];
	if (detector_altitudes.size()>0){
		if (altitude > detector_altitudes.front()){
	            	detector_altitudes.insert(detector_altitudes.begin(),
		                              	  altitude);
		}					  
	  	else if (altitude < detector_altitudes.back()){
	        	detector_altitudes.push_back(altitude);
		}	   
	  	else{
			G4bool out;
	      		G4int id =myfunc::locate(detector_altitudes, altitude,out); 
	      		detector_altitudes.insert(detector_altitudes.begin()+id+1,
	                                          altitude);
	      
	     	}
	}
        else {
	  	detector_altitudes.push_back(flux_detector_altitudes[i]);
      	}
  }
  
  for (unsigned int i=0; i<flux_detector_depths.size();i++){
  	G4double altitude;
        altitude=0.;
        if (flux_detector_depths[i]<= depth_max-depth_min){
		altitude=InterpolateAltitudeAtDepth(atmosphere_altitude,
                                                    atmosphere_depth,
    				                    atmosphere_mass_density,
					            flux_detector_depths[i]
						    +depth_min);
	 		
	 	unsigned int ndet;
	 	ndet = detector_altitudes.size();
	 
	 	if (ndet>0){
	    		if (altitude > detector_altitudes.front())
	            		detector_altitudes.insert(detector_altitudes.begin(),
		                              altitude);
	    		else if  (altitude < detector_altitudes.back())
	            		detector_altitudes.push_back(altitude);   
	    		else {
	     			G4bool out;
	      			G4int id =myfunc::locate(detector_altitudes, altitude,out); 
	      			detector_altitudes.insert(detector_altitudes.begin()+id+1,
	                                                  altitude);
	      
	     		}
	   	        
	   	}
          	else
	  		detector_altitudes.push_back(altitude);
	}
  } 
 
 //computation of depth and altitude of detection  boundaries 
 ////////////////////////////////////////////////////////////
  for (unsigned int i=0; i<detector_altitudes.size();i++){
  	G4double depth;
       	depth=0.;
	if ( detector_altitudes[i] < Altitudes[0]){
		depth = InterpolateDepthAtAltitude(detector_altitudes[i])-depth_min;
	  	unsigned int index = 0;                                                        
	  	while (detector_altitudes[i] <  Altitudes[index] 
	                                 && index < Altitudes.size()) index++;
						   
	  	if (index < Altitudes.size()){
	   		if (Altitudes[index-1] - detector_altitudes[i] <= 0.0001*mm){
	    			//detection layer coincides with upper bondary 
	     			IndexDetectionBoundaries.push_back(index-1);
	    		}
	    		else if (detector_altitudes[i] - Altitudes[index]  <= 0.0001*mm){
	    			//detection layer coincides with lower bondary 
	     			IndexDetectionBoundaries.push_back(index);
	    		} 
	    		else{
				//insert a new layer 
	     			if (index <Altitudes.size()){
					Altitudes.insert(Altitudes.begin()+index,detector_altitudes[i]);
	       				Depths.insert(Depths.begin()+index,depth);
	      			}
	     			else {
					Altitudes.push_back(detector_altitudes[i]);
	       				Depths.push_back(depth);
	      			}  
	     			IndexDetectionBoundaries.push_back(index);
	    		}
		}
		else{
			// below minimum altitude  the bounadry is set at this minimum
	    		//altitude
             		IndexDetectionBoundaries.push_back(Altitudes.size()-1);
	     		G4cout<<"The flux detection at this altitude will be set at the";
	     		G4cout<<" bottom of the atmosphere"<<std::endl;  
		}   
 	}
	else {	
		G4cout<<"The detection altitude "<<detector_altitudes[i]/km;
	  	G4cout<<" km is above the top of the atmosphere at"<<Altitudes[0]/km;
	  	G4cout<<" km "<<std::endl;
	  	G4cout<<"The flux detection at this altitude will be set"<<
	  	G4cout<<" at the top of the atmosphere"<<std::endl; 
	  	IndexDetectionBoundaries.push_back(0);
	} 
  } 
 

  if (WithMaterial) {
  
  	//construct the vector of materials
  	//----------------------------------
  	for (unsigned int i =0;i<Altitudes.size()-1;i++){
   		G4double density=(Depths[i+1]-Depths[i])/(Altitudes[i]-Altitudes[i+1]);
    		G4double temperature= myfunc::MeanValueOfY(atmosphere_altitude,
                                       		   atmosphere_temperature,
				       		   Altitudes[i],Altitudes[i+1]);
				       
    		// G4cout<<"mean temperature "<<temperature/kelvin<<std::endl;				       
    		G4double pressure= myfunc::MeanValueOfY(atmosphere_altitude,
                                       		atmosphere_pressure,
				       		Altitudes[i],Altitudes[i+1]);
   		// G4cout<<"mean temperature "<<pressure<<std::endl; 				       
    		G4String str_n0;
    		std::stringstream astream;
    		astream<<i+1;
    		astream>>str_n0;
    		G4String name_material ="AtmoLayer"+str_n0;
    
    		std::vector<G4double> mass_ratio;
    		std::vector<int> id;
    		mass_ratio.clear();
    		id.clear();
    		for (unsigned int j=0;j<elements_in_atmosphere.size(); j++){ 
        		G4double ratio=myfunc::IntegrationOfY_exp(atmosphere_altitude,
                                   			  atmosphere_mass_density_elements[j],
				   			  Altitudes[i],
							  Altitudes[i+1])/(Depths[i+1]-Depths[i]);
			// G4cout<<ratio<<" test mass ratio"<<std::endl;
	 		if (ratio>1.e-9 && ratio<=1.) { //avoid NAN and 0 density element 
	   			mass_ratio.push_back(ratio);
	    			id.push_back(int (j));
	   		}	 
		} 
    
    
    		Materials.push_back(new G4Material(name_material,
                            density,
			    int (id.size()), //number of valid elements
		            kStateGas,
		            temperature,
		            pressure));
   		//G4cout<<Altitudes[i]/m<<'\t'<<Altitudes[i+1]/m<<std::endl;
   		//G4cout<<temperature<<std::endl;
   		// G4cout<<pressure/pascal<<std::endl;	
   		// G4cout<<density*cm3/mg<<std::endl;	
    		       
   	 	G4double sum_ratio =0.;


    
    		for (unsigned int j=0;j<id.size(); j++){ 
         		//G4cout<<mass_ratio[j]<<'\t'; 			   
	 		Materials[i]->AddElement(elements_in_atmosphere[id[j]],mass_ratio[j]); 
	 		sum_ratio+=mass_ratio[j];
		} 
		// G4cout<<sum_ratio<<'\t'<<std::endl;
    		mass_ratio.clear();
    		id.clear();
	}	
   	   
        					
 }
}
////////////////////////////////////////////////////////////////////////////////
//   
void PlanetAtmosphereModel
          ::ComputeAtmosphericLayersFromTable1( G4double alt_min,G4double alt_max,
                                   G4double perc_depth, G4double min_thickness,
				   G4double max_thickness,
				   const std::vector<double> flux_detector_altitudes, 
                                   const std::vector<double> flux_detector_depths,
				   std::vector< double >& Altitudes,
				   std::vector< double >& Depths,
				   std::vector< int >& IndexDetectionLayers,
				   bool detect_at_top,
				   bool detect_at_ground,
				   std::vector< G4Material* >& Materials)
{
 
  std::vector< int > IndexDetectionBoundaries;
  WithMaterial =false;
  ComputeAtmosphericLayersFromTable(  alt_min,alt_max,perc_depth, min_thickness,
                                   max_thickness, flux_detector_altitudes, flux_detector_depths, Altitudes,
				   Depths, IndexDetectionBoundaries, Materials);
  detect_at_top=false;
  detect_at_ground=false;
  WithMaterial =true; 
  Materials.clear();
  IndexDetectionLayers.clear(); 
  
  
  //variable determining the increasing sense of depth
  unsigned int f1=0;
 
  if (atmosphere_altitude[0] <atmosphere_altitude[atmosphere_altitude.size()-1]){
 	f1=atmosphere_altitude.size()-1;
  }

  
  G4double altitude_top =std::min(atmosphere_altitude[f1],alt_max);
 
 
  //find depth_min corresponding at Altitudes[0] and alt_min
  
  G4double depth_min = InterpolateDepthAtAltitude(altitude_top);
  G4double depth_max;
  depth_max = InterpolateDepthAtAltitude(alt_min);
  
  // Go for the last time trough the vector todefine detection layer insteead of boundary layer
  // at detection layer
  	G4int incr = 0;
  for (unsigned int i=0; i<IndexDetectionBoundaries.size();i++){
  	unsigned int  index= IndexDetectionBoundaries[i]+incr;
  	if (index!= 0 && index !=Altitudes.size()-1){
		G4double alt0,alt1,alt2; 
		G4double depth0,depth1,depth2;
		G4double h1,h2,dd1,dd2;
		alt0=Altitudes[index-1];
		alt1=Altitudes[index];
		alt2=Altitudes[index+1];
		depth0=Depths[index-1];
		depth1=Depths[index];
		depth2=Depths[index+1];
		dd1=depth1-depth0;
		dd2=depth2-depth1;
		h1=alt0-alt1;
		h2=alt1-alt2;
		G4double dd=std::min(1.*g/cm2,dd1/2.);
		dd=std::min(dd,dd2/2.);
		IndexDetectionLayers.push_back(index);
		Depths[index]=depth1+dd;
		Depths.insert(Depths.begin()+index,depth1-dd);
		Altitudes.insert(Altitudes.begin()+index,InterpolateAltitudeAtDepth(atmosphere_altitude,
                                                    atmosphere_depth,
    				                    atmosphere_mass_density,
					            Depths[index]
						    +depth_min)
						    );
		Altitudes[index+1]=InterpolateAltitudeAtDepth(atmosphere_altitude,
                                                    atmosphere_depth,
    				                    atmosphere_mass_density,
					            Depths[index+1]
						    +depth_min);
		incr=incr+1;
	}
	else if  (index == 0)  detect_at_top=true;
	else detect_at_ground=true;
  }

//construct the vector of materials
//----------------------------------
  for (unsigned int i =0;i<Altitudes.size()-1;i++){
   	G4double density=(Depths[i+1]-Depths[i])/(Altitudes[i]-Altitudes[i+1]);
    	G4double temperature= myfunc::MeanValueOfY(atmosphere_altitude,
                                       		   atmosphere_temperature,
				       		   Altitudes[i],Altitudes[i+1]);
				       
    	// G4cout<<"mean temperature "<<temperature/kelvin<<std::endl;				       
    	G4double pressure= myfunc::MeanValueOfY(atmosphere_altitude,
                                       		atmosphere_pressure,
				       		Altitudes[i],Altitudes[i+1]);
   	// G4cout<<"mean temperature "<<pressure<<std::endl; 				       
    	G4String str_n0;
    	std::stringstream astream;
    	astream<<i+1;
    	astream>>str_n0;
    	G4String name_material ="AtmoLayer"+str_n0;
    
    	std::vector<G4double> mass_ratio;
    	std::vector<int> id;
    	mass_ratio.clear();
    	id.clear();
    	for (unsigned int j=0;j<elements_in_atmosphere.size(); j++){ 
        	G4double ratio=myfunc::IntegrationOfY_exp(atmosphere_altitude,
                                   			  atmosphere_mass_density_elements[j],
				   			  Altitudes[i],
							  Altitudes[i+1])/(Depths[i+1]-Depths[i]);
		// G4cout<<ratio<<" test mass ratio"<<std::endl;
	 	if (ratio>1.e-9 && ratio<=1.) { //avoid NAN and 0 density element 
	   		mass_ratio.push_back(ratio);
	    		id.push_back(int (j));
	   	} 
	} 
    
    
    	Materials.push_back(new G4Material(name_material,
                            density,
			    int (id.size()), //number of valid elements
		            kStateGas,
		            temperature,
		            pressure));
   	//G4cout<<Altitudes[i]/m<<'\t'<<Altitudes[i+1]/m<<std::endl;
   	//G4cout<<temperature<<std::endl;
   	// G4cout<<pressure/pascal<<std::endl;	
   	// G4cout<<density*cm3/mg<<std::endl;	
    		       
   	 G4double sum_ratio =0.;


    
    	for (unsigned int j=0;j<id.size(); j++){ 
         	//G4cout<<mass_ratio[j]<<'\t'; 			   
	 	Materials[i]->AddElement(elements_in_atmosphere[id[j]],mass_ratio[j]); 
	 	sum_ratio+=mass_ratio[j];
	} 
	G4cout<<name_material<<sum_ratio<<'\t'<<std::endl;
    	mass_ratio.clear();
    	id.clear();
   	   
        					
 }
}
////////////////////////////////////////////////////////////////////////////////
//
G4double PlanetAtmosphereModel::InterpolateAltitudeAtDepth(const std::vector< double> altitudes,
                                        const std::vector< double> depths,
					const std::vector< double> densities,
					G4double aDepth)
{ G4double depth1,depth2,alt1,alt2,density1,density2;
  alt1=altitudes[0];
  alt2=altitudes[1];
  depth1=depths[0];
  depth2=depths[1];
  density1=densities[0];
  density2=densities[1];
 
  for (unsigned i=0;i<altitudes.size()-1;i++){
  	alt1=altitudes[i];
   	alt2=altitudes[i+1];
   	depth1=depths[i];
   	depth2=depths[i+1];
   	density1=densities[i];
   	density2=densities[i+1];
   	if ((depth1-aDepth)/(depth1-depth2) <1) i=altitudes.size();
  }
  
  // compute the depth 
  if ( density2 != density1){
  	G4double dh=(depth2-depth1)/(density2-density1);
  	G4double dh1=(alt2-alt1)/std::log(density2/density1);
  	G4double density= density1 + (aDepth- depth1)/dh;
  	G4double altitude= std::log(density/density1)*dh1  + alt1;
  	//G4cout<<"alt "<<altitude/m<<std::endl;
	return altitude;
 }
 else { 
 	G4double dh=(alt2-alt1)/(depth2-depth1);
 	G4double altitude = alt1 + dh * (aDepth -depth1);
	//G4cout<<"alt1 "<<alt1/m<<std::endl;
	//G4cout<<"alt2 "<<alt2/m<<std::endl;
	//G4cout<<"depth2 "<<depth2*cm2/g<<std::endl;
	//G4cout<<"depth1 "<<depth1*cm2/g<<std::endl;
	//G4cout<<"aDepth "<<aDepth*cm2/g<<std::endl;
	//G4cout<<"alt "<<altitude/m<<std::endl;
	return altitude;
 }	
}
////////////////////////////////////////////////////////////////////////////////
//
G4double PlanetAtmosphereModel::InterpolateDepthAtAltitude(G4double anAltitude)
{ G4double depth1,depth2,alt1,alt2,density1,density2,exp_val;
  exp_val=1.;
  exp_val=std::pow(2.71828,0.);
 
  alt1=atmosphere_altitude[0];
  alt2=atmosphere_altitude[1];
  depth1=atmosphere_depth[0];
  depth2=atmosphere_depth[1];
  density1=atmosphere_mass_density[0];
  density2=atmosphere_mass_density[1];
 
  for (unsigned i=0;i<atmosphere_altitude.size()-1;i++){
  	alt1=atmosphere_altitude[i];
   	alt2=atmosphere_altitude[i+1];
   	depth1=atmosphere_depth[i];
   	depth2=atmosphere_depth[i+1];
   	density1=atmosphere_mass_density[i];
   	density2=atmosphere_mass_density[i+1];
   	if ((alt1-anAltitude)/(alt1-alt2) <1) i=atmosphere_altitude.size();
   
  }
  
 // compute the depth 
  if (density1 != density2){
    	G4double density =
    	    density1*std::pow(density2/density1,(anAltitude-alt1)/(alt2-alt1));
   	G4double dh=(depth2-depth1)/(density2-density1);
   	G4double depth=depth1+dh*(density-density1);
	return depth;
  }
  else{ 
 	G4double dh=(depth2-depth1)/(alt2-alt1);
       	G4double depth=depth1+dh*(anAltitude-alt1); 
       	return depth;
  }   
}
