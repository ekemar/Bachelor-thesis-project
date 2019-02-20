
#include "PlanetSoil.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4UImanager.hh"
#include"G4RunManager.hh"
#include"PlanetSoilMessenger.hh"
#include"PlanetManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
PlanetSoil::PlanetSoil()
{ SoilMaterials.clear();
  SoilThickness.clear();
  SoilNbSubLayers.clear();
  nel_of_new_layer =0;
  nel_already_defined =0;
  theElementTable = const_cast<G4ElementTable* > (G4Element::GetElementTable());
  ElementsOfNewLayer.clear();
  ConcentrationOfElementsInNewLayer.clear();
  theMessenger= new PlanetSoilMessenger(this);
}
////////////////////////////////////////////////////////////////////////
//
PlanetSoil::~PlanetSoil()
{; 
}
////////////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::AddMonoElementLayerAndSetThickness(G4String el_name_or_symbol, G4double density, G4double thickness,G4int nb_sub_layers)
{ if (nel_of_new_layer != nel_already_defined){ 
	G4cout<<"The definition of the composition of the last "
              <<"layer is not finished"<<std::endl;
  	return;
  }
  G4Element* anElement = G4Element::GetElement( el_name_or_symbol);
  if (!anElement){ //check if symbol existG4cout<<"PlanetSoil 351"<<std::endl;
  	for (unsigned int i=0;i<theElementTable->size();i++){
	  	if (el_name_or_symbol == (*theElementTable)[i]->GetSymbol()) {
			anElement = (*theElementTable)[i];
			i=theElementTable->size();
		}
		
	}  
  	
  }
  std::stringstream astream;
  G4String str_i;
  astream<<SoilMaterials.size()+1;
  astream>>str_i;
  G4String material_name="Soil"+str_i;
  if (anElement) {
  	G4Material* aMaterial = new G4Material(material_name,density,1);
	aMaterial->AddElement(anElement,1.);
	SoilMaterials.push_back(aMaterial);
	SoilThickness.push_back(thickness);	
  }
  else {// try chemical formula
        std::vector<G4Element*> aVectorOfElement;
  	std::vector<G4int> occurences;
        ComputeCompositionFromChemicalFormula(el_name_or_symbol,
					      aVectorOfElement,
					      occurences);
  	unsigned int nb_elements = aVectorOfElement.size();
	if (nb_elements >0 ){
		G4double total_mass =0.;
		for (unsigned i=0;i<nb_elements;i++){
			G4double mass =double(occurences[i])*
						aVectorOfElement[i]->GetA();
			total_mass +=mass;
		}
		G4Material* aMaterial = new G4Material(material_name,density,nb_elements);
		for (unsigned i=0;i<nb_elements;i++){
			G4double mass =double(occurences[i])*
						aVectorOfElement[i]->GetA();
			aMaterial->AddElement(aVectorOfElement[i],mass/total_mass); 
		}
		SoilMaterials.push_back(aMaterial);
	        SoilThickness.push_back(thickness);
		SoilNbSubLayers.push_back(nb_sub_layers);
		
		
	}
	else {
		G4cout<<"Impossible to build the material made of "
		      <<el_name_or_symbol<<std::endl;
	}
	
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::AddMonoElementLayerAndSetDepth(G4String el_name_or_symbol, G4double density, G4double depth,G4int nb_sub_layers)
{AddMonoElementLayerAndSetThickness(el_name_or_symbol, density, depth/density,nb_sub_layers); 
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::AddLayerAndSetThickness(G4int nel, G4double density, G4double thickness,G4int nb_sub_layers)
{ if (nel_of_new_layer != nel_already_defined){ 
	G4cout<<"The definition of the composition of the last "
              <<"layer is not finished"<<std::endl;
  	return;
  } 	      
		     
  nel_of_new_layer =nel;
  nel_already_defined =0;
  ElementsOfNewLayer.clear();
  ConcentrationOfElementsInNewLayer.clear();
  density_of_newlayer = density;
  thickness_of_newlayer =thickness;
  nb_sub_layers_of_newlayer =nb_sub_layers;
  
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::AddLayerAndSetDepth(G4int nel, G4double density, G4double depth,G4int nb_sub_layers)
{AddLayerAndSetThickness( nel, density, depth/density,nb_sub_layers);
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::AddElementToLayer(G4String el_name_or_symbol, G4double weight_concentration)
{ 
	if (nel_already_defined == nel_of_new_layer ) {
		G4cout<<"All elements of the layer have been already defined!"<<std::endl;
  		return;
  }
 
 
  G4Element* anElement = G4Element::GetElement( el_name_or_symbol);
  if (!anElement){ //check if symbol exist
  	for (unsigned int i=0;i<theElementTable->size();i++){
	  	if (el_name_or_symbol == (*theElementTable)[i]->GetSymbol()) {
			anElement = (*theElementTable)[i];
			i=theElementTable->size();
		}
		
	}  
  	
  }
  if (anElement) {
  	unsigned int index=ElementsOfNewLayer.size();
  	for (unsigned int i=0; i<ElementsOfNewLayer.size();i++){
		if (anElement == ElementsOfNewLayer[i]) {
			ConcentrationOfElementsInNewLayer[i]+=weight_concentration;
			index =i;
			i = ElementsOfNewLayer.size();
		
		}
	}	
	if (index == ElementsOfNewLayer.size()) {
			ElementsOfNewLayer.push_back(anElement);
			ConcentrationOfElementsInNewLayer.push_back(weight_concentration);
	}
	nel_already_defined++;
	
  		
  }
  else {// try chemical formula
        std::vector<G4Element*> aVectorOfElement;
  	std::vector<G4int> occurences;
        ComputeCompositionFromChemicalFormula(el_name_or_symbol,
					      aVectorOfElement,
					      occurences);
  	unsigned int nb_elements = aVectorOfElement.size();
	if (nb_elements >0 ){
		G4double total_mass =0.;
		for (unsigned i=0;i<nb_elements;i++){
			G4double mass =double(occurences[i])*
						aVectorOfElement[i]->GetA();
			total_mass +=mass;
		}
		for (unsigned i=0;i<nb_elements;i++){
			G4double mass =double(occurences[i])*
						aVectorOfElement[i]->GetA();
			G4double concentration =mass*weight_concentration/total_mass;			
			unsigned int index=ElementsOfNewLayer.size();
  			for (unsigned int j=0; j<ElementsOfNewLayer.size();j++){
				if (aVectorOfElement[i] == ElementsOfNewLayer[j]) {
					ConcentrationOfElementsInNewLayer[j]+=concentration;
					index =j;
					j = ElementsOfNewLayer.size();
		
				}
			}	
			if (index == ElementsOfNewLayer.size()) {
				ElementsOfNewLayer.push_back(aVectorOfElement[i]);
				ConcentrationOfElementsInNewLayer.push_back(concentration);
			} 
		}
		nel_already_defined++;
	}
	else {
		G4cout<<"Impossible to add an element  made of "
		      <<el_name_or_symbol<<std::endl;
		return;
	}
	
  }
  if (nel_already_defined == nel_of_new_layer ) {
		std::stringstream astream;
             	G4String str_i;
  		astream<<SoilMaterials.size()+1;
  		astream>>str_i;
  		G4String material_name="Soil"+str_i;
		G4int nb_elements = int(ElementsOfNewLayer.size());
		G4Material* aMaterial = new G4Material(material_name,density_of_newlayer,nb_elements);
		G4double total_concentration=0;
		for (unsigned int i=0;i<ElementsOfNewLayer.size();i++){
			total_concentration+=ConcentrationOfElementsInNewLayer[i]; 
		}
		for (unsigned int i=0;i<ElementsOfNewLayer.size();i++){
			aMaterial->AddElement(ElementsOfNewLayer[i],
					      ConcentrationOfElementsInNewLayer[i]/total_concentration); 
		}
		SoilMaterials.push_back(aMaterial);
	        SoilThickness.push_back(thickness_of_newlayer);
		SoilNbSubLayers.push_back(nb_sub_layers_of_newlayer);
		
		
		
  }
  
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::ComputeCompositionFromChemicalFormula(G4String  aFormula,
						      std::vector <G4Element* >& aVectorOfElement,
					              std::vector <G4int>& occurences)
{ G4String lower_formula =aFormula;
  lower_formula.toLower();
  std::vector< G4String > substrings;
  aVectorOfElement.clear();
  
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
		    	aVectorOfElement.clear();
			occurences.clear();
		    	substrings.clear();
		    	return ;
		   	}
		}   
	}
	      
	// find if the element exist
	const G4ElementTable* theElementTable = G4Element::GetElementTable();
	unsigned int j=0;
	unsigned int index=(*theElementTable).size();
	while (j< (*theElementTable).size()){ 
		if (element_name == (*theElementTable)[j]->GetSymbol()) {
			aVectorOfElement.push_back((*theElementTable)[j]);
			occurences.push_back(nb_of_elements);
		     	index=j;
		     	j=(*theElementTable).size(); 
		}
		j++; 
	}  
	
	if (index == (*theElementTable).size()){
		G4cout<<"the element "<<element_name<<" was not found in the element table"
	                                            <<std::endl;
		aVectorOfElement.clear();
		occurences.clear();
		substrings.clear();
		return ;				      
	}	
  } 
  substrings.clear();
  return;
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::ResetLayers()
{ SoilMaterials.clear();
  SoilThickness.clear();
  SoilNbSubLayers.clear();
  nel_of_new_layer =0;
  nel_already_defined =0;
  ElementsOfNewLayer.clear();
  ConcentrationOfElementsInNewLayer.clear();
}
////////////////////////////////////////////////////////////////////////
//
void PlanetSoil::SetDefault()
{ ResetLayers();
  G4String planet_name = PlanetManager::GetInstance()->GetPlanetName();
  if (planet_name == "Mars"){
  	AddLayerAndSetThickness(9, 1.5 *g/cm3, 10.*m,1);
  	AddElementToLayer("Na2O", 1.5);
	AddElementToLayer("MgO",7.7);
	AddElementToLayer("Al2O3",8.1);
	AddElementToLayer("SiO2",46.8);
	AddElementToLayer("SO3",6.);
	AddElementToLayer("K2O",0.2);
	AddElementToLayer("CaO",6.2);
	AddElementToLayer("TiO2",1.1);
	AddElementToLayer("Fe2O3",18.8);
  }
  else if (planet_name == "Mercury"){
  	AddLayerAndSetThickness(7, 1.3 *g/cm3, 10.*m,1);
	AddElementToLayer("MgO",35.);
	AddElementToLayer("Al2O3",7);
	AddElementToLayer("SiO2",45.);
	AddElementToLayer("Na2O",0.7);
	AddElementToLayer("CaO",7);
	AddElementToLayer("TiO2",0.3);
	AddElementToLayer("FeO",5.);
  	
  }
  else if (planet_name == "Earth"){
  	AddMonoElementLayerAndSetThickness("SiO2",1.7 *g/cm3, 10.*m,1);
  	
  }
  // by default for Jupiter no soil
  
   	
}
