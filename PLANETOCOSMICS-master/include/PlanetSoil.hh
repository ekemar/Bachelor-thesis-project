#ifndef PlanetSoil_HH
#define PlanetSoil_HH 

#include "globals.hh"
#include"G4ios.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"
#include "G4Element.hh"
#include "G4Material.hh"

class PlanetSoilMessenger;
class PlanetSoil 
{
public:
		       
         PlanetSoil();
	 ~PlanetSoil();
	 void AddMonoElementLayerAndSetThickness(G4String el_name_or_symbol, G4double density, G4double thickness, G4int nb_sub_layers);
	 void AddMonoElementLayerAndSetDepth(G4String el_name_or_symbol, G4double density, G4double depth, G4int nb_sub_layers);
	 void AddLayerAndSetThickness(G4int nel, G4double density, G4double thickness, G4int nb_sub_layers);
	 void AddLayerAndSetDepth(G4int nel, G4double density, G4double depth, G4int nb_sub_layers);
	 void AddElementToLayer(G4String el_name_or_symbol, G4double weight_concentration);
         void ComputeCompositionFromChemicalFormula(G4String  aFormula,
	 					    std::vector <G4Element* >&,
						    std::vector <G4int > &);
	 void ResetLayers();
	 
	 void SetDefault();
	 
	 //GetMethods
	 inline std::vector<G4Material*> GetSoilMaterials(){return SoilMaterials;}
	 inline std::vector<G4double> GetSoilThickness(){return SoilThickness;}
	 inline std::vector<G4int> GetSoilNbSubLayers(){return SoilNbSubLayers;}

private:
	PlanetSoilMessenger* theMessenger;
	std::vector<G4Material*> SoilMaterials;
	std::vector<G4double> SoilThickness;
	std::vector<G4Element*> ElementsOfNewLayer;
	std::vector<G4double> ConcentrationOfElementsInNewLayer;
	std::vector<G4int> SoilNbSubLayers;
	G4int nel_of_new_layer;
	G4int nel_already_defined; 
	G4ElementTable* theElementTable;
	G4double density_of_newlayer;
	G4double thickness_of_newlayer;
	G4int nb_sub_layers_of_newlayer; 
	
	

};

#endif
