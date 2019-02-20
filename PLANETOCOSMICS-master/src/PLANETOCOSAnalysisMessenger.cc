#include "PLANETOCOSAnalysisMessenger.hh"
#include "PLANETOCOSAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserStackingAction.hh"
#include "PLANETOCOSStackingAction.hh"
#include "G4RunManager.hh"


PLANETOCOSAnalysisMessenger::PLANETOCOSAnalysisMessenger
                                        (PLANETOCOSAnalysisManager* theManager)
                                                  :pAnalysisManager(theManager)
{
 //parameters
  //-----------
  
  G4UIparameter* histo_label = new G4UIparameter("HistoLabel",'s',false);
  
 
 
  G4UIparameter* particleName = new G4UIparameter("ParticleName",'s',false);
  particleName->SetDefaultValue("all");
  
  G4UIparameter* nbin_alt = new G4UIparameter("nbin_alt",'i',false);
  nbin_alt->SetDefaultValue("200");
  
  G4UIparameter* nbin_depth = new G4UIparameter("nbin_depth",'i',false);
  nbin_depth->SetDefaultValue("200");
  
  G4UIparameter* nbin_energy = new G4UIparameter("nbin_energy",'i',false);
  nbin_energy->SetDefaultValue("20");
  
  G4UIparameter* nbin_latitude = new G4UIparameter("nbin_latitude",'i',false);
  nbin_latitude->SetDefaultValue("20");
  
  G4UIparameter* nbin_longitude = new G4UIparameter("nbin_longitude",'i',false);
  nbin_longitude->SetDefaultValue("20");
  
  G4UIparameter* energy_min = new G4UIparameter("energy_min",'d',false);
  energy_min->SetDefaultValue("0.01");
 
  G4UIparameter* energy_max = new G4UIparameter("energy_max",'d',false);
  energy_max->SetDefaultValue("1000.");
  
  
  G4UIparameter* latitude_min = new G4UIparameter("latitude_min",'d',false);
  latitude_min->SetDefaultValue("-90");
 
  G4UIparameter* latitude_max = new G4UIparameter("latitude_max",'d',false);
  latitude_max->SetDefaultValue("90.");
  
  G4UIparameter* longitude_min = new G4UIparameter("longitude_min",'d',false);
  longitude_min->SetDefaultValue("-180");
 
  G4UIparameter* longitude_max = new G4UIparameter("longitude_max",'d',false);
  longitude_max->SetDefaultValue("180");
  
  
  G4UIparameter* energy_unit = new G4UIparameter("energy_unit",'s',false);
  energy_unit->SetDefaultValue("MeV ");
  energy_unit->SetParameterCandidates("GeV keV MeV eV");
  
  G4UIparameter* nbin_cos = new G4UIparameter("nbin_cos",'i',false);
  nbin_cos->SetDefaultValue("20");
  
  G4UIparameter* cos_min = new G4UIparameter("cos_min",'d',false);
  cos_min->SetDefaultValue("-1.");
  cos_min->SetParameterRange(" cos_min >= -1 && cos_min <= 1");
 
  G4UIparameter* cos_max = new G4UIparameter("cos_max",'d',false);
  cos_max->SetDefaultValue("1.");
  cos_max->SetParameterRange(" cos_max >= -1 && cos_max <=1");
  
  G4UIparameter* nbin_angle = new G4UIparameter("nbin_angle",'i',false);
  nbin_angle->SetDefaultValue("180");
  
  G4UIparameter* nbin_x = new G4UIparameter("nbin_x",'i',false);
  nbin_x->SetParameterRange(" nbin_x  >0 ");
  nbin_x->SetDefaultValue("100");
  
  G4UIparameter* nbin_y = new G4UIparameter("nbin_y",'i',false);
  nbin_y->SetParameterRange(" nbin_y  >0 ");
  nbin_y->SetDefaultValue("100");
  
  
  
  
  G4UIparameter* angle_min = new G4UIparameter("angle_min",'d',false);
  angle_min->SetDefaultValue("0");
 
  G4UIparameter* angle_max = new G4UIparameter("angle_max",'d',false);
  angle_max->SetDefaultValue("360");
  
  
  G4UIparameter* scale_type = new G4UIparameter("scale_type",'s',false);
  scale_type->SetDefaultValue("linear");
  scale_type->SetParameterCandidates("Linear linear LINEAR lin Lin LIN Log log LOG");
  
  G4UIparameter* tree_filename = new G4UIparameter("tree_filename",'s',false);

  
  G4UIparameter* tree_store_type = new G4UIparameter("tree_store_type",'s',false);
  
  G4UIparameter* tree_dir = new G4UIparameter("tree_dir",'s',true);
  tree_dir->SetDefaultValue("/");
  
  G4UIparameter* tree_filename_opt = new G4UIparameter("tree_filename",'s',true);
  tree_filename_opt->SetDefaultValue("nofile"); 
  
  G4UIparameter* tree_store_type_opt = new G4UIparameter("tree_store_type",'s',true);
  tree_store_type_opt->SetDefaultValue("xml");   

  G4UIparameter* tree_dir_opt = new G4UIparameter("tree_dir",'s',true);
  tree_dir_opt->SetDefaultValue("/"); 
  
  G4UIparameter* tree_path = new G4UIparameter("tree_path",'s',true);
  tree_path->SetDefaultValue("/");
  
  G4UIparameter* true_bool = new G4UIparameter("true_bool",'b',true);
  true_bool->SetDefaultValue(true);
  
  G4UIparameter* false_bool = new G4UIparameter("false_bool",'b',true);
  false_bool->SetDefaultValue(false);
  
  G4UIparameter* normalisation_type = 
                    new G4UIparameter("type of normalisation",'s',true);
  normalisation_type->SetParameterCandidates(" NONE TO_PRIMARY_FLUX PER_PRIMARY");
  normalisation_type->SetDefaultValue("NONE");
  
  
  G4UIparameter* nbin_nb_turn = new G4UIparameter("nbin_nb_turn",'i',false);
  G4UIparameter* max_nb_turn =  new G4UIparameter("max_nb_turn",'d',false);
  max_nb_turn->SetParameterRange("max_nb_turn >0.");
  
  G4UIparameter* max_nb_equator_crossing =  new G4UIparameter("max_nb_eq_cross",'i',false);
  max_nb_equator_crossing->SetParameterRange("max_nb_eq_cross > 0"); 
  
  
  G4UIparameter* nbin_time = new G4UIparameter("nbin_time",'i',false);
  nbin_time->SetDefaultValue("20");
  
  G4UIparameter* time_min = new G4UIparameter("time_min",'d',false);
  time_min->SetDefaultValue("0.");
 
  G4UIparameter* time_max = new G4UIparameter("time_max",'d',false);
  time_max->SetDefaultValue("1000000.");
  
  
  G4UIparameter* time_unit = new G4UIparameter("time_unit",'s',false);
  time_unit->SetDefaultValue("s ");
  time_unit->SetParameterCandidates("s hour minute day");
  //command directory 
  //-------------------
  
  analysisDir = new G4UIdirectory("/PLANETOCOS/ANALYSIS/");
  analysisDir->SetGuidance("Analysis control");
  
  primfluxDir= new G4UIdirectory("/PLANETOCOS/ANALYSIS/PRIMARY/");
  primfluxDir->SetGuidance("Control of primary flux histograms");
   
  detfluxDir= new G4UIdirectory("/PLANETOCOS/ANALYSIS/SECONDARY/");
  detfluxDir->SetGuidance("Control of secondary flux histograms");
  
  completeEventDir= new G4UIdirectory("/PLANETOCOS/ANALYSIS/COMPLETE_EVENT/");
  completeEventDir->SetGuidance("Control of complete Event root-File.");
    
  edepDir= new G4UIdirectory("/PLANETOCOS/ANALYSIS/EDEP/");
  edepDir->SetGuidance("Control of energy deposited histograms");
  
  pseudotrappingDir = new G4UIdirectory("/PLANETOCOS/ANALYSIS/QUASITRAPPED/");
  pseudotrappingDir->SetGuidance("Control of  histograms for registering quasi trapped particle information");
  
  //commands
  //--------
  G4String title;
  SelectDetectorCmd =  new G4UIcmdWithAnInteger("/PLANETOCOS/ANALYSIS/SECONDARY/SelectDetector",this);
  SelectDetectorCmd->SetGuidance("Select a detector for creating flux histogram");
  SelectDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DeselectDetectorCmd =  new G4UIcmdWithAnInteger("/PLANETOCOS/ANALYSIS/SECONDARY/UnselectDetector",this);
  DeselectDetectorCmd->SetGuidance("Unselect a detector");
  DeselectDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SelectAllDetectorsCmd =  new G4UIcmdWithoutParameter("/PLANETOCOS/ANALYSIS/SECONDARY/SelectAllDetectors",this);
  SelectAllDetectorsCmd->SetGuidance("Select all detectors for creating flux histogram");
  SelectAllDetectorsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DeselectAllDetectorsCmd =  new G4UIcmdWithoutParameter("/PLANETOCOS/ANALYSIS/SECONDARY/UnselectAllDetectors",this);
  DeselectAllDetectorsCmd->SetGuidance("Unselect all detectors");
  DeselectAllDetectorsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
 
  
  //Edep in atmosphere histo creation
  //-------------------
  
  
  
  CreateEdepVsDepthHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/EDEP/AtmoEdepVsDepthHisto",this);
  CreateEdepVsDepthHistoCmd
         ->SetGuidance("Create an histogram to register the deposited energy in the atmosphere  in function of depth [g/cm2]");
  CreateEdepVsDepthHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateEdepVsDepthHistoCmd->SetParameter(histo_label);
  CreateEdepVsDepthHistoCmd->SetParameter(nbin_depth);
  
  CreateEdepVsAltitudeHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/EDEP/AtmoEdepVsAltitudeHisto",this);
  CreateEdepVsAltitudeHistoCmd
         ->SetGuidance("Create an histogram to register the deposited energy in the atmosphere in function of altitude [km]");
  CreateEdepVsAltitudeHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateEdepVsAltitudeHistoCmd->SetParameter(histo_label);
  CreateEdepVsAltitudeHistoCmd->SetParameter(nbin_alt);
  
   //Edep in soil histo creation
  //-------------------
  
  
  
  CreateSoilEdepVsDepthHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/EDEP/SoilEdepVsDepthInGPerCm2Histo",this);
  CreateSoilEdepVsDepthHistoCmd
         ->SetGuidance("Create an histogram to register the deposited energy in the soil in function of depth [g/cm2]");
  CreateSoilEdepVsDepthHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateSoilEdepVsDepthHistoCmd->SetParameter(histo_label);
  CreateSoilEdepVsDepthHistoCmd->SetParameter(nbin_depth);
  
  CreateSoilEdepVsThicknessHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/EDEP/SoilEdepVsDepthInKmHisto",this);
  CreateSoilEdepVsThicknessHistoCmd
         ->SetGuidance("Create an histogram to register the deposited energy in the soil in function of depth [km] ");
  CreateSoilEdepVsThicknessHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateSoilEdepVsThicknessHistoCmd->SetParameter(histo_label);
  CreateSoilEdepVsThicknessHistoCmd->SetParameter(nbin_alt);
  
  //root file creation for complete event
  //------------------------------------- 
  
  CreateFluxRootFileCmd = new G4UIcommand("/PLANETOCOS/ANALYSIS/COMPLETE_EVENT/CreateFluxRootFile",this);
  CreateFluxRootFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  G4UIparameter*param = new G4UIparameter("filename",'s',false);
  CreateFluxRootFileCmd->SetParameter(param);
  
  
  WriteFluxRootFileCmd = new G4UIcmdWithoutParameter("/PLANETOCOS/ANALYSIS/COMPLETE_EVENT/WriteFluxRoot",this);
  WriteFluxRootFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  //flux histo creation
  //-------------------
  
  
  G4String guidance1;
  CreateDownwardFluxVsEkinHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/DownwardFluxHisto",this);
  guidance1 ="Create an Histogram to register the downward "; 
  guidance1 += " flux of a particle at the selected detectors";
  CreateDownwardFluxVsEkinHistoCmd->SetGuidance(guidance1);
  CreateDownwardFluxVsEkinHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(particleName);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(histo_label);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(nbin_energy);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(energy_min);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(energy_max);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(energy_unit);
  CreateDownwardFluxVsEkinHistoCmd->SetParameter(scale_type);
  
  
  
  CreateUpwardFluxVsEkinHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/UpwardFluxHisto",this);
  guidance1 ="Create an Histogram to register the Upward omnidirectional"; 
  guidance1 += " flux of a particle at the selected detectors";
  CreateUpwardFluxVsEkinHistoCmd->SetGuidance(guidance1);
  CreateUpwardFluxVsEkinHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(particleName);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(histo_label);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(nbin_energy);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(energy_min);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(energy_max);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(energy_unit);
  CreateUpwardFluxVsEkinHistoCmd->SetParameter(scale_type);
 
  CreateEnergyVsCosZenithHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/CosZenithVsEkinHisto",this);
  guidance1 ="Create a 2DHistogram to register the cos_zenith vs Ekin distribution "; 
  guidance1 += " of the flux of  particle at the selected detectors";
  CreateEnergyVsCosZenithHistoCmd->SetGuidance(guidance1);
  CreateEnergyVsCosZenithHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(particleName);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(histo_label);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(nbin_energy);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(energy_min);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(energy_max);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(energy_unit);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(scale_type);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(nbin_cos);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(cos_min);
  CreateEnergyVsCosZenithHistoCmd->SetParameter(cos_max);
  
  CreateCosZenithHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/CosZenithHisto",this);
  guidance1 ="Create an histogram to register the angular zenith distribution"; 
  guidance1 += " of the flux of a particle at the selected detectors";
  CreateCosZenithHistoCmd->SetGuidance(guidance1);
  CreateCosZenithHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateCosZenithHistoCmd->SetParameter(particleName);
  CreateCosZenithHistoCmd->SetParameter(histo_label);
  CreateCosZenithHistoCmd->SetParameter(nbin_cos);
  CreateCosZenithHistoCmd->SetParameter(cos_min);
  CreateCosZenithHistoCmd->SetParameter(cos_max);
  
 
  CreateAzimuthHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/AzimuthHisto",this);
  guidance1 ="Create an histogram to register the angular azimuth distribution"; 
  guidance1 += " of the flux of a particle at the selected detectors";
  CreateAzimuthHistoCmd->SetGuidance(guidance1);
  CreateAzimuthHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateAzimuthHistoCmd->SetParameter(particleName);
  CreateAzimuthHistoCmd->SetParameter(histo_label);
  CreateAzimuthHistoCmd->SetParameter(nbin_angle);
  CreateAzimuthHistoCmd->SetParameter(angle_min);
  CreateAzimuthHistoCmd->SetParameter(angle_max);
 
  
  CreateDownwardFluxVsPosHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/DownwardIntegralFluxVsPosHisto",this);
  guidance1 ="Create an Histogram to register the downward omnidirectional"; 
  guidance1 += " flux of a particle at the selected detectors in function of position";
  CreateDownwardFluxVsPosHistoCmd->SetGuidance(guidance1);
  CreateDownwardFluxVsPosHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(particleName);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(histo_label);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(nbin_longitude);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(longitude_min);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(longitude_max);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(nbin_latitude);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(latitude_min);
  CreateDownwardFluxVsPosHistoCmd->SetParameter(latitude_max);
 
  
  CreateUpwardFluxVsPosHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/UpwardIntegralFluxVsPosHisto",this);
  guidance1 ="Create an Histogram to register the Upward omnidirectional"; 
  guidance1 += " flux of a particle at the selected detectors in function of position";
  CreateUpwardFluxVsPosHistoCmd->SetGuidance(guidance1);
  CreateUpwardFluxVsPosHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(particleName);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(histo_label);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(nbin_longitude);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(longitude_min);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(longitude_max);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(nbin_latitude);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(latitude_min);
  CreateUpwardFluxVsPosHistoCmd->SetParameter(latitude_max);
  
  
  
  CreateDownwardFluxVsPosFlatHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/DownwardIntegralFluxVsFlatPosHisto",this);
  guidance1 ="Create an Histogram to register the downward omnidirectional"; 
  guidance1 += " flux of a particle at the selected detectors in function of xy position";
  CreateDownwardFluxVsPosFlatHistoCmd->SetGuidance(guidance1);
  CreateDownwardFluxVsPosFlatHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateDownwardFluxVsPosFlatHistoCmd->SetParameter(particleName);
  CreateDownwardFluxVsPosFlatHistoCmd->SetParameter(histo_label);
  CreateDownwardFluxVsPosFlatHistoCmd->SetParameter(nbin_x);
  CreateDownwardFluxVsPosFlatHistoCmd->SetParameter(nbin_y);
 
  
  CreateUpwardFluxVsPosFlatHistoCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/UpwardIntegralFluxVsFlatPosHisto",this);
  guidance1 ="Create an Histogram to register the Upward omnidirectional"; 
  guidance1 += " flux of a particle at the selected detectors in function of xy position";
  CreateUpwardFluxVsPosFlatHistoCmd->SetGuidance(guidance1);
  CreateUpwardFluxVsPosFlatHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateUpwardFluxVsPosFlatHistoCmd->SetParameter(particleName);
  CreateUpwardFluxVsPosFlatHistoCmd->SetParameter(histo_label);
  CreateUpwardFluxVsPosFlatHistoCmd->SetParameter(nbin_x);
  CreateUpwardFluxVsPosFlatHistoCmd->SetParameter(nbin_y);
  
  
  SetEnergyLimitsCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/SetIntegralFluxEkinRange",this);
  guidance1 ="Define the  energy range over which integral flux  are integrated  for flux vs position histos"; 
  SetEnergyLimitsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetEnergyLimitsCmd->SetGuidance(guidance1);
  SetEnergyLimitsCmd->SetParameter(energy_min);
  SetEnergyLimitsCmd->SetParameter(energy_max);
  SetEnergyLimitsCmd->SetParameter(energy_unit);
  
  
  SetLatLongLimitsCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/SECONDARY/SetLatitudeLongitudeLimits",this);
  guidance1 ="Define the longitude and latitude limits for the next created secondary histo "; 
  SetLatLongLimitsCmd->SetGuidance(guidance1);
  SetLatLongLimitsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetLatLongLimitsCmd->SetParameter(latitude_min);
  SetLatLongLimitsCmd->SetParameter(latitude_max);
  SetLatLongLimitsCmd->SetParameter(longitude_min);
  SetLatLongLimitsCmd->SetParameter(longitude_max);
  
  
  SetDivideByCosThCmd = new
             G4UIcmdWithAString("/PLANETOCOS/ANALYSIS/SECONDARY/SetTypeOfWeight",this);
  guidance1 ="If the parameter is INVERSE_COSTH (ONE) the weight of registered particles is porportiobnla to 1/cos_th (1) in the next created secondary histograms "; 
  SetDivideByCosThCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetDivideByCosThCmd->SetGuidance(guidance1);
  
  //histo primary
  //-------------
  
  CreateEnergyVsZenithPrimaryCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/PRIMARY/CosZenVsEnergyHisto",this);
  guidance1 ="Create an histogram  to register the flux of primary particle ";
  guidance1 += "  cos_zenith vs energy";
  CreateEnergyVsZenithPrimaryCmd->SetGuidance(guidance1);
  CreateEnergyVsZenithPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(particleName);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(histo_label);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(nbin_energy);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(energy_min);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(energy_max);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(energy_unit);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(scale_type);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(nbin_cos);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(cos_min);
  CreateEnergyVsZenithPrimaryCmd->SetParameter(cos_max);
  
  CreateZenithVsAzimuthPrimaryCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/PRIMARY/CosZenVsAzimuthHisto",this);
  guidance1 ="Create an histogram to register the integral flux of primaries";
  guidance1 += "  cos_zenith vs azimuth";
  CreateZenithVsAzimuthPrimaryCmd->SetGuidance(guidance1);
  CreateZenithVsAzimuthPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(particleName);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(histo_label);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(nbin_angle);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(angle_min);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(angle_max);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(nbin_cos);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(cos_min);
  CreateZenithVsAzimuthPrimaryCmd->SetParameter(cos_max);
  
  CreateOmniFluxPrimaryCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/PRIMARY/FluxHisto",this);
  guidance1 ="Create an histogram to register the omnidirectional"; 
  guidance1 += " flux of a primaries";
  CreateOmniFluxPrimaryCmd->SetGuidance(guidance1);
  CreateOmniFluxPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateOmniFluxPrimaryCmd->SetParameter(particleName);
  CreateOmniFluxPrimaryCmd->SetParameter(histo_label);
  CreateOmniFluxPrimaryCmd->SetParameter(nbin_energy);
  CreateOmniFluxPrimaryCmd->SetParameter(energy_min);
  CreateOmniFluxPrimaryCmd->SetParameter(energy_max);
  CreateOmniFluxPrimaryCmd->SetParameter(energy_unit);
  CreateOmniFluxPrimaryCmd->SetParameter(scale_type);
  
  CreateZenithPrimaryCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/PRIMARY/CosZenithHisto",this);
  guidance1 ="Create an histogram to register the angular zenith distribution"; 
  guidance1 += " of primaries";
  CreateZenithPrimaryCmd->SetGuidance(guidance1);
  CreateZenithPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateZenithPrimaryCmd->SetParameter(particleName);
  CreateZenithPrimaryCmd->SetParameter(histo_label);
  CreateZenithPrimaryCmd->SetParameter(nbin_cos);
  CreateZenithPrimaryCmd->SetParameter(cos_min);
  CreateZenithPrimaryCmd->SetParameter(cos_max);
 
  CreateAzimuthPrimaryCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/PRIMARY/AzimuthHisto",this);
  guidance1 ="Create an histogram to register the angular azimuth distribution"; 
  guidance1 += " of primaries";
  CreateAzimuthPrimaryCmd->SetGuidance(guidance1);
  CreateAzimuthPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateAzimuthPrimaryCmd->SetParameter(particleName);
  CreateAzimuthPrimaryCmd->SetParameter(histo_label);
  CreateAzimuthPrimaryCmd->SetParameter(nbin_angle);
  CreateAzimuthPrimaryCmd->SetParameter(angle_min);
  CreateAzimuthPrimaryCmd->SetParameter(angle_max);
  
  //histogram tree commands
  //-----------------------
  
  ResetHistoCmd = new
      G4UIcmdWithoutParameter("/PLANETOCOS/ANALYSIS/ResetHistograms", this);
  ResetHistoCmd->SetGuidance("Reset all the histograms");
  ResetHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ResetTreeCmd = new
               G4UIcmdWithoutParameter("/PLANETOCOS/ANALYSIS/ResetTree", this);
  ResetTreeCmd->SetGuidance("Remove all the histograms from the tree");
  ResetTreeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   
  SaveTreeCmd = new G4UIcommand("/PLANETOCOS/ANALYSIS/SaveTree",this);
  SaveTreeCmd->SetGuidance("Save the tree in your selected file");
  guidance1 = "[usage] /PLANETOCOS/ANALYSIS/SaveTree file_name store_type";
  guidance1 += " tree_dir normalisation_type";
  SaveTreeCmd->SetGuidance(guidance1); 
  SaveTreeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SaveTreeCmd->SetParameter(tree_filename);
  SaveTreeCmd->SetParameter(tree_store_type);
  SaveTreeCmd->SetParameter(tree_dir);
  SaveTreeCmd->SetParameter(normalisation_type);
  
  
  PrintTreeCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/WriteTreeInASCIIFile",this);
  guidance1= "Write histograms from a selected subtree  of a selected ";
  guidance1+="histogram tree in your selected ASCII file";
  PrintTreeCmd->SetGuidance(guidance1);
  guidance1="[usage] /PLANETOCOS/ANALYSIS/WriteTreeInASCIIFile file_name";
  guidance1+=" subtree_path [normalisation_type tree_file_name  tree_file_type]";
  PrintTreeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  PrintTreeCmd->SetParameter(tree_filename);
  PrintTreeCmd->SetParameter(tree_path);
  PrintTreeCmd->SetParameter(normalisation_type);
  PrintTreeCmd->SetParameter(tree_filename_opt);
  PrintTreeCmd->SetParameter(tree_store_type_opt);

#ifndef USE_ANALYSIS_ROOT 
  AddTreeFilesCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/AddTreeFiles",this);
  guidance1="Add histograms from two diffrent tree files";
  AddTreeFilesCmd->SetGuidance(guidance1);
  guidance1="[usage] /PLANETOCOS/ANALYSIS/AddTreeFiles file1 store_type1 dir1";
  guidance1= guidance1 +" file2 store_type2 dir2";
  AddTreeFilesCmd->SetGuidance(guidance1);
  AddTreeFilesCmd->SetParameter(tree_filename);
  AddTreeFilesCmd->SetParameter(tree_store_type);
  AddTreeFilesCmd->SetParameter(tree_dir);
  AddTreeFilesCmd->SetParameter(tree_filename);
  AddTreeFilesCmd->SetParameter(tree_store_type);
  AddTreeFilesCmd->SetParameter(tree_dir);
  
  AddFileToTreeCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/AddFileToTree",this);
  guidance1="Add histograms of a tree file to the tree";
  AddFileToTreeCmd->SetGuidance(guidance1);
  guidance1="[usage] /PLANETOCOS/ANALYSIS/AddFileToTree file store_type dir";
  AddFileToTreeCmd->SetGuidance(guidance1);
  AddFileToTreeCmd->SetParameter(tree_filename);
  AddFileToTreeCmd->SetParameter(tree_store_type);
  AddFileToTreeCmd->SetParameter(tree_dir);
  
  AddTreeToFileCmd = new
             G4UIcommand("/PLANETOCOS/ANALYSIS/AddTreeToFile",this);
  guidance1="Add histograms of the tree to the histograms of a file tree";
  AddTreeToFileCmd->SetGuidance(guidance1);
  guidance1="[usage] /PLANETOCOS/ANALYSIS/AddTreeToFile file store_type dir";
  AddTreeToFileCmd->SetGuidance(guidance1);
  AddTreeToFileCmd->SetParameter(tree_filename);
  AddTreeToFileCmd->SetParameter(tree_store_type);
  AddTreeToFileCmd->SetParameter(tree_dir);
#endif  
  
  //post track hit histo creation
  //-------------------

  guidance1 = "Create an histogram that registers the number of planet";
  guidance1 += "turns of particle during the simulation vs start ekin";
  
  CreateNbOfPlanetTurnVsStartEkinHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/NbOfPlanetTurnVsStartEkinHisto",this);
  
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetGuidance(guidance1);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(particleName);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(histo_label);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(nbin_energy);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(energy_min);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(energy_max);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(energy_unit);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(scale_type);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(nbin_nb_turn);
  CreateNbOfPlanetTurnVsStartEkinHistoCmd->SetParameter(max_nb_turn);
  
  guidance1 = "Create an histogram that registers the number of planet";
  guidance1 += "turns of particle during the simulation vs life time";
  
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/NbOfPlanetTurnVsLifeTimeHisto",this);
  
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetGuidance(guidance1);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(particleName);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(histo_label);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(nbin_time);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(time_min);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(time_max);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(time_unit);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(scale_type);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(nbin_nb_turn);
  CreateNbOfPlanetTurnVsLifeTimeHistoCmd->SetParameter(max_nb_turn);
  
  guidance1 = "Create an histogram that registers the number of equator crossing";
  guidance1 += "of particle during the simulation vs start ekin";
  
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/NbOfEquatorCrossingVsStartEkinHisto",this);
  
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetGuidance(guidance1);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(particleName);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(histo_label);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(nbin_energy);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(energy_min);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(energy_max);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(energy_unit);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(scale_type);
  CreateNbOfEquatorCrossingVsStartEkinHistoCmd->SetParameter(max_nb_equator_crossing);
  
  guidance1 = "Create an histogram that registers the number of equator crossing";
  guidance1 += "of particle during the simulation vs life time";
  
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/NbOfEquatorCrossingVsLifeTimeHisto",this);
  
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetGuidance(guidance1);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(particleName);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(histo_label);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(nbin_time);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(time_min);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(time_max);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(time_unit);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(scale_type);
  CreateNbOfEquatorCrossingVsLifeTimeHistoCmd->SetParameter(max_nb_equator_crossing);
  
  guidance1 = "Create an histogram that registers the number of equator crossing";
  guidance1 += "of particle during the simulation vs nb of planet turn";
  
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/NbOfEquatorCrossingVsNbOfPlanetTurnHisto",this);
  
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetGuidance(guidance1);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetParameter(particleName);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetParameter(histo_label);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetParameter(nbin_nb_turn);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetParameter(max_nb_turn);
  CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd->SetParameter(max_nb_equator_crossing);
  
  guidance1 = "Create an histogram that registers the life time of a particle";
  guidance1 += "turning around the planer in function of start kinertic energy";
  
  CreateLifeTimeVsStartEkinHistoCmd
  	= new G4UIcommand("/PLANETOCOS/ANALYSIS/QUASITRAPPED/LifeTimeVsStartEkinHisto",this);
  
  CreateLifeTimeVsStartEkinHistoCmd->SetGuidance(guidance1);
  CreateLifeTimeVsStartEkinHistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(particleName);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(histo_label);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(nbin_energy);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(energy_min);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(energy_max);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(energy_unit);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(scale_type);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(nbin_time);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(time_min);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(time_max);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(time_unit);
  CreateLifeTimeVsStartEkinHistoCmd->SetParameter(scale_type);
  
  
  
       
  //Security save cmd
  //-----------------
  
  SetSecuritySaveCmd= new
             G4UIcmdWithABool("/PLANETOCOS/ANALYSIS/SetSecuritySave",this);
  SetSecuritySaveCmd->
    SetGuidance("If true the tree is saved in your selected file after your selected number of events ");
  SetSecuritySaveCmd->SetParameterName("SecuritySave",true);
  SetSecuritySaveCmd->SetDefaultValue(true);
  SetSecuritySaveCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetSecurityFileNameCmd = 
            new G4UIcommand("/PLANETOCOS/ANALYSIS/SetSecurityFile",this);
  guidance1 = "Set the type and name of the file where the simulation results";
  guidance1 += " should be saved after a given  number";
  guidance1 += " of events fort security reason.";
  SetSecurityFileNameCmd->SetGuidance(guidance1);
  SetSecurityFileNameCmd->SetParameter(tree_filename);
  SetSecurityFileNameCmd->SetParameter(tree_store_type);
  
  
  SetSecurityNbOfEventsCmd =  new G4UIcmdWithAnInteger("/PLANETOCOS/ANALYSIS/SetSecurityNbOfEvents",this);
  SetSecurityNbOfEventsCmd->SetGuidance("Set number of events after which a security save should be performed");
  SetSecurityNbOfEventsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  //Normalisation flux command
  
 /* SetPrimaryIncidentFluxCmd = new G4UIcmdWithADoubleAndUnit("/PLANETOCOS/ANALYSIS/SetIncidentFlux",this);
  SetPrimaryIncidentFluxCmd->SetGuidance("Set the flux of incident primary for normalisation purpose"); 
  SetPrimaryIncidentFluxCmd->SetParameterName("Incident flux",false);
  SetPrimaryIncidentFluxCmd->SetUnitCategory("Integral omni flux");*/
  

  
  
  //bug verbosity
  SetBugVerboseCmd = new
  	 G4UIcmdWithAnInteger("/PLANETOCOS/ANALYSIS/SetBugVerbosity",this);
  G4String guidance = "If n>0 the user is informed when an event is aborted due to a bug";	 
  SetBugVerboseCmd->SetGuidance(guidance);
  SetBugVerboseCmd->SetParameterName("verbosity",false);
  SetBugVerboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
   
  
  
  
  
 
 
  
}

PLANETOCOSAnalysisMessenger::~PLANETOCOSAnalysisMessenger()
{delete analysisDir;
}

void PLANETOCOSAnalysisMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  //Select detector
  if (command == SelectDetectorCmd)  
     pAnalysisManager->GetFluxDetectionAnalyser()->SelectDetector(SelectDetectorCmd->GetNewIntValue(newValues));			      			      
  
  else if (command == DeselectDetectorCmd)  
     pAnalysisManager->GetFluxDetectionAnalyser()->DeselectDetector(DeselectDetectorCmd->GetNewIntValue(newValues));			      			      
  
  else if (command == SelectAllDetectorsCmd)	      			      
                     pAnalysisManager->GetFluxDetectionAnalyser()->SelectAllDetectors();
  
  else if (command == DeselectAllDetectorsCmd)	      			      
                     pAnalysisManager->GetFluxDetectionAnalyser()->DeselectAllDetectors();		     
  
  //edep histo			      
  else if (command == CreateEdepVsDepthHistoCmd){
  	const char* paramString=newValues;
     	G4String histo_label;
     	G4int nbin_depth;
     	std::istringstream is((char*)paramString);
     	is >> histo_label>> nbin_depth ;
     	pAnalysisManager->GetEdepAnalyser()->CreateEdepVsDepth(histo_label,nbin_depth);
                                           
  }
  else if (command == CreateEdepVsAltitudeHistoCmd){
  	const char* paramString=newValues;
     	G4String histo_label;
     	G4int nbin_alt;
     	std::istringstream is((char*)paramString);
     	is >> histo_label>> nbin_alt ;
     	pAnalysisManager->GetEdepAnalyser()->CreateEdepVsAltitude(histo_label,nbin_alt);
                                           
  } 
  
  else if (command == CreateSoilEdepVsDepthHistoCmd){
  	const char* paramString=newValues;
     	G4String histo_label;
     	G4int nbins;
     	std::istringstream is((char*)paramString);
     	is >> histo_label>> nbins ;
     	pAnalysisManager->GetSoilEdepAnalyser()->CreateEdepVsDepth(histo_label,nbins);
                                           
  }
  else if (command == CreateSoilEdepVsThicknessHistoCmd){
  	const char* paramString=newValues;
     	G4String histo_label;
     	G4int nbins;
     	std::istringstream is((char*)paramString);
     	is >> histo_label>> nbins ;
     	pAnalysisManager->GetSoilEdepAnalyser()->CreateEdepVsThickness(histo_label,nbins);
                                           
  }
  
  //secondary histo cmd	
  else if (command == CreateDownwardFluxVsEkinHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type;
     	G4double energy_min,energy_max;
     	G4int nbin_energy;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;

     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateDownwardFluxVsEkinHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,scale_type);
   }
 
   else if (command == CreateUpwardFluxVsEkinHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type;
     	G4double energy_min,energy_max;
     	G4int nbin_energy;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
     	
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateUpwardFluxVsEkinHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type); 
   }					  
   else if (command == CreateEnergyVsCosZenithHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type;
     	G4double energy_min,energy_max;
     	G4int nbin_energy;
	G4double cos_min,cos_max;
     	G4int nbin_cos;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
     	is >> nbin_cos >> cos_min >>cos_max;
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateEnergyVsCosZenithHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type,
					  nbin_cos,cos_min,cos_max); 
   }   
   
  
   else if (command == CreateCosZenithHistoCmd){
   	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double cos_min,cos_max;
     	G4int nbin_cos;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_cos >> cos_min >>cos_max;
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateCosZenithHisto(particle_name, histo_label,
       		                                   nbin_cos,cos_min,cos_max); 
   }
   else if (command == CreateAzimuthHistoCmd){
   	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double angle_min,angle_max;
     	G4int nbin_angle;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_angle >> angle_min >>angle_max;
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateAzimuthHisto(particle_name, histo_label,
       	                                   nbin_angle,angle_min,angle_max); 
   }
   
   else if (command == CreateDownwardFluxVsPosHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type;
     	G4double lon_min,lon_max,lat_min,lat_max;
     	G4int nbin_lon,nbin_lat;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
	is >> nbin_lon>> lon_min>>lon_max
	   >> nbin_lat>> lat_min>>lat_max;
     	
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateDownwardFluxVsPosHisto(particle_name, histo_label,
                                                       nbin_lon,lon_min,lon_max,
					               nbin_lat,lat_min,lat_max); 
   }
   else if (command == CreateUpwardFluxVsPosHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type;
     	G4double lon_min,lon_max,lat_min,lat_max;
     	G4int nbin_lon,nbin_lat;
     	
	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
	is >> nbin_lon>> lon_min>>lon_max
	   >> nbin_lat>> lat_min>>lat_max;
     	
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateUpwardFluxVsPosHisto(particle_name, histo_label,
                                                       nbin_lon,lon_min,lon_max,
					               nbin_lat,lat_min,lat_max); 
   }
   
   else if (command == CreateDownwardFluxVsPosFlatHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4int nbin_x,nbin_y;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
	is >> nbin_x>>  nbin_y;
     
     	pAnalysisManager->GetFluxDetectionAnalyser()->CreateDownwardFluxVsPosFlatHisto(particle_name, histo_label,
                                                       nbin_x, nbin_y); 
   }
   else if (command == CreateUpwardFluxVsPosFlatHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4int nbin_x,nbin_y;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
	is >> nbin_x>>  nbin_y;
	
	pAnalysisManager->GetFluxDetectionAnalyser()->CreateUpwardFluxVsPosFlatHisto(particle_name, histo_label,
                                                       nbin_x, nbin_y); 
   }
   else if (command == SetEnergyLimitsCmd ){
       
        const char* paramString=newValues;
     	G4String energy_unit;
     	G4double energy_min,energy_max;
     	std::istringstream is((char*)paramString);
     	is >> energy_min >>energy_max >> energy_unit;
     	
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->SetEkinMin(energy_min);
	pAnalysisManager->GetFluxDetectionAnalyser()->SetEkinMax(energy_max);
   
   }
   else if (command == SetLatLongLimitsCmd ){
       
        const char* paramString=newValues;
     	G4double lon_min,lon_max,lat_min,lat_max;
     	
	std::istringstream is((char*)paramString);
     
	is>>lat_min>>lat_max>> lon_min>>lon_max;
	
     	pAnalysisManager->GetFluxDetectionAnalyser()->SetMinLat(lat_min);
	pAnalysisManager->GetFluxDetectionAnalyser()->SetMaxLat(lat_max);
	pAnalysisManager->GetFluxDetectionAnalyser()->SetMinLong(lon_min);
	pAnalysisManager->GetFluxDetectionAnalyser()->SetMaxLong(lon_max);
   
   }
   else if (command == SetDivideByCosThCmd ){
        G4bool aBool=false;
	if (newValues =="ONE") aBool = false;
	else if  (newValues =="INVERSE_COSTH") aBool = true;
	else {
		G4cout<<"Your input parameter is not a good candidate"<<std::endl; 
   		return;
	}
	pAnalysisManager->GetFluxDetectionAnalyser()->SetDivideByCosTh(aBool);		
   }	
   //primary histo cmd
   else if (command == CreateEnergyVsZenithPrimaryCmd){
   	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type;
     	G4double energy_min,energy_max, cos_min,cos_max;
     	G4int nbin_energy, nbin_cos;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
     	//G4cout<<scale_type<<std::endl;
     	is >> nbin_cos >> cos_min >>cos_max;
	
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
     	pAnalysisManager->GetPrimaryFluxAnalyser()->CreateEnergyVsZenithPrimaryHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type,
					  nbin_cos,cos_min,cos_max);	 
   }	
   else if (command == CreateZenithVsAzimuthPrimaryCmd){
   	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double azim_min,azim_max, cos_min,cos_max;
     	G4int nbin_azim, nbin_cos;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_azim >> azim_min >>azim_max;
     	is >> nbin_cos >> cos_min >>cos_max;
     	pAnalysisManager->GetPrimaryFluxAnalyser()->CreateZenithVsAzimuthPrimaryHisto(particle_name, histo_label,
                                          nbin_cos,cos_min,cos_max,
					  nbin_azim,azim_min,azim_max);	 
   }
   else if (command == CreateOmniFluxPrimaryCmd){
   	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type;
     	G4double energy_min,energy_max;
     	G4int nbin_energy;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
     	//G4cout<<scale_type<<std::endl;
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
     	pAnalysisManager->GetPrimaryFluxAnalyser()->OmnifluxFluxPrimaryHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type); 
   }
   else if (command == CreateZenithPrimaryCmd){
    	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double cos_min,cos_max;
     	G4int nbin_cos;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_cos >> cos_min >>cos_max ;
     	pAnalysisManager->GetPrimaryFluxAnalyser()->CreateZenithPrimaryHisto(particle_name, histo_label,
                                          nbin_cos,cos_min,cos_max); 
   }
   
   else if (command == CreateAzimuthPrimaryCmd) {
   	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double angle_min,angle_max;
     	G4int nbin_angle;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_angle >> angle_min >>angle_max ;
     	pAnalysisManager->GetPrimaryFluxAnalyser()->CreateAzimuthPrimaryHisto(particle_name, histo_label,
                                          nbin_angle,angle_min,angle_max); 
   }
   
// PVD -----------------------------------------------------
   
   else if (command == CreateFluxRootFileCmd) {
   	const char* paramString=newValues;
     	G4String filename;
     	std::istringstream is((char*)paramString);
     	is >> filename;  

	pAnalysisManager->GetFluxDetectionAnalyser()->SetDetectFlux();	
	pAnalysisManager->SetBaseFilenameRootFile(filename);
	pAnalysisManager->CreateFluxRootFile();
   }
   
   else if (command == WriteFluxRootFileCmd) {    			      
     	pAnalysisManager->WriteFluxRoot();
   }
   
// PVD -----------------------------------------------------   
  
   
   //pseudo trapping
   
   else if (command == CreateNbOfPlanetTurnVsLifeTimeHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type,  time_unit;
     	G4double time_min, time_max,turn_max;
     	G4int nbin_time, nbin_turn;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_time >> time_min >>time_max >> time_unit>>scale_type;
        is>> nbin_turn >>turn_max;	
     	
	time_min *= G4UnitDefinition::GetValueOf(time_unit);
     	time_max *= G4UnitDefinition::GetValueOf(time_unit);
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateNbOfPlanetTurnVsLifeTimeHisto(particle_name, histo_label,
                                          nbin_time,time_min,time_max,
					  scale_type,
					  nbin_turn,turn_max); 
   }
   else if (command == CreateNbOfPlanetTurnVsStartEkinHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type,  energy_unit;
     	G4double energy_min, energy_max,turn_max;
     	G4int nbin_energy, nbin_turn;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
        is>> nbin_turn >>turn_max;	
     	
	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateNbOfPlanetTurnVsStartEkinHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type,
					  nbin_turn,turn_max); 
   } 
   else if (command == CreateNbOfEquatorCrossingVsLifeTimeHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type,  time_unit;
     	G4double time_min, time_max;
     	G4int nbin_time, max_nb_crossing;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_time >> time_min >>time_max >> time_unit>>scale_type;
        is>> max_nb_crossing;	
     	
	time_min *= G4UnitDefinition::GetValueOf(time_unit);
     	time_max *= G4UnitDefinition::GetValueOf(time_unit);
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateNbOfEquatorCrossingVsLifeTimeHisto(particle_name, histo_label,
                                          nbin_time,time_min,time_max,
					  scale_type,
					  max_nb_crossing); 
   }
   else if (command == CreateNbOfEquatorCrossingVsStartEkinHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label, scale_type,  energy_unit;
     	G4double energy_min, energy_max;
     	G4int nbin_energy,max_nb_crossing;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type;
        is>> max_nb_crossing;	
     	
	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateNbOfEquatorCrossingVsStartEkinHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type,
					  max_nb_crossing); 
   } 
   else if (command == CreateNbOfEquatorCrossingVsNbOfPlanetTurnHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, histo_label;
     	G4double turn_max;
     	G4int  max_nb_crossing, nbin_turn;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_turn >>turn_max;
        is>> max_nb_crossing;	
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateNbOfEquatorCrossingVsNbOfPlanetTurnHisto(particle_name, histo_label,
                                               nbin_turn, turn_max,
					       max_nb_crossing); 
   }
   
   else if (command == CreateLifeTimeVsStartEkinHistoCmd){
  	const char* paramString=newValues;
     	G4String particle_name, energy_unit, histo_label, scale_type1, 
	         scale_type2, time_unit;
     	G4double energy_min,energy_max,time_min, time_max;
     	G4int nbin_energy, nbin_time;
     	std::istringstream is((char*)paramString);
     	is >> particle_name >> histo_label;
     	is >> nbin_energy >> energy_min >>energy_max >> energy_unit>>scale_type1
	   >> nbin_time >> time_min >>time_max >> time_unit>>scale_type2;
     	
     	energy_min *= G4UnitDefinition::GetValueOf(energy_unit);
     	energy_max *= G4UnitDefinition::GetValueOf(energy_unit);
	time_min *= G4UnitDefinition::GetValueOf(time_unit);
     	time_max *= G4UnitDefinition::GetValueOf(time_unit);
     	
	
	pAnalysisManager->GetPseudoTrappingAnalyser()->CreateLifeTimeVsStartEkinHisto(particle_name, histo_label,
                                          nbin_energy,energy_min,energy_max,
					  scale_type1,
					  nbin_time,time_min,time_max,
					  scale_type2 ); 
   } 
   
   //general cmds      
   else if (command == ResetHistoCmd) pAnalysisManager->ResetHistograms(); 
   
   else if (command == ResetTreeCmd) pAnalysisManager->ResetTree();  
   
   else if (command == SaveTreeCmd){
   	const char* paramString=newValues;
     	G4String file_name, store_type, tree_dir,normalisation_type;
     	std::istringstream is((char*)paramString);
     	is >> file_name >> store_type >>tree_dir>>normalisation_type;
     	if (normalisation_type == "NONE") {
     		pAnalysisManager->SetTypeOfNormalisation(normalisation_type);
		pAnalysisManager->SaveTree(file_name,store_type,tree_dir,false);
	}
	else {  pAnalysisManager->SetTypeOfNormalisation(normalisation_type);
		pAnalysisManager->SaveTree(file_name,store_type,tree_dir,true);
	}
    }
  
   
   else if (command == PrintTreeCmd){
   	const char* paramString=newValues;
     	G4String file_name, path, tree_file, store_type, normalisation_type;
     	std::istringstream is((char*)paramString);
     	is >> file_name >> path >> normalisation_type>>tree_file >> store_type;
     	if (tree_file == "nofile"){
	        pAnalysisManager->SetTypeOfNormalisation(normalisation_type);
		pAnalysisManager->PrintHistogramTree(file_name,path);
	}	
     	else 
       	pAnalysisManager->PrintHistogramTree(file_name, path, 
                                            tree_file, store_type);	   
   }
#ifndef USE_ANALYSIS_ROOT   
 
   else if (command == AddTreeFilesCmd){
   	const char* paramString=newValues;
     	G4String file_name1, store_type1, tree_dir1;
     	G4String file_name2, store_type2, tree_dir2;
     	std::istringstream is((char*)paramString);
     	is >> file_name1 >> store_type1 >>tree_dir1;
     	is >> file_name2 >> store_type2 >>tree_dir2;
     	pAnalysisManager->AddTreeFiles(file_name1,store_type1,tree_dir1,
                                    		file_name2,store_type2,tree_dir2);
   }   
   else if (command == AddFileToTreeCmd){
   	const char* paramString=newValues;
     	G4String file_name, store_type, tree_dir;
     	std::istringstream is((char*)paramString);
     	is >> file_name >> store_type >>tree_dir;
     	pAnalysisManager->AddFileToTree(file_name,store_type,tree_dir);
   } 
   
   else if (command == AddTreeToFileCmd){
   	const char* paramString=newValues;
     	G4String file_name, store_type, tree_dir;
     	std::istringstream is((char*)paramString);
     	is >> file_name >> store_type >>tree_dir;
     	pAnalysisManager->AddTreeToFile(file_name,store_type,tree_dir);
   }    
#endif      
   else if (command ==   SetSecuritySaveCmd){
   	pAnalysisManager
           ->SetSecuritySave(SetSecuritySaveCmd->GetNewBoolValue(newValues));
   }
   else if  (command ==   SetSecurityFileNameCmd){
   	const char* paramString=newValues;
     	G4String file_name, store_type;
     	std::istringstream is((char*)paramString);
     	is>>file_name>>store_type;
     	pAnalysisManager->SetSecuritySaveFilename(file_name);
     	pAnalysisManager->SetSecuritySaveType(store_type);
   }
   else if  (command ==   SetSecurityNbOfEventsCmd){
   	pAnalysisManager
        	->SetNeventsSecuritySave(SetSecurityNbOfEventsCmd->
	                                     	GetNewIntValue(newValues));
   }
  /* else if (command ==   SetPrimaryIncidentFluxCmd){
  	pAnalysisManager->SetPrimaryIntegralFlux(SetPrimaryIncidentFluxCmd->
                                                  GetNewDoubleValue(newValues));
   }*/

   else if (command == SetBugVerboseCmd){
   	G4int verbosity = SetBugVerboseCmd->GetNewIntValue(newValues);
	pAnalysisManager->SetBugVerbose(verbosity);
	const G4UserStackingAction* theStackingAction = 
     		G4RunManager::GetRunManager()->GetUserStackingAction();
  	G4UserStackingAction* theStackingAction1 =    
           const_cast<G4UserStackingAction*> (theStackingAction);
 	(dynamic_cast<PLANETOCOSStackingAction*> (theStackingAction1))
						->SetBugVerbose(verbosity);	 

	
   }
  
        
}

