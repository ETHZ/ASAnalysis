#ifndef DiPhotonMiniTree_hh
#define DiPhotonMiniTree_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

#include "TLorentzVector.h"
#include "TParameter.h"

#include "EnergyCorrection.hh"

#include <map>
#include <utility> 
#include <vector>
#include <set>

#include "TRandom3.h"
#include "helper/Utilities.hh"

#include "TTree.h"
#include "TFile.h"
#include "Math/Vector3D.h"
#include "TVector3.h"

#include "TGeoPara.h"
#include "TGeoTube.h"

enum SigBkgMode {
  kSigSig = 0,
  kSigBkg = 1,
  kBkgSig = 2,
  kBkgBkg = 3,
  kSig = 4,
  kBkg = 5
};

enum ChoiceMixingTemplates {
  k1Event = 0,
  k2Events = 1
};

enum MatchingStatus {
  kSignal = 0,
  kBackground = 1,
  kElectron = 2
};

typedef enum {
  k2Dstandard_selection=0,
  k1Drandomcone_template=1,
  k1Dsideband_template=2,
  k2DZee_pixelvetoreversed_selection=3,
  k1Dpreselection=4,
  k1Dselection=5,
  k2Drandomcone_template=6,
  k2Drandomconesideband_template=7,
  k2Dsideband_template=8,
  k2Dstandard_preselection=9,
  k2DZmumu_selection=10,
  k1Dsignal_template=11,
  k1Dbackground_template=12,
  k2Dtruesigsig_template=13,
  k2Dtruesigbkg_template=14,
  k2Dtruebkgbkg_template=15,
  k2Drconeplusgenfake_template=16,
  k2Dgenpromptplussideband_template=17,
  kLightDefault=100,
  kLightJECup=101,
  kLightJECdown=102,
  kLightJERup=103,
  kLightJERdown=104,
  kLightESCALEup=105,
  kLightESCALEdown=106,
  kLightESMEARup=107,
  kLightESMEARdown=108,
  kLightDiElectron=200
} EnumSel;

typedef enum {
  kTurnOff=-1,
  kData=0,
  kMC=1,
  kBoth=2
} DataOrMCSwitch;

typedef enum {
  kUncorrected = -999,
  kCorrectedNoShift = 0,
  kShiftUp = 1,
  kShiftDown = -1
} StatusScaleUpScaleDown;

class Selection {

public:

  Selection(EnumSel id_, TString name_, DataOrMCSwitch dataormc_, int islight_, TFile *OutputFile, TFile *OutputExtraFile=NULL): id(id_), name(name_){
    switch (dataormc_){
    case kTurnOff:
      runondata=false;
      runonmc=false;
      break;
    case kData:
      runondata=true;
      runonmc=false;
      break;
    case kMC:
      runondata=false;
      runonmc=true;
      break;
    case kBoth:
      runondata=true;
      runonmc=true;
      break;
    default:
      runondata=false;
      runonmc=false;
      break;
    }
    isfullsel=!islight_;
    islightsel=islight_;
    tree_full=NULL;
    tree_extra=NULL;
    tree_light=NULL;

    if (isfullsel){
      if (OutputFile){
	OutputFile->cd();
	TString tname = Form("Tree_%s",name.Data());
	tree_full = new TTree(tname.Data(),tname.Data());
      }
      if (OutputExtraFile){
	OutputExtraFile->cd();
	TString tname = Form("Tree_%s_EXTRA",name.Data());
	tree_extra = new TTree(tname.Data(),tname.Data());
      }
    }
    if (islightsel){
      if (OutputFile){
        OutputFile->cd();
	TString tname = Form("LightTreeGenReco_%s",name.Data());
	tree_light = new TTree(tname.Data(),tname.Data());
      }
    }

    is2d = (name.First("2")>-1);
  };
  ~Selection(){};

  bool operator==(const EnumSel& rhs){ return (this->id==rhs); }
  bool operator!=(const EnumSel& rhs){ return !(*this==rhs); }

  EnumSel id;
  TString name;
  bool runondata;
  bool runonmc;
  bool isfullsel;
  bool islightsel;
  TTree *tree_full;
  TTree *tree_extra;
  TTree *tree_light;
  bool is2d;

};


const int global_maxN_photonpfcandidates = 2000;
const int global_maxN_vetoobjects = 200;
const int global_maxN_jets = 200;

const int global_size_pfcandarrays = 30;

const float global_dR_cut_acceptance = 0.45;
// const float global_mindR_photon_jet = 1.0; unused, is implemented at a later stage to keep all jets for the inclusive case

const bool do_only_light_tree = false; // set this to true to run only the efficiency/unfolding light tree
const int dy_dataset_id = 30; // be SURE to know where to propagate this if you change it
const int sherpa_dataset_id = -30; // be SURE to know where to propagate this if you change it

const int nclosest = 5;
const int nclosest_inputmatching = 40;

// JET PRE-SELECTION (additional selection applied in template_production)
const float min_jet_pt_alljets = 20;
const float max_eta_jets = 4.7;

// CAREFUL if you set this to true: the selection of saved pfcandidates,
// for instance, is decided at ntuple production level.
// Check that you really don't have to re-run the ntuple as well.
const bool do_recalc_isolation = false;
// ------------------------------------

typedef struct {
  float photon;
  float charged;
  float neutral;
  int nphotoncand;
  int nchargedcand;
  int nneutralcand;
  std::vector<float> photoncandenergies;
  std::vector<float> chargedcandenergies;
  std::vector<float> neutralcandenergies;
  std::vector<float> photoncandets;
  std::vector<float> chargedcandets;
  std::vector<float> neutralcandets;
  std::vector<float> photoncanddetas;
  std::vector<float> chargedcanddetas;
  std::vector<float> neutralcanddetas;
  std::vector<float> photoncanddphis;
  std::vector<float> chargedcanddphis;
  std::vector<float> neutralcanddphis;
  float newphi;
} isolations_struct;

typedef struct {
  float dR;
  float dEta;
  float dPhi;
} angular_distances_struct;

typedef struct {
  std::vector<float> PfCandPt;
  std::vector<float> PfCandEta;
  std::vector<float> PfCandPhi;
  std::vector<float> PfCandVx;
  std::vector<float> PfCandVy;
  std::vector<float> PfCandVz;
} pfcandidates_struct;

typedef struct {
  int m_jet;
  float phopt_footprint_total;
  float phopt_footprint_m_frac;
  float jetpt_pf;
  float jetpt_m_frac;
  float jetpt_m_frac_PhoComp;
} jetmatching_struct;

typedef struct {
  float meta;
  float Meta;
  float mr9;
  float Mr9;
  int mrun;
  int Mrun;
  float val;
  float err;
} struct_escale_item;


class DiPhotonMiniTree : public UserAnalysisBase{
public:
  DiPhotonMiniTree(TreeReader *tr = NULL, std::string dataType="data", Float_t aw=-999, Float_t* _kfac=NULL, Float_t _minthrpfphotoncandEB=0, Float_t _minthrpfphotoncandEE=0, bool _isstep2 = false, TString _input_filename = "", UInt_t _uuid=0, int year=-1, int dataset_id=0);
  virtual ~DiPhotonMiniTree();

  void Begin();
  void Analyze();
  void End();

private:

  int debug;

  float global_linkbyrechit_enlargement;
  float global_minthrpfphotoncandEB;
  float global_minthrpfphotoncandEE;

  int dataset_id;
  int year;
  bool global_is2011;
  bool global_is2012;

  TGeoPara eegeom;

  Float_t* kfactors;

  float getEtaCorrectionBarrel(float eta);
  void CorrPhoton(TreeReader *fTR, int i, std::vector<float> const &unscaled_r9, std::vector<float> const &unscaled_energy, StatusScaleUpScaleDown status_escale, StatusScaleUpScaleDown status_esmear);
  void CorrectPhotonEnergy(StatusScaleUpScaleDown status_escale, StatusScaleUpScaleDown status_esmear, std::vector<std::pair<int,float> > &neworder, std::vector<float> const &unscaled_r9, std::vector<float> const &unscaled_energy);
  StatusScaleUpScaleDown CorrectPhotonR9(StatusScaleUpScaleDown status, std::vector<float> &unscaled_r9);
  StatusScaleUpScaleDown CorrectPhotonSieie(StatusScaleUpScaleDown status, std::vector<float> &unscaled_sieie);
  void CorrPhoton(TreeReader *fTR, int i, std::vector<float> *unscaled_r9);
	
  std::vector<int> ApplyPixelVeto(TreeReader *fTR, vector<int> passing, bool forelectron=0);
  std::vector<int> PhotonPreSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> GenLevelIsolationCut(TreeReader *fTR, vector<int> passing);
  std::vector<int> PhotonSelection(TreeReader *fTR, vector<int> passing, TString mode="");
  std::vector<int> SignalSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> BackgroundSelection(TreeReader *fTR, vector<int> passing);
  std::vector<int> ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert=false);
  std::vector<int> GetPFCandRemovals(TreeReader *fTR, int phoqi);
  bool SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> &passing);
  bool StandardEventSelection(TreeReader *fTR, std::vector<int> &passing, std::vector<int> &passing_jets);
  bool VetoJetPhotonOverlap(std::vector<int> &passing, std::vector<int> &passing_jets);
  void JetSelection(std::vector<int> &passing_jets);
  void GenJetSelection(std::vector<int> &passing_gen_jets);
  bool TriggerSelection();
  int Count_part_isrfsr_gamma(TreeReader *fTR);
  void ResetVars();
  void Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex);  
  std::set<int> GetPFCandInsideFootprint(TreeReader *fTR, int phoqi, float rotation_phi, TString component);
  std::set<int> GetPFCandInsideFootprint(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, TString component);
  std::set<int> GetPrecalculatedFootprintPhoEl(int phoqi);
  std::set<int> GetPFCandWithFootprintRemoval(TreeReader *fTR, int phoqi, float rotation_phi, bool outoffootprint, TString component);
  std::set<int> GetPFCandWithFootprintRemoval(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, bool outoffootprint, TString component);
  TVector3 PropagatePFCandToEcal(TreeReader *fTR, int pfcandindex, float position, bool isbarrel);
  TVector3 PropagatePFCandToEcal(pfcandidates_struct *pfcands, int pfcandindex, float position, bool isbarrel);
  bool FindImpingingTrack(TreeReader *fTR, int phoqi, int &reference_index_found, bool dofootprintremoval = false, std::set<int> removals = std::set<int>());
  float PFIsolation(int phoqi, float rotation_phi, TString component, int *counter = NULL, std::vector<float> *energies = NULL, std::vector<float> *ets = NULL,  std::vector<float> *detas = NULL, std::vector<float> *dphis = NULL, float *newphi = NULL, std::set<int> removals = std::set<int>());
  std::pair<float,float> PFPhotonIsolationFromMinitree(int phoqi1, int phoqi2, pfcandidates_struct *pfcands, bool doremoval1 = true, bool doremoval2 = true, float matched_eta1=-999, float matched_eta2=-999);
  angular_distances_struct GetPFCandDeltaRFromSC(TreeReader *fTR, int phoqi, int pfindex, float rotation_phi = 0);
  angular_distances_struct GetPFCandDeltaRFromSC(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, int pfindex, float matched_eta=-999);
  bool FindCloseJetsAndPhotons(TreeReader *fTR, float rotation_phi, int phoqi, TString mod="");
  bool FindCloseJetsAndPhotons(std::vector<std::pair<float,float> > obj, float eta, float phi);
  std::set<int> GetPFCandIDedRemovals(TreeReader *fTR, int phoqi);

  void FillPhoIso_NewTemplates(TreeReader *fTR, Int_t *n1_arr, Int_t *n2_arr, std::vector<int> passing, SigBkgMode mode, ChoiceMixingTemplates mixing);
  void FillVetoObjects(TreeReader *fTR, int phoqi, TString mod);
  void InitInputTree();

  std::vector<jetmatching_struct> photon_jet_matching;

  std::vector<int> DiPhotonInvariantMassCutSelection(TreeReader *fTR, std::vector<int> passing);

  Int_t Choose_bin_eta(float eta, int region);

  float SieieRescale(float sieie, bool isbarrel);
  float R9Rescale(float r9, bool isbarrel);
  float CalculateSCArea(TreeReader *fTR, int scindex);
  float GetPUEnergy(TreeReader *fTR, TString mode, float eta);

  std::vector<struct_escale_item> energyScaleDatabase;
  std::vector<struct_escale_item> energySmearingDatabase;
  float EnergyScaleOffset(float eta, float r9, int run, float sigmas);
  float EnergySmearingCorrection(float eta, float r9, int run, float sigmas);
  void InitEnergyScalesAndSmearingsDatabase();
  void InsertEnergyScaleItem(float meta, float Meta, float mr9, float Mr9, int mrun, int Mrun, float val, float err);
  void InsertEnergySmearingItem(float meta, float Meta, float mr9, float Mr9, int mrun, int Mrun, float val, float err);

  void FillLead(int index, std::vector<int> passing_jets);
  void FillTrail(int index, std::vector<int> passing_jets);
  void FillJetsInfo(std::vector<int> passing, std::vector<int> passing_jets);
  void FillGenJetsInfo(std::vector<int> passing_gen, std::vector<int> passing_gen_jets);

  //  double etaTransformation(float EtaParticle, float Zvertex);
  double phiNorm(float phi);
  isolations_struct RandomConeIsolation(TreeReader *fTR, int phoqi, TString mod="");
  isolations_struct PFConeIsolation(TreeReader *fTR, int phoqi);
  int FindPFCandType(int id);

  MatchingStatus check_matching_status(int phoindex);
  MatchingStatus determine_matchingstatus(int status, float geniso);
  void FixMatchingStatusElectrons(int phoindex);
  void FixMatchingStatusSherpa(int phoindex);

  int CountChargedHadronsInCone(TreeReader *fTR, int phoqi, std::set<int> removals = std::set<int>(), bool skipvetocones=false);
  std::set<int> NChargedHadronsInConeSelection(TreeReader *fTR, std::vector<int> passing, int minimum=0, int maximum=9999);

  std::vector<int> MuonSelection(TreeReader *fTR, std::vector<int> passing);
  bool DiMuonFromZSelection(TreeReader *fTR, std::vector<int> &passing);
  float PFPhotonIsolationAroundMuon(int muqi, int *counter, std::vector<float> *energies = NULL, std::vector<float> *ets = NULL,  std::vector<float> *detas = NULL, std::vector<float> *dphis = NULL);
  void FillMuonInfo(int index);

  bool PassPrimaryVertexFilter();

  float DeltaPhiSigned(float phi1, float phi2);

  typedef std::pair<int,float> OrderPair;
  struct IndexByPt {
    const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
      return j1.second > j2.second;
    }
  };
  IndexByPt indexComparator;

  jetmatching_struct PFMatchPhotonToJet(int phoqi);

  OnTheFlyCorrections *MyJetCorrector;
  void JECJERCorrection(StatusScaleUpScaleDown e_scale, StatusScaleUpScaleDown e_smear, vector<pair<int,float> > &ordering_jets);

  TRandom3 *randomgen;
  TRandom3 *randomgen_forEsmearing;

  static const int n_templates_EB=7;
  static const int n_templates_EE=5;

  std::vector<float> binsdef_single_gamma_EB_eta;
  std::vector<float> binsdef_single_gamma_EE_eta;
  std::vector<float> eff_areas_EB_photon_data;
  std::vector<float> eff_areas_EE_photon_data;
  std::vector<float> eff_areas_EB_charged_data;
  std::vector<float> eff_areas_EE_charged_data;
  std::vector<float> eff_areas_EB_neutral_data;
  std::vector<float> eff_areas_EE_neutral_data;
  std::vector<float> eff_areas_EB_photon_MC;
  std::vector<float> eff_areas_EE_photon_MC;
  std::vector<float> eff_areas_EB_charged_MC;
  std::vector<float> eff_areas_EE_charged_MC;
  std::vector<float> eff_areas_EB_neutral_MC;
  std::vector<float> eff_areas_EE_neutral_MC;

  bool inputtree_isinitialized;

  float scarea[100];
  float scareaSF[100];
  
  std::string fDataType_;
  bool isdata;

  TH1F *histo_PFPhotonDepositAroundImpingingTrack;

  Float_t AddWeight;

  map<EnumSel,Selection> sels;
  TFile* fOutputExtraFile;

  Float_t event_luminormfactor;
  Float_t event_Kfactor;
  Float_t event_weight;
  Float_t event_weight3D;

  Float_t event_rho;
  Float_t event_sigma;
  Int_t event_nPU;
  Int_t event_nPUtrue;
  Int_t event_nRecVtx;
  Int_t event_pass12whoissiglike;
  
  Float_t dipho_mgg_photon;
  
  Float_t pholead_eta, photrail_eta;
  Float_t pholead_phi, photrail_phi;
  Float_t pholead_pt, photrail_pt;
  Float_t pholead_energy, photrail_energy;

  Float_t pholead_SCeta, photrail_SCeta;
  Float_t pholead_SCphi, photrail_SCphi;
 
  Int_t pholead_PhoHasPixSeed, photrail_PhoHasPixSeed;

  Float_t pholead_r9, photrail_r9;
  Float_t pholead_sieie, photrail_sieie;
  Float_t pholead_hoe, photrail_hoe;

  Float_t pholead_PhoSCRemovalPFIsoCharged;
  Float_t photrail_PhoSCRemovalPFIsoCharged;
  Float_t pholead_PhoSCRemovalPFIsoNeutral;
  Float_t photrail_PhoSCRemovalPFIsoNeutral;
  Float_t pholead_PhoSCRemovalPFIsoPhoton;
  Float_t photrail_PhoSCRemovalPFIsoPhoton;
  Float_t pholead_PhoSCRemovalPFIsoCombined;
  Float_t photrail_PhoSCRemovalPFIsoCombined;

  Float_t pholead_PhoIso03Ecal;
  Float_t pholead_PhoIso03Hcal;
  Float_t pholead_PhoIso03TrkHollow;
  Float_t photrail_PhoIso03Ecal;
  Float_t photrail_PhoIso03Hcal;
  Float_t photrail_PhoIso03TrkHollow;
  
  Int_t pholead_PhoPassConversionVeto;
  Int_t photrail_PhoPassConversionVeto;
  Float_t pholead_GenPhotonIsoDR04;
  Float_t photrail_GenPhotonIsoDR04;
 
  Int_t pholead_PhoMCmatchexitcode, photrail_PhoMCmatchexitcode;

  Float_t pholead_scarea, photrail_scarea;
  Float_t pholead_scareaSF, photrail_scareaSF;

  Float_t pholead_m_jet_ptcorr;
  Float_t pholead_m_jet_dR;
  Float_t pholead_phopt_footprint_total;
  Float_t pholead_phopt_footprint_m_frac;
  Float_t pholead_jetpt_pf;
  Float_t pholead_jetpt_m_frac;
  Float_t pholead_jetpt_m_frac_PhoComp;
  Float_t photrail_m_jet_ptcorr;
  Float_t photrail_m_jet_dR;
  Float_t photrail_phopt_footprint_total;
  Float_t photrail_phopt_footprint_m_frac;
  Float_t photrail_jetpt_pf;
  Float_t photrail_jetpt_m_frac;
  Float_t photrail_jetpt_m_frac_PhoComp;

  Float_t pholead_pt_closestjet;
  Float_t pholead_dR_closestjet;
  Float_t photrail_pt_closestjet;
  Float_t photrail_dR_closestjet;

  Int_t pholead_Npfcandphotonincone;
  Int_t pholead_Npfcandchargedincone;
  Int_t pholead_Npfcandneutralincone;
  Int_t photrail_Npfcandphotonincone;
  Int_t photrail_Npfcandchargedincone;
  Int_t photrail_Npfcandneutralincone;

  Float_t pholead_photonpfcandenergies[global_size_pfcandarrays];
  Float_t pholead_chargedpfcandenergies[global_size_pfcandarrays];
  Float_t pholead_neutralpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_photonpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_chargedpfcandenergies[global_size_pfcandarrays];
  Float_t photrail_neutralpfcandenergies[global_size_pfcandarrays];

  Float_t pholead_photonpfcandets[global_size_pfcandarrays];
  Float_t pholead_chargedpfcandets[global_size_pfcandarrays];
  Float_t pholead_neutralpfcandets[global_size_pfcandarrays];
  Float_t photrail_photonpfcandets[global_size_pfcandarrays];
  Float_t photrail_chargedpfcandets[global_size_pfcandarrays];
  Float_t photrail_neutralpfcandets[global_size_pfcandarrays];

  Float_t pholead_photonpfcanddetas[global_size_pfcandarrays];
  Float_t pholead_chargedpfcanddetas[global_size_pfcandarrays];
  Float_t pholead_neutralpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_photonpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_chargedpfcanddetas[global_size_pfcandarrays];
  Float_t photrail_neutralpfcanddetas[global_size_pfcandarrays];

  Float_t pholead_photonpfcanddphis[global_size_pfcandarrays];
  Float_t pholead_chargedpfcanddphis[global_size_pfcandarrays];
  Float_t pholead_neutralpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_photonpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_chargedpfcanddphis[global_size_pfcandarrays];
  Float_t photrail_neutralpfcanddphis[global_size_pfcandarrays];

  TH1F *fHNumPU;
  TH1F *fHNumPU_noweight;
  TH1F *fHNumPUTrue;
  TH1F *fHNumPUTrue_noweight;
  TH1F *fHNumVtx;
  TParameter<double> *fPNumEvents;
  TParameter<double> *fPSumWeights;
  TParameter<double> *fPSumWeights2;

  UInt_t event_fileuuid;

  Int_t event_run;
  Int_t event_lumi;
  UInt_t event_number;
  Int_t mydataset_id;

  Float_t pholead_GEN_eta, photrail_GEN_eta;
  Float_t pholead_GEN_phi, photrail_GEN_phi;
  Float_t pholead_GEN_px, photrail_GEN_px;
  Float_t pholead_GEN_py, photrail_GEN_py;
  Float_t pholead_GEN_pz, photrail_GEN_pz;
  Float_t pholead_GEN_pt, photrail_GEN_pt;
  Float_t pholead_GEN_energy, photrail_GEN_energy;

  Int_t n_GEN_jets;
  Float_t jet_GEN_pt[global_maxN_jets];
  Float_t jet_GEN_eta[global_maxN_jets];
  Float_t jet_GEN_phi[global_maxN_jets];
  Float_t jet_GEN_energy[global_maxN_jets];

  Bool_t tree_gen_in_acc;
  Bool_t tree_reco_in_acc;
  Bool_t tree_matched;

  Float_t pholead_test_rotatedphotoniso[50];
  Float_t pholead_test_rotatedwithcheckphotoniso[50];

  Int_t allphotonpfcand_count;
  Float_t allphotonpfcand_pt[global_maxN_photonpfcandidates];
  Float_t allphotonpfcand_eta[global_maxN_photonpfcandidates];
  Float_t allphotonpfcand_phi[global_maxN_photonpfcandidates];
  Float_t allphotonpfcand_vx[global_maxN_photonpfcandidates];
  Float_t allphotonpfcand_vy[global_maxN_photonpfcandidates];
  Float_t allphotonpfcand_vz[global_maxN_photonpfcandidates];

  Float_t phoiso_template_1event_sigsig_1[nclosest];
  Float_t phoiso_template_1event_sigsig_2[nclosest];
  Float_t phoiso_template_1event_sigbkg_1[nclosest];
  Float_t phoiso_template_1event_sigbkg_2[nclosest];
  Float_t phoiso_template_1event_bkgsig_1[nclosest];
  Float_t phoiso_template_1event_bkgsig_2[nclosest];
  Float_t phoiso_template_1event_bkgbkg_1[nclosest];
  Float_t phoiso_template_1event_bkgbkg_2[nclosest];
  Float_t phoiso_template_2events_sigsig_1[nclosest];
  Float_t phoiso_template_2events_sigsig_2[nclosest];
  Float_t phoiso_template_2events_sigbkg_1[nclosest];
  Float_t phoiso_template_2events_sigbkg_2[nclosest];
  Float_t phoiso_template_2events_bkgsig_1[nclosest];
  Float_t phoiso_template_2events_bkgsig_2[nclosest];
  Float_t phoiso_template_2events_bkgbkg_1[nclosest];
  Float_t phoiso_template_2events_bkgbkg_2[nclosest];

  Float_t rewinfo_template_1event_sigsig_1[nclosest*6];
  Float_t rewinfo_template_1event_sigsig_2[nclosest*6];
  Float_t rewinfo_template_1event_sigbkg_1[nclosest*6];
  Float_t rewinfo_template_1event_sigbkg_2[nclosest*6];
  Float_t rewinfo_template_1event_bkgsig_1[nclosest*6];
  Float_t rewinfo_template_1event_bkgsig_2[nclosest*6];
  Float_t rewinfo_template_1event_bkgbkg_1[nclosest*6];
  Float_t rewinfo_template_1event_bkgbkg_2[nclosest*6];
  Float_t rewinfo_template_2events_sigsig_1[nclosest*6];
  Float_t rewinfo_template_2events_sigsig_2[nclosest*6];
  Float_t rewinfo_template_2events_sigbkg_1[nclosest*6];
  Float_t rewinfo_template_2events_sigbkg_2[nclosest*6];
  Float_t rewinfo_template_2events_bkgsig_1[nclosest*6];
  Float_t rewinfo_template_2events_bkgsig_2[nclosest*6];
  Float_t rewinfo_template_2events_bkgbkg_1[nclosest*6];
  Float_t rewinfo_template_2events_bkgbkg_2[nclosest*6];

  Int_t   vetoobjects_count;
  Float_t vetoobjects_pt[global_maxN_vetoobjects];
  Float_t vetoobjects_eta[global_maxN_vetoobjects];
  Float_t vetoobjects_phi[global_maxN_vetoobjects];
  Int_t vetoobjects_type[global_maxN_vetoobjects];

  Int_t n_jets;
  Float_t jet_pt[global_maxN_jets];
  Float_t jet_eta[global_maxN_jets];
  Float_t jet_phi[global_maxN_jets];
  Float_t jet_energy[global_maxN_jets];

  UInt_t  input_event_number;
  Int_t   input_allphotonpfcand_count;
  Float_t input_allphotonpfcand_pt[global_maxN_photonpfcandidates];
  Float_t input_allphotonpfcand_eta[global_maxN_photonpfcandidates];
  Float_t input_allphotonpfcand_phi[global_maxN_photonpfcandidates];
  Float_t input_allphotonpfcand_vx[global_maxN_photonpfcandidates];
  Float_t input_allphotonpfcand_vy[global_maxN_photonpfcandidates];
  Float_t input_allphotonpfcand_vz[global_maxN_photonpfcandidates];
  Float_t input_pholead_SCeta;
  Float_t input_pholead_SCphi;
  Float_t input_pholead_pt;
  Float_t input_photrail_SCeta;
  Float_t input_photrail_SCphi;
  Float_t input_photrail_pt;
  Float_t input_event_rho;
  Float_t input_event_sigma;
  Int_t   input_event_pass12whoissiglike;
  Int_t   input_vetoobjects_count;
  Float_t input_vetoobjects_eta[global_maxN_vetoobjects];
  Float_t input_vetoobjects_phi[global_maxN_vetoobjects];

  TBranch *b_input_event_number;
  TBranch *b_input_allphotonpfcand_count;
  TBranch *b_input_allphotonpfcand_pt   ;
  TBranch *b_input_allphotonpfcand_eta  ;
  TBranch *b_input_allphotonpfcand_phi  ;
  TBranch *b_input_allphotonpfcand_vx   ;
  TBranch *b_input_allphotonpfcand_vy   ;
  TBranch *b_input_allphotonpfcand_vz   ;
  TBranch *b_input_pholead_SCeta;
  TBranch *b_input_pholead_SCphi;
  TBranch *b_input_pholead_pt;
  TBranch *b_input_photrail_SCeta;
  TBranch *b_input_photrail_SCphi;
  TBranch *b_input_photrail_pt;
  TBranch *b_input_event_rho;
  TBranch *b_input_event_sigma;
  TBranch *b_input_event_pass12whoissiglike;
  TBranch *b_input_vetoobjects_count;
  TBranch *b_input_vetoobjects_eta;
  TBranch *b_input_vetoobjects_phi;

  TBranch *b_matchingtree_event_fileuuid;

  TBranch *b_matchingtree_event_run;
  TBranch *b_matchingtree_event_lumi;
  TBranch *b_matchingtree_event_number;
  TBranch *b_matchingtree_index_1event_sigsig_1;
  TBranch *b_matchingtree_index_1event_sigsig_2;
  TBranch *b_matchingtree_index_1event_sigbkg_1;
  TBranch *b_matchingtree_index_1event_sigbkg_2;
  TBranch *b_matchingtree_index_1event_bkgsig_1;
  TBranch *b_matchingtree_index_1event_bkgsig_2;
  TBranch *b_matchingtree_index_1event_bkgbkg_1;
  TBranch *b_matchingtree_index_1event_bkgbkg_2;
  TBranch *b_matchingtree_index_2events_sigsig_1;
  TBranch *b_matchingtree_index_2events_sigsig_2;
  TBranch *b_matchingtree_index_2events_sigbkg_1;
  TBranch *b_matchingtree_index_2events_sigbkg_2;
  TBranch *b_matchingtree_index_2events_bkgsig_1;
  TBranch *b_matchingtree_index_2events_bkgsig_2;
  TBranch *b_matchingtree_index_2events_bkgbkg_1;
  TBranch *b_matchingtree_index_2events_bkgbkg_2;


  bool isstep2;
  TFile *f_input;
  TString input_filename;
  UInt_t uuid;
  TTree *InputTree[2];
  TTree *matchingtree;
  Int_t matchingtree_event_run;
  Int_t matchingtree_event_lumi;
  UInt_t matchingtree_event_number;
  Int_t matchingtree_index_1event_sigsig_1[nclosest_inputmatching];
  Int_t matchingtree_index_1event_sigsig_2[nclosest_inputmatching];
  Int_t matchingtree_index_1event_sigbkg_1[nclosest_inputmatching];
  Int_t matchingtree_index_1event_sigbkg_2[nclosest_inputmatching];
  Int_t matchingtree_index_1event_bkgsig_1[nclosest_inputmatching];
  Int_t matchingtree_index_1event_bkgsig_2[nclosest_inputmatching];
  Int_t matchingtree_index_1event_bkgbkg_1[nclosest_inputmatching];
  Int_t matchingtree_index_1event_bkgbkg_2[nclosest_inputmatching];
  Int_t matchingtree_index_2events_sigsig_1[nclosest_inputmatching];
  Int_t matchingtree_index_2events_sigsig_2[nclosest_inputmatching];
  Int_t matchingtree_index_2events_sigbkg_1[nclosest_inputmatching];
  Int_t matchingtree_index_2events_sigbkg_2[nclosest_inputmatching];
  Int_t matchingtree_index_2events_bkgsig_1[nclosest_inputmatching];
  Int_t matchingtree_index_2events_bkgsig_2[nclosest_inputmatching];
  Int_t matchingtree_index_2events_bkgbkg_1[nclosest_inputmatching];
  Int_t matchingtree_index_2events_bkgbkg_2[nclosest_inputmatching];



};
#endif
