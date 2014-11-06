#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType, Float_t aw, Float_t* _kfac, Float_t _minthrpfphotoncandEB, Float_t _minthrpfphotoncandEE, bool _isstep2, TString _input_filename, UInt_t _uuid, int _year, int dataset_id_) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw), kfactors(_kfac), global_minthrpfphotoncandEB(_minthrpfphotoncandEB), global_minthrpfphotoncandEE(_minthrpfphotoncandEE), isstep2(_isstep2), input_filename(_input_filename), uuid(_uuid), year(_year), dataset_id(dataset_id_){
  Util::SetStyle();	
  if (fDataType_ == "mc") isdata=false;
  else if (fDataType_ == "data") isdata=true; 
  else {
    std::cout << "wrong data type" << std::endl;
    assert(1==0);
  }
  randomgen = new TRandom3(0);
  randomgen_forEsmearing = new TRandom3(0);

  global_linkbyrechit_enlargement = 0.25; // xtal_size_eff = (1+global_linkbyrechit_enlargement)*xtal_size

  eegeom = TGeoPara(1,1,1,0,0,0);

  global_is2011=false;
  global_is2012=false;
  if (year==2011) global_is2011=true;
  if (year==2012) global_is2012=true;

  if (isdata) assert(dataset_id==0);
  else assert(dataset_id!=0);

  debug = 0;

}

DiPhotonMiniTree::~DiPhotonMiniTree(){
  delete randomgen;
  delete randomgen_forEsmearing;
}

void DiPhotonMiniTree::Begin(){

  cout << "Begin" << endl;

  fOutputExtraFile = (isdata && !isstep2) ? new TFile(Form("extrainfo_%u.root",uuid),"recreate") : NULL;

  sels.insert(pair<EnumSel,Selection>(k2Dstandard_selection,Selection(k2Dstandard_selection,"2Dstandard_selection",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Drandomcone_template,Selection(k1Drandomcone_template,"1Drandomcone_template",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Dsideband_template,Selection(k1Dsideband_template,"1Dsideband_template",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2DZee_pixelvetoreversed_selection,Selection(k2DZee_pixelvetoreversed_selection,"2DZee_pixelvetoreversed_selection",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Dpreselection,Selection(k1Dpreselection,"1Dpreselection",kTurnOff,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Dselection,Selection(k1Dselection,"1Dselection",kTurnOff,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Drandomcone_template,Selection(k2Drandomcone_template,"2Drandomcone_template",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Drandomconesideband_template,Selection(k2Drandomconesideband_template,"2Drandomconesideband_template",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dsideband_template,Selection(k2Dsideband_template,"2Dsideband_template",kBoth,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dstandard_preselection,Selection(k2Dstandard_preselection,"2Dstandard_preselection",kTurnOff,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2DZmumu_selection,Selection(k2DZmumu_selection,"2DZmumu_selection",kTurnOff,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Dsignal_template,Selection(k1Dsignal_template,"1Dsignal_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k1Dbackground_template,Selection(k1Dbackground_template,"1Dbackground_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dtruesigsig_template,Selection(k2Dtruesigsig_template,"2Dtruesigsig_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dtruesigbkg_template,Selection(k2Dtruesigbkg_template,"2Dtruesigbkg_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dtruebkgbkg_template,Selection(k2Dtruebkgbkg_template,"2Dtruebkgbkg_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Drconeplusgenfake_template,Selection(k2Drconeplusgenfake_template,"2Drconeplusgenfake_template",kMC,false,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(k2Dgenpromptplussideband_template,Selection(k2Dgenpromptplussideband_template,"2Dgenpromptplussideband_template",kMC,false,fOutputFile,fOutputExtraFile)));

  sels.insert(pair<EnumSel,Selection>(kLightDefault,Selection(kLightDefault,"Default",kBoth,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightJECup,Selection(kLightJECup,"JECup",kBoth,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightJECdown,Selection(kLightJECdown,"JECdown",kBoth,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightJERup,Selection(kLightJERup,"JERup",kMC,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightJERdown,Selection(kLightJERdown,"JERdown",kMC,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightESCALEup,Selection(kLightESCALEup,"ESCALEup",kData,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightESCALEdown,Selection(kLightESCALEdown,"ESCALEdown",kData,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightESMEARup,Selection(kLightESMEARup,"ESMEARup",kMC,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightESMEARdown,Selection(kLightESMEARdown,"ESMEARdown",kMC,true,fOutputFile,fOutputExtraFile)));
  sels.insert(pair<EnumSel,Selection>(kLightDiElectron,Selection(kLightDiElectron,"DiElectron",kMC,true,fOutputFile,fOutputExtraFile)));


  for (map<EnumSel,Selection>::iterator it = sels.begin(); it!=sels.end(); it++){

    if (!it->second.isfullsel) continue;

    TTree *thistree = it->second.tree_full;
    TTree *thisextratree = it->second.tree_extra;

    thistree->Branch("event_fileuuid",&event_fileuuid,"event_fileuuid/i");
    thistree->Branch("event_run",&event_run,"event_run/I");
    thistree->Branch("event_lumi",&event_lumi,"event_lumi/I");
    thistree->Branch("event_number",&event_number,"event_number/i");
    if (thisextratree) thisextratree->Branch("event_fileuuid",&event_fileuuid,"event_fileuuid/i");
    if (thisextratree) thisextratree->Branch("event_run",&event_run,"event_run/I");
    if (thisextratree) thisextratree->Branch("event_lumi",&event_lumi,"event_lumi/I");
    if (thisextratree) thisextratree->Branch("event_number",&event_number,"event_number/i");
    thistree->Branch("dataset_id",&mydataset_id,"dataset_id/I");
    if (thisextratree) thisextratree->Branch("dataset_id",&mydataset_id,"dataset_id/I");

    thistree->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
    if (thisextratree) thisextratree->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
    thistree->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");
    if (thisextratree) thisextratree->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");

    thistree->Branch("event_weight",&event_weight,"event_weight/F");
    if (thisextratree) thisextratree->Branch("event_weight",&event_weight,"event_weight/F");
    thistree->Branch("event_rho",&event_rho,"event_rho/F");
    if (thisextratree) thisextratree->Branch("event_rho",&event_rho,"event_rho/F");
    thistree->Branch("event_sigma",&event_sigma,"event_sigma/F");
    if (thisextratree) thisextratree->Branch("event_sigma",&event_sigma,"event_sigma/F");
    thistree->Branch("event_nPU",&event_nPU,"event_nPU/I");
    thistree->Branch("event_nPUtrue",&event_nPUtrue,"event_nPUtrue/I");
    thistree->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");
    thistree->Branch("event_pass12whoissiglike",&event_pass12whoissiglike,"event_pass12whoissiglike/I");
    if (thisextratree) thisextratree->Branch("event_pass12whoissiglike",&event_pass12whoissiglike,"event_pass12whoissiglike/I");

    thistree->Branch("dipho_mgg_photon",&dipho_mgg_photon,"dipho_mgg_photon/F");

    thistree->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
    if (thisextratree) thisextratree->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
    thistree->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
    if (thisextratree) thisextratree->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
    thistree->Branch("pholead_phi",&pholead_phi,"pholead_phi/F");
    if (thisextratree) thisextratree->Branch("pholead_phi",&pholead_phi,"pholead_phi/F");
    thistree->Branch("photrail_phi",&photrail_phi,"photrail_phi/F");
    if (thisextratree) thisextratree->Branch("photrail_phi",&photrail_phi,"photrail_phi/F");
    thistree->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
    if (thisextratree) thisextratree->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
    thistree->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
    if (thisextratree) thisextratree->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
    thistree->Branch("pholead_energy",&pholead_energy,"pholead_energy/F");
    if (thisextratree) thisextratree->Branch("pholead_energy",&pholead_energy,"pholead_energy/F");
    thistree->Branch("photrail_energy",&photrail_energy,"photrail_energy/F");
    if (thisextratree) thisextratree->Branch("photrail_energy",&photrail_energy,"photrail_energy/F");
    thistree->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
    if (thisextratree) thisextratree->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
    thistree->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
    if (thisextratree) thisextratree->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
    thistree->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
    if (thisextratree) thisextratree->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
    thistree->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");
    if (thisextratree) thisextratree->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");
  
    thistree->Branch("pholead_PhoHasPixSeed",&pholead_PhoHasPixSeed,"pholead_PhoHasPixSeed/I");
    thistree->Branch("photrail_PhoHasPixSeed",&photrail_PhoHasPixSeed,"photrail_PhoHasPixSeed/I");

    thistree->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
    thistree->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");
    thistree->Branch("pholead_sieie",&pholead_sieie,"pholead_sieie/F");
    thistree->Branch("photrail_sieie",&photrail_sieie,"photrail_sieie/F");
    thistree->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");
    thistree->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");

    thistree->Branch("pholead_PhoSCRemovalPFIsoCharged",&pholead_PhoSCRemovalPFIsoCharged,"pholead_PhoSCRemovalPFIsoCharged/F");
    thistree->Branch("photrail_PhoSCRemovalPFIsoCharged",&photrail_PhoSCRemovalPFIsoCharged,"photrail_PhoSCRemovalPFIsoCharged/F");
    thistree->Branch("pholead_PhoSCRemovalPFIsoNeutral",&pholead_PhoSCRemovalPFIsoNeutral,"pholead_PhoSCRemovalPFIsoNeutral/F");
    thistree->Branch("photrail_PhoSCRemovalPFIsoNeutral",&photrail_PhoSCRemovalPFIsoNeutral,"photrail_PhoSCRemovalPFIsoNeutral/F");
    thistree->Branch("pholead_PhoSCRemovalPFIsoPhoton",&pholead_PhoSCRemovalPFIsoPhoton,"pholead_PhoSCRemovalPFIsoPhoton/F");
    thistree->Branch("photrail_PhoSCRemovalPFIsoPhoton",&photrail_PhoSCRemovalPFIsoPhoton,"photrail_PhoSCRemovalPFIsoPhoton/F");
    thistree->Branch("pholead_PhoSCRemovalPFIsoCombined",&pholead_PhoSCRemovalPFIsoCombined,"pholead_PhoSCRemovalPFIsoCombined/F");
    thistree->Branch("photrail_PhoSCRemovalPFIsoCombined",&photrail_PhoSCRemovalPFIsoCombined,"photrail_PhoSCRemovalPFIsoCombined/F");

    thistree->Branch("pholead_PhoIso03Ecal",&pholead_PhoIso03Ecal,"pholead_PhoIso03Ecal/F");
    thistree->Branch("pholead_PhoIso03Hcal",&pholead_PhoIso03Hcal,"pholead_PhoIso03Hcal/F");
    thistree->Branch("pholead_PhoIso03TrkHollow",&pholead_PhoIso03TrkHollow,"pholead_PhoIso03TrkHollow/F");
    thistree->Branch("photrail_PhoIso03Ecal",&photrail_PhoIso03Ecal,"photrail_PhoIso03Ecal/F");
    thistree->Branch("photrail_PhoIso03Hcal",&photrail_PhoIso03Hcal,"photrail_PhoIso03Hcal/F");
    thistree->Branch("photrail_PhoIso03TrkHollow",&photrail_PhoIso03TrkHollow,"photrail_PhoIso03TrkHollow/F");

    thistree->Branch("pholead_PhoPassConversionVeto",&pholead_PhoPassConversionVeto,"pholead_PhoPassConversionVeto/I");
    thistree->Branch("photrail_PhoPassConversionVeto",&photrail_PhoPassConversionVeto,"photrail_PhoPassConversionVeto/I");

    thistree->Branch("pholead_GenPhotonIsoDR04",&pholead_GenPhotonIsoDR04,"pholead_GenPhotonIsoDR04/F");
    thistree->Branch("photrail_GenPhotonIsoDR04",&photrail_GenPhotonIsoDR04,"photrail_GenPhotonIsoDR04/F");

    thistree->Branch("pholead_PhoMCmatchexitcode",&pholead_PhoMCmatchexitcode,"pholead_PhoMCmatchexitcode/I");
    thistree->Branch("photrail_PhoMCmatchexitcode",&photrail_PhoMCmatchexitcode,"photrail_PhoMCmatchexitcode/I");

    //  thistree->Branch("pholead_scarea",&pholead_scarea,"pholead_scarea/F");
    //  thistree->Branch("pholead_scareaSF",&pholead_scareaSF,"pholead_scareaSF/F");
    //  thistree->Branch("photrail_scarea",&photrail_scarea,"photrail_scarea/F");
    //  thistree->Branch("photrail_scareaSF",&photrail_scareaSF,"photrail_scareaSF/F");

    thistree->Branch("pholead_m_jet_ptcorr",&pholead_m_jet_ptcorr,"pholead_m_jet_ptcorr/F");
    thistree->Branch("pholead_m_jet_dR",&pholead_m_jet_dR,"pholead_m_jet_dR/F");
    thistree->Branch("pholead_phopt_footprint_total",&pholead_phopt_footprint_total,"pholead_phopt_footprint_total/F");
    thistree->Branch("pholead_phopt_footprint_m_frac",&pholead_phopt_footprint_m_frac,"pholead_phopt_footprint_m_frac/F");
    thistree->Branch("pholead_jetpt_pf",&pholead_jetpt_pf,"pholead_jetpt_pf/F");
    thistree->Branch("pholead_jetpt_m_frac",&pholead_jetpt_m_frac,"pholead_jetpt_m_frac/F");
    thistree->Branch("pholead_jetpt_m_frac_PhoComp",&pholead_jetpt_m_frac_PhoComp,"pholead_jetpt_m_frac_PhoComp/F");
    thistree->Branch("photrail_m_jet_ptcorr",&photrail_m_jet_ptcorr,"photrail_m_jet_ptcorr/F");
    thistree->Branch("photrail_m_jet_dR",&photrail_m_jet_dR,"photrail_m_jet_dR/F");
    thistree->Branch("photrail_phopt_footprint_total",&photrail_phopt_footprint_total,"photrail_phopt_footprint_total/F");
    thistree->Branch("photrail_phopt_footprint_m_frac",&photrail_phopt_footprint_m_frac,"photrail_phopt_footprint_m_frac/F");
    thistree->Branch("photrail_jetpt_pf",&photrail_jetpt_pf,"photrail_jetpt_pf/F");
    thistree->Branch("photrail_jetpt_m_frac",&photrail_jetpt_m_frac,"photrail_jetpt_m_frac/F");
    thistree->Branch("photrail_jetpt_m_frac_PhoComp",&photrail_jetpt_m_frac_PhoComp,"photrail_jetpt_m_frac_PhoComp/F");

    thistree->Branch("pholead_pt_closestjet",&pholead_pt_closestjet,"pholead_pt_closestjet/F");
    thistree->Branch("pholead_dR_closestjet",&pholead_dR_closestjet,"pholead_dR_closestjet/F");
    thistree->Branch("photrail_pt_closestjet",&photrail_pt_closestjet,"photrail_pt_closestjet/F");
    thistree->Branch("photrail_dR_closestjet",&photrail_dR_closestjet,"photrail_dR_closestjet/F");

    if (do_recalc_isolation){
      thistree->Branch("pholead_Npfcandphotonincone",&pholead_Npfcandphotonincone,"pholead_Npfcandphotonincone/I");
      thistree->Branch("pholead_Npfcandchargedincone",&pholead_Npfcandchargedincone,"pholead_Npfcandchargedincone/I");
      thistree->Branch("pholead_Npfcandneutralincone",&pholead_Npfcandneutralincone,"pholead_Npfcandneutralincone/I");
      thistree->Branch("photrail_Npfcandphotonincone",&photrail_Npfcandphotonincone,"photrail_Npfcandphotonincone/I");
      thistree->Branch("photrail_Npfcandchargedincone",&photrail_Npfcandchargedincone,"photrail_Npfcandchargedincone/I");
      thistree->Branch("photrail_Npfcandneutralincone",&photrail_Npfcandneutralincone,"photrail_Npfcandneutralincone/I");
      thistree->Branch("pholead_photonpfcandenergies",&pholead_photonpfcandenergies,"pholead_photonpfcandenergies[pholead_Npfcandphotonincone]/F");
      thistree->Branch("pholead_chargedpfcandenergies",&pholead_chargedpfcandenergies,"pholead_chargedpfcandenergies[pholead_Npfcandchargedincone]/F");
      thistree->Branch("pholead_neutralpfcandenergies",&pholead_neutralpfcandenergies,"pholead_neutralpfcandenergies[pholead_Npfcandneutralincone]/F");
      thistree->Branch("pholead_photonpfcandets",&pholead_photonpfcandets,"pholead_photonpfcandets[pholead_Npfcandphotonincone]/F");
      thistree->Branch("pholead_chargedpfcandets",&pholead_chargedpfcandets,"pholead_chargedpfcandets[pholead_Npfcandchargedincone]/F");
      thistree->Branch("pholead_neutralpfcandets",&pholead_neutralpfcandets,"pholead_neutralpfcandets[pholead_Npfcandneutralincone]/F");
      thistree->Branch("pholead_photonpfcanddetas",&pholead_photonpfcanddetas,"pholead_photonpfcanddetas[pholead_Npfcandphotonincone]/F");
      thistree->Branch("pholead_chargedpfcanddetas",&pholead_chargedpfcanddetas,"pholead_chargedpfcanddetas[pholead_Npfcandchargedincone]/F");
      thistree->Branch("pholead_neutralpfcanddetas",&pholead_neutralpfcanddetas,"pholead_neutralpfcanddetas[pholead_Npfcandneutralincone]/F");
      thistree->Branch("pholead_photonpfcanddphis",&pholead_photonpfcanddphis,"pholead_photonpfcanddphis[pholead_Npfcandphotonincone]/F");
      thistree->Branch("pholead_chargedpfcanddphis",&pholead_chargedpfcanddphis,"pholead_chargedpfcanddphis[pholead_Npfcandchargedincone]/F");
      thistree->Branch("pholead_neutralpfcanddphis",&pholead_neutralpfcanddphis,"pholead_neutralpfcanddphis[pholead_Npfcandneutralincone]/F");
      thistree->Branch("photrail_photonpfcandenergies",&photrail_photonpfcandenergies,"photrail_photonpfcandenergies[photrail_Npfcandphotonincone]/F");
      thistree->Branch("photrail_chargedpfcandenergies",&photrail_chargedpfcandenergies,"photrail_chargedpfcandenergies[photrail_Npfcandchargedincone]/F");
      thistree->Branch("photrail_neutralpfcandenergies",&photrail_neutralpfcandenergies,"photrail_neutralpfcandenergies[photrail_Npfcandneutralincone]/F");
      thistree->Branch("photrail_photonpfcandets",&photrail_photonpfcandets,"photrail_photonpfcandets[photrail_Npfcandphotonincone]/F");
      thistree->Branch("photrail_chargedpfcandets",&photrail_chargedpfcandets,"photrail_chargedpfcandets[photrail_Npfcandchargedincone]/F");
      thistree->Branch("photrail_neutralpfcandets",&photrail_neutralpfcandets,"photrail_neutralpfcandets[photrail_Npfcandneutralincone]/F");
      thistree->Branch("photrail_photonpfcanddetas",&photrail_photonpfcanddetas,"photrail_photonpfcanddetas[photrail_Npfcandphotonincone]/F");
      thistree->Branch("photrail_chargedpfcanddetas",&photrail_chargedpfcanddetas,"photrail_chargedpfcanddetas[photrail_Npfcandchargedincone]/F");
      thistree->Branch("photrail_neutralpfcanddetas",&photrail_neutralpfcanddetas,"photrail_neutralpfcanddetas[photrail_Npfcandneutralincone]/F");
      thistree->Branch("photrail_photonpfcanddphis",&photrail_photonpfcanddphis,"photrail_photonpfcanddphis[photrail_Npfcandphotonincone]/F");
      thistree->Branch("photrail_chargedpfcanddphis",&photrail_chargedpfcanddphis,"photrail_chargedpfcanddphis[photrail_Npfcandchargedincone]/F");
      thistree->Branch("photrail_neutralpfcanddphis",&photrail_neutralpfcanddphis,"photrail_neutralpfcanddphis[photrail_Npfcandneutralincone]/F");
    }

    //  // IF SCAN ROTATED CONES
    //  thistree->Branch("pholead_test_rotatedphotoniso",&pholead_test_rotatedphotoniso,"pholead_test_rotatedphotoniso[50]/F");
    //  thistree->Branch("pholead_test_rotatedwithcheckphotoniso",&pholead_test_rotatedwithcheckphotoniso,"pholead_test_rotatedwithcheckphotoniso[50]/F");

    if (thisextratree) thisextratree->Branch("allphotonpfcand_count",&allphotonpfcand_count,"allphotonpfcand_count/I");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_pt",&allphotonpfcand_pt,"allphotonpfcand_pt[allphotonpfcand_count]/F");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_eta",&allphotonpfcand_eta,"allphotonpfcand_eta[allphotonpfcand_count]/F");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_phi",&allphotonpfcand_phi,"allphotonpfcand_phi[allphotonpfcand_count]/F");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_vx",&allphotonpfcand_vx,"allphotonpfcand_vx[allphotonpfcand_count]/F");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_vy",&allphotonpfcand_vy,"allphotonpfcand_vy[allphotonpfcand_count]/F");
    if (thisextratree) thisextratree->Branch("allphotonpfcand_vz",&allphotonpfcand_vz,"allphotonpfcand_vz[allphotonpfcand_count]/F");

    thistree->Branch("phoiso_template_1event_sigsig_1",&phoiso_template_1event_sigsig_1,Form("phoiso_template_1event_sigsig_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_sigsig_2",&phoiso_template_1event_sigsig_2,Form("phoiso_template_1event_sigsig_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_sigbkg_1",&phoiso_template_1event_sigbkg_1,Form("phoiso_template_1event_sigbkg_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_sigbkg_2",&phoiso_template_1event_sigbkg_2,Form("phoiso_template_1event_sigbkg_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_bkgsig_1",&phoiso_template_1event_bkgsig_1,Form("phoiso_template_1event_bkgsig_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_bkgsig_2",&phoiso_template_1event_bkgsig_2,Form("phoiso_template_1event_bkgsig_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_bkgbkg_1",&phoiso_template_1event_bkgbkg_1,Form("phoiso_template_1event_bkgbkg_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_1event_bkgbkg_2",&phoiso_template_1event_bkgbkg_2,Form("phoiso_template_1event_bkgbkg_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_sigsig_1",&phoiso_template_2events_sigsig_1,Form("phoiso_template_2events_sigsig_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_sigsig_2",&phoiso_template_2events_sigsig_2,Form("phoiso_template_2events_sigsig_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_sigbkg_1",&phoiso_template_2events_sigbkg_1,Form("phoiso_template_2events_sigbkg_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_sigbkg_2",&phoiso_template_2events_sigbkg_2,Form("phoiso_template_2events_sigbkg_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_bkgsig_1",&phoiso_template_2events_bkgsig_1,Form("phoiso_template_2events_bkgsig_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_bkgsig_2",&phoiso_template_2events_bkgsig_2,Form("phoiso_template_2events_bkgsig_2[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_bkgbkg_1",&phoiso_template_2events_bkgbkg_1,Form("phoiso_template_2events_bkgbkg_1[%d]/F",nclosest));
    thistree->Branch("phoiso_template_2events_bkgbkg_2",&phoiso_template_2events_bkgbkg_2,Form("phoiso_template_2events_bkgbkg_2[%d]/F",nclosest));

    // rewinfo = {eta1, eta2, pt1, pt2, rho, sigma}
    thistree->Branch("rewinfo_template_1event_sigsig_1",&rewinfo_template_1event_sigsig_1,Form("rewinfo_template_1event_sigsig_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_sigsig_2",&rewinfo_template_1event_sigsig_2,Form("rewinfo_template_1event_sigsig_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_sigbkg_1",&rewinfo_template_1event_sigbkg_1,Form("rewinfo_template_1event_sigbkg_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_sigbkg_2",&rewinfo_template_1event_sigbkg_2,Form("rewinfo_template_1event_sigbkg_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_bkgsig_1",&rewinfo_template_1event_bkgsig_1,Form("rewinfo_template_1event_bkgsig_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_bkgsig_2",&rewinfo_template_1event_bkgsig_2,Form("rewinfo_template_1event_bkgsig_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_bkgbkg_1",&rewinfo_template_1event_bkgbkg_1,Form("rewinfo_template_1event_bkgbkg_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_1event_bkgbkg_2",&rewinfo_template_1event_bkgbkg_2,Form("rewinfo_template_1event_bkgbkg_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_sigsig_1",&rewinfo_template_2events_sigsig_1,Form("rewinfo_template_2events_sigsig_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_sigsig_2",&rewinfo_template_2events_sigsig_2,Form("rewinfo_template_2events_sigsig_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_sigbkg_1",&rewinfo_template_2events_sigbkg_1,Form("rewinfo_template_2events_sigbkg_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_sigbkg_2",&rewinfo_template_2events_sigbkg_2,Form("rewinfo_template_2events_sigbkg_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_bkgsig_1",&rewinfo_template_2events_bkgsig_1,Form("rewinfo_template_2events_bkgsig_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_bkgsig_2",&rewinfo_template_2events_bkgsig_2,Form("rewinfo_template_2events_bkgsig_2[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_bkgbkg_1",&rewinfo_template_2events_bkgbkg_1,Form("rewinfo_template_2events_bkgbkg_1[%d]/F",nclosest*6));
    thistree->Branch("rewinfo_template_2events_bkgbkg_2",&rewinfo_template_2events_bkgbkg_2,Form("rewinfo_template_2events_bkgbkg_2[%d]/F",nclosest*6));

    if (thisextratree) thisextratree->Branch("vetoobjects_count",&vetoobjects_count,"vetoobjects_count/I");
    if (thisextratree) thisextratree->Branch("vetoobjects_pt",&vetoobjects_pt,"vetoobjects_pt[vetoobjects_count]/F");
    if (thisextratree) thisextratree->Branch("vetoobjects_eta",&vetoobjects_eta,"vetoobjects_eta[vetoobjects_count]/F");
    if (thisextratree) thisextratree->Branch("vetoobjects_phi",&vetoobjects_phi,"vetoobjects_phi[vetoobjects_count]/F");
    if (thisextratree) thisextratree->Branch("vetoobjects_type",&vetoobjects_type,"vetoobjects_type[vetoobjects_count]/I");

    thistree->Branch("n_jets",&n_jets,"n_jets/I");
    thistree->Branch("jet_pt",&jet_pt,"jet_pt[n_jets]/F");
    thistree->Branch("jet_eta",&jet_eta,"jet_eta[n_jets]/F");
    thistree->Branch("jet_phi",&jet_phi,"jet_phi[n_jets]/F");
    thistree->Branch("jet_energy",&jet_energy,"jet_energy[n_jets]/F");

  }

  inputtree_isinitialized = false;

  for (map<EnumSel,Selection>::iterator it = sels.begin(); it!=sels.end(); it++){

    if (!it->second.islightsel) continue;

    TTree *thislighttree = it->second.tree_light;

    thislighttree->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
    thislighttree->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");
    thislighttree->Branch("event_weight",&event_weight,"event_weight/F");
    thislighttree->Branch("dataset_id",&mydataset_id,"dataset_id/I");
    thislighttree->Branch("event_nPUtrue",&event_nPUtrue,"event_nPUtrue/I");

    thislighttree->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
    thislighttree->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
    thislighttree->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
    thislighttree->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
    thislighttree->Branch("pholead_phi",&pholead_phi,"pholead_phi/F");
    thislighttree->Branch("photrail_phi",&photrail_phi,"photrail_phi/F");
    thislighttree->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
    thislighttree->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
    thislighttree->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
    thislighttree->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");

    thislighttree->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
    thislighttree->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");
    thislighttree->Branch("dipho_mgg_photon",&dipho_mgg_photon,"dipho_mgg_photon/F");

    thislighttree->Branch("pholead_GEN_pt",&pholead_GEN_pt,"pholead_GEN_pt/F");
    thislighttree->Branch("photrail_GEN_pt",&photrail_GEN_pt,"photrail_GEN_pt/F");
    thislighttree->Branch("pholead_GEN_eta",&pholead_GEN_eta,"pholead_GEN_eta/F");
    thislighttree->Branch("photrail_GEN_eta",&photrail_GEN_eta,"photrail_GEN_eta/F");
    thislighttree->Branch("pholead_GEN_phi",&pholead_GEN_phi,"pholead_GEN_phi/F");
    thislighttree->Branch("photrail_GEN_phi",&photrail_GEN_phi,"photrail_GEN_phi/F");

    thislighttree->Branch("n_jets",&n_jets,"n_jets/I");
    thislighttree->Branch("jet_pt",&jet_pt,"jet_pt[n_jets]/F");
    thislighttree->Branch("jet_eta",&jet_eta,"jet_eta[n_jets]/F");
    thislighttree->Branch("jet_phi",&jet_phi,"jet_phi[n_jets]/F");
    thislighttree->Branch("jet_energy",&jet_energy,"jet_energy[n_jets]/F");

    thislighttree->Branch("n_GEN_jets",&n_GEN_jets,"n_GEN_jets/I");
    thislighttree->Branch("jet_GEN_pt",&jet_GEN_pt,"jet_GEN_pt[n_GEN_jets]/F");
    thislighttree->Branch("jet_GEN_eta",&jet_GEN_eta,"jet_GEN_eta[n_GEN_jets]/F");
    thislighttree->Branch("jet_GEN_phi",&jet_GEN_phi,"jet_GEN_phi[n_GEN_jets]/F");
    thislighttree->Branch("jet_GEN_energy",&jet_GEN_energy,"jet_GEN_energy[n_GEN_jets]/F");

    thislighttree->Branch("gen_in_acc",&tree_gen_in_acc,"gen_in_acc/O");
    thislighttree->Branch("reco_in_acc",&tree_reco_in_acc,"reco_in_acc/O");
    thislighttree->Branch("matched",&tree_matched,"matched/O");

  }


  fHNumPU = new TH1F("NumPU_rew","NumPU_rew",100,0,100);
  fHNumPU_noweight = new TH1F("NumPU_noweight","NumPU_noweight",100,0,100);
  fHNumPUTrue = new TH1F("NumPUTrue_rew","NumPUTrue_rew",100,0,100);
  fHNumPUTrue_noweight = new TH1F("NumPUTrue_noweight","NumPUTrue_noweight",100,0,100);
  fHNumVtx = new TH1F("NumVtx_rew","NumVtx_rew",100,0,100);
  fPNumEvents = new TParameter<double>("NumEvents",0);
  fPSumWeights = new TParameter<double>("SumGenWeights",0);
  fPSumWeights2 = new TParameter<double>("SumGenWeights2",0);

  cout << "Trees and histos created" << endl;

  InitEnergyScalesAndSmearingsDatabase();

  std::string mygtag="";
  if (year==2011) mygtag = (isdata) ? "FT_R_53_LV3" : "START53_LV4";
  else if (year==2012) mygtag="TODO";
  std::string mypath = "/shome/peruzzi/JetCorrectionFiles/";
  MyJetCorrector = new OnTheFlyCorrections(mygtag,isdata,mypath);
  cout << "Using GTAG " << mygtag << " for JEC/JER" << endl;

}

void DiPhotonMiniTree::Analyze(){

  //  cout << endl << "-----------------------" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumTrueInteractions); // TRUE
  //  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions); // OBSERVED
  else weight=1;
  
  event_luminormfactor=AddWeight;
  if (!isdata) event_luminormfactor*=fTR->GenWeight;

  fPNumEvents->SetVal(fPNumEvents->GetVal()+1);
  fPSumWeights->SetVal(fPSumWeights->GetVal()+fTR->GenWeight);
  fPSumWeights2->SetVal(fPSumWeights2->GetVal()+fTR->GenWeight*fTR->GenWeight);

  event_fileuuid = uuid;

  event_run = fTR->Run;
  event_lumi = fTR->LumiSection;
  event_number = fTR->Event;
  mydataset_id = dataset_id;

  if (!isdata) {
    fHNumPU->Fill(fTR->PUnumInteractions,weight);
    fHNumPU_noweight->Fill(fTR->PUnumInteractions);
    fHNumPUTrue->Fill(fTR->PUnumTrueInteractions,weight);
    fHNumPUTrue_noweight->Fill(fTR->PUnumTrueInteractions);
  }
  fHNumVtx->Fill(fTR->NVrtx,weight);

  //  return; // RUNNING FOR PU FILE


  event_weight = weight;
  event_rho = fTR->Rho;
  event_sigma = fTR->Sigma;
  if (!isdata) {
    event_nPU = fTR->PUnumInteractions;
    event_nPUtrue = fTR->PUnumTrueInteractions;
  }
  event_nRecVtx = fTR->NVrtx;

  if (!isdata) {
    int npart_isrfsr_gamma=Count_part_isrfsr_gamma(fTR); // TO FIX WITH PROCESS ID
    event_Kfactor = kfactors[2-npart_isrfsr_gamma];
  }
  else event_Kfactor=1;
  

//  // sc area and scale factor calculation for isolation // TO BE IMPROVED WITH EFFECTIVE CONE AREA
//  {
//    float const conearea = TMath::Pi()*0.4*0.4;
//    //    std::cout << "nsc " << fTR->NSuperClusters << std::endl;
//    for (int i=0; i<100; i++){
//      if (i<fTR->NSuperClusters){
//	scarea[i] = CalculateSCArea(fTR,i); 
//	scareaSF[i] = conearea/(conearea-scarea[i]);
//	if (scareaSF[i]<0) std::cout << "SC area larger than isolation cone!!!" << std::endl;
//      }
//      else {
//	scarea[i] = -999;
//	scareaSF[i] = -999;
//      }
//    }
//  }

  bool passtrigger = TriggerSelection();

  std::vector<std::pair<int,float> > ordering;
  std::vector<std::pair<int,float> > ordering_jets;
  std::vector<float> unscaled_energy;
  std::vector<float> unscaled_r9;
  std::vector<float> unscaled_sieie;

  StatusScaleUpScaleDown StatusScaleUpScaleDown_Photons_R9=kUncorrected;
  StatusScaleUpScaleDown StatusScaleUpScaleDown_Photons_Sieie=kUncorrected;
  StatusScaleUpScaleDown StatusScaleUpScaleDown_Photons_EnergyScale=kUncorrected;
  StatusScaleUpScaleDown StatusScaleUpScaleDown_Photons_EnergySmear=kUncorrected;
  StatusScaleUpScaleDown StatusScaleUpScaleDown_Jets_JEC=kUncorrected;
  StatusScaleUpScaleDown StatusScaleUpScaleDown_Jets_JER=kUncorrected;


  // rescale shower shapes and correct energy (once and for all, these functions should never be called twice on the same photon)
  for (int i=0; i<fTR->NPhotons; i++){
    if (!isdata) FixMatchingStatusElectrons(i); // correct match exit code for photons matched to gen electrons
    if (dataset_id==sherpa_dataset_id) FixMatchingStatusSherpa(i); // correct match exit code for SHERPA
    unscaled_energy.push_back(fTR->PhoEnergy[i]);
  }

  if (debug) {for (int i=0; i<fTR->NPhotons; i++) cout << fTR->PhoPt[i] << " "; cout<<endl;}

  StatusScaleUpScaleDown_Photons_R9 = CorrectPhotonR9(kCorrectedNoShift,unscaled_r9);
  StatusScaleUpScaleDown_Photons_Sieie = CorrectPhotonSieie(kCorrectedNoShift,unscaled_sieie);
  CorrectPhotonEnergy(kCorrectedNoShift,kCorrectedNoShift,ordering,unscaled_r9,unscaled_energy);
  StatusScaleUpScaleDown_Photons_EnergyScale=kCorrectedNoShift;
  StatusScaleUpScaleDown_Photons_EnergySmear=kCorrectedNoShift;

  if (debug) {for (int i=0; i<fTR->NPhotons; i++) cout << fTR->PhoPt[i] << " "; cout<<endl;}

  JECJERCorrection(kCorrectedNoShift,kCorrectedNoShift,ordering_jets);
  StatusScaleUpScaleDown_Jets_JEC=kCorrectedNoShift;
  StatusScaleUpScaleDown_Jets_JER=kCorrectedNoShift;

  
  photon_jet_matching.clear();
  for (int i=0; i<fTR->NPhotons; i++) photon_jet_matching.push_back(PFMatchPhotonToJet(i));


  map<EnumSel,vector<int> > passing_selection;
  map<EnumSel,vector<int> > passing_selection_jets;
  map<EnumSel,bool> pass;
  map<EnumSel,int> pass12_whoissiglike;
  map<EnumSel,bool> donethisfullsel;

  if (debug) cout << "start selection loop" << endl;

  for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++){

    const EnumSel sel_cat = sel->second.id;
    donethisfullsel[sel_cat]=false;

    if (!(sel->second.isfullsel)) continue;
    if (do_only_light_tree) continue;
    if (isstep2 && sel_cat!=0) continue;
    if (isdata && !(sel->second.runondata)) continue;
    if (!isdata && !(sel->second.runonmc)) continue;

    if (sel_cat!=k2DZmumu_selection && !passtrigger) continue; // no trigger for Zmumu selection

    pass[sel_cat]=false;
    donethisfullsel[sel_cat]=true;

    std::vector<int> passing;
    std::vector<int> passing_jets;

    if (sel_cat==k2DZmumu_selection){
      for (int i=0; i<fTR->NMus; i++){
        passing.push_back(i);
      }
      passing = MuonSelection(fTR,passing);
    }
    else {
      for (int i=0; i<fTR->NPhotons; i++){
	passing.push_back(ordering.at(i).first);
      }
      for (int i=0; i<(int)(passing.size())-1; i++){
	assert(fTR->PhoPt[passing.at(i)]>=fTR->PhoPt[passing.at(i+1)]);
      }
      passing = PhotonPreSelection(fTR,passing);

      for (int i=0; i<(int)(ordering_jets.size()); i++) passing_jets.push_back(ordering_jets.at(i).first);
      for (int i=0; i<(int)(passing_jets.size())-1; i++){
	if (!(fTR->JPt[passing_jets.at(i)]>=fTR->JPt[passing_jets.at(i+1)])) cout << event_run << ":" << event_lumi << ":" << event_number << " " << fTR->JPt[passing_jets.at(i)] << " " << fTR->JPt[passing_jets.at(i+1)] << endl;
        assert(fTR->JPt[passing_jets.at(i)]>=fTR->JPt[passing_jets.at(i+1)]);
      }
      JetSelection(passing_jets);

    }


    if (sel_cat==k2Dstandard_selection){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k1Drandomcone_template){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k1Dsideband_template){
      //      std::vector<int> passing_sel = PhotonSelection(fTR,passing);
      //      std::vector<int> passing_presel = passing;
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      //      pass[sel_cat] = (passing_sel.size()>=1) && SinglePhotonEventSelection(fTR,passing);
      //      pass[sel_cat] = (passing_presel.size()>=2) && SinglePhotonEventSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k2DZee_pixelvetoreversed_selection){
      passing = PhotonSelection(fTR,passing,"revert_pixel_veto"); // revert pixel veto already done
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k1Dpreselection){
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k1Dselection){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k2Drandomcone_template){
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Drandomconesideband_template){
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> passing_sel = PhotonSelection(fTR,passing);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sel.size()>=1){
	int fondo = passing_bkg[0];
	int forcone = passing_sel[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike[sel_cat]=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike[sel_cat]=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Dsideband_template){
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Dstandard_preselection){
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2DZmumu_selection){
      pass[sel_cat] = DiMuonFromZSelection(fTR,passing);
    }
    else if (sel_cat==k1Dsignal_template){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k1Dbackground_template){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==k2Dtruesigsig_template){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Dtruesigbkg_template){
      passing = PhotonSelection(fTR,passing);
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike[sel_cat]=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike[sel_cat]=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Dtruebkgbkg_template){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Drconeplusgenfake_template){
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      passing_bkg = PhotonSelection(fTR,passing_bkg);
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      passing_sig = PhotonSelection(fTR,passing_sig);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int forcone = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike[sel_cat]=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike[sel_cat]=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);
    }
    else if (sel_cat==k2Dgenpromptplussideband_template){
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      passing_sig = PhotonSelection(fTR,passing_sig);
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike[sel_cat]=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike[sel_cat]=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing,passing_jets);      
    }

    passing_selection[sel_cat] = passing;
    passing_selection_jets[sel_cat] = passing_jets;

  }

  if (debug) cout << "start filling loop" << endl;

  for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++){

    const EnumSel sel_cat = sel->second.id;
    if (!donethisfullsel[sel_cat]) continue;
    if (!pass[sel_cat]) continue;

    if (debug) cout << sel_cat << endl;

    std::vector<int> passing = passing_selection[sel_cat];
    std::vector<int> passing_jets = passing_selection_jets[sel_cat];
    int minsize = (sel->second.is2d) ? 2 : 1;

    if (passing_selection[sel_cat].size()<minsize){
      std::cout << "Error!!!" << std::endl;
      continue;
    }

    if (sel_cat==k2DZmumu_selection){
      for (int i=0; i<passing.size(); i++){
	ResetVars();
	FillMuonInfo(passing.at(i));
	sel->second.tree_full->Fill();
      }
    }
    else if (sel->second.is2d){
      ResetVars();
      FillLead(passing.at(0),passing_jets);
      FillTrail(passing.at(1),passing_jets);
      FillJetsInfo(passing,passing_jets);
      bool dofill=true;

      if (sel_cat==k2Drandomconesideband_template || sel_cat==k2Dtruesigbkg_template || sel_cat==k2Drconeplusgenfake_template || sel_cat==k2Dgenpromptplussideband_template) event_pass12whoissiglike=pass12_whoissiglike[sel_cat];

      if (sel_cat==k2Drandomcone_template || ((sel_cat==k2Drandomconesideband_template || sel_cat==k2Drconeplusgenfake_template) && pass12_whoissiglike[sel_cat]==0)) {

	  if (do_recalc_isolation){
	    isolations_struct rcone_isos;
	    rcone_isos = RandomConeIsolation(fTR,passing.at(0),"");
	    pholead_phi = rcone_isos.newphi;
	    pholead_SCphi = rcone_isos.newphi;
	    pholead_PhoSCRemovalPFIsoPhoton = rcone_isos.photon;
	    pholead_PhoSCRemovalPFIsoNeutral = rcone_isos.neutral;
	    pholead_PhoSCRemovalPFIsoCharged = rcone_isos.charged;
	    pholead_PhoSCRemovalPFIsoCombined = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	    pholead_Npfcandphotonincone = rcone_isos.nphotoncand;
	    pholead_Npfcandchargedincone = rcone_isos.nchargedcand;
	    pholead_Npfcandneutralincone = rcone_isos.nneutralcand;
	    for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	    for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	    for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	    for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	    for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	    for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	    for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	    for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	    for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	    for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	    for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	    for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	    if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
	  }
	  else{
	    pholead_phi = fTR->PhoSCRemovalRConePhi[passing.at(0)];
	    pholead_SCphi = fTR->PhoSCRemovalRConePhi[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedPrimVtxRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoCombined = pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
	    if (pholead_PhoSCRemovalPFIsoCharged==999 || pholead_PhoSCRemovalPFIsoNeutral==999 || pholead_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	  }
      }
      if (sel_cat==k2Drandomcone_template || ((sel_cat==k2Drandomconesideband_template || sel_cat==k2Drconeplusgenfake_template) && pass12_whoissiglike[sel_cat]==1)) {

	if (do_recalc_isolation){
	  isolations_struct rcone_isos;
	  rcone_isos = RandomConeIsolation(fTR,passing.at(1),"");
	  photrail_phi= rcone_isos.newphi;
	  photrail_SCphi= rcone_isos.newphi;
	  photrail_PhoSCRemovalPFIsoPhoton = rcone_isos.photon;
	  photrail_PhoSCRemovalPFIsoNeutral = rcone_isos.neutral;
	  photrail_PhoSCRemovalPFIsoCharged = rcone_isos.charged;
	  photrail_PhoSCRemovalPFIsoCombined = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	  photrail_Npfcandphotonincone = rcone_isos.nphotoncand;
	  photrail_Npfcandchargedincone = rcone_isos.nchargedcand;
	  photrail_Npfcandneutralincone = rcone_isos.nneutralcand;
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	  if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
	}
	else{
	  photrail_phi = fTR->PhoSCRemovalRConePhi[passing.at(1)];
	  photrail_SCphi = fTR->PhoSCRemovalRConePhi[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedPrimVtxRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoCombined = photrail_PhoSCRemovalPFIsoCharged+photrail_PhoSCRemovalPFIsoNeutral+photrail_PhoSCRemovalPFIsoPhoton;
	  if (photrail_PhoSCRemovalPFIsoCharged==999 || photrail_PhoSCRemovalPFIsoNeutral==999 || photrail_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	}

	}


      if (sel_cat==k2Dstandard_selection && isstep2){ // new templates from event mixing

	if (!inputtree_isinitialized) InitInputTree();

	Long64_t index_matchingtree = matchingtree->GetEntryNumberWithIndex(event_number%2,event_number>>1);
	if (index_matchingtree<0) {
	  cout << "NO MATCHING FOUND (including under/overflow in templ. variable)" << endl; cout << event_run << " " << event_lumi << " " << event_number << endl; 
	  dofill=false;
	}
	else{
	  matchingtree->GetEntry(index_matchingtree);
	  if (event_run!=matchingtree_event_run || event_lumi!=matchingtree_event_lumi || event_number!=matchingtree_event_number){
	    cout << "WRONG MATCHING" << endl;
	    cerr << "ERROR IN MATCHING" << endl;
	    dofill=false;
	  }
	}

	if (dofill){
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_1event_sigsig_1,matchingtree_index_1event_sigsig_2,passing,kSigSig,k1Event);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_1event_sigbkg_1,matchingtree_index_1event_sigbkg_2,passing,kSigBkg,k1Event);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_1event_bkgsig_1,matchingtree_index_1event_bkgsig_2,passing,kBkgSig,k1Event);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_1event_bkgbkg_1,matchingtree_index_1event_bkgbkg_2,passing,kBkgBkg,k1Event);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_2events_sigsig_1,matchingtree_index_2events_sigsig_2,passing,kSigSig,k2Events);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_2events_sigbkg_1,matchingtree_index_2events_sigbkg_2,passing,kSigBkg,k2Events);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_2events_bkgsig_1,matchingtree_index_2events_bkgsig_2,passing,kBkgSig,k2Events);
	  FillPhoIso_NewTemplates(fTR,matchingtree_index_2events_bkgbkg_1,matchingtree_index_2events_bkgbkg_2,passing,kBkgBkg,k2Events);
	}

//	for (int l=0; l<nclosest; l++) cout << matchingtree_index_sigsig_1[l] << endl;
//	for (int l=0; l<nclosest; l++) cout << phoiso_template_sigsig_1[l] << endl;


      }

      if (sel_cat==k2Drandomconesideband_template) {
	std::set<int> removals = (do_recalc_isolation) ? GetPFCandInsideFootprint(fTR,passing.at(!pass12_whoissiglike[sel_cat]),0,"photon") : GetPrecalculatedFootprintPhoEl(passing.at(!pass12_whoissiglike[sel_cat]));
	int index=0;
	for (int k=0; k<fTR->NPfCand; k++){
	  if (index==global_maxN_photonpfcandidates) {std::cout << "Too many pfcandidates" << std::endl; dofill=false; break;}
	  if (fTR->PfCandPdgId[k]!=22) continue;
	  float eta = fabs(fTR->PfCandEta[k]);
	  if (eta>1.4442 && eta<1.566) continue;
	  if (eta>2.5) continue;
	  if (fTR->PhoMatchedPFPhotonCand.size()>0) if (fTR->PhoMatchedPFPhotonCand[passing.at(!pass12_whoissiglike[sel_cat])]==k) continue;	
	  if (fTR->PhoMatchedPFElectronCand.size()>0) if (fTR->PhoMatchedPFElectronCand[passing.at(!pass12_whoissiglike[sel_cat])]==k) continue;	
	  if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) if (fTR->PhoMatchedPFPhotonOrElectronCand[passing.at(!pass12_whoissiglike[sel_cat])]==k) continue;	
	  bool removed = false;
	  for (set<int>::iterator j=removals.begin(); j!=removals.end(); j++) if (k==*j) removed=true;
	  if (removed) continue;
	  allphotonpfcand_pt[index] = fTR->PfCandPt[k];
	  allphotonpfcand_eta[index] = fTR->PfCandEta[k];
	  allphotonpfcand_phi[index] = fTR->PfCandPhi[k];
	  allphotonpfcand_vx[index] = fTR->PfCandVx[k];
	  allphotonpfcand_vy[index] = fTR->PfCandVy[k];
	  allphotonpfcand_vz[index] = fTR->PfCandVz[k];
	  index++;
	}
	allphotonpfcand_count = index;

	FillVetoObjects(fTR,passing.at(!pass12_whoissiglike[sel_cat]),TString("exclude_object_itself"));

      }


      TLorentzVector pho1(fTR->PhoPx[passing.at(0)],fTR->PhoPy[passing.at(0)],fTR->PhoPz[passing.at(0)],fTR->PhoEnergy[passing.at(0)]);
      TLorentzVector pho2(fTR->PhoPx[passing.at(1)],fTR->PhoPy[passing.at(1)],fTR->PhoPz[passing.at(1)],fTR->PhoEnergy[passing.at(1)]);
      float invmass = (pho1+pho2).M();
      dipho_mgg_photon = invmass;
      if (dofill) sel->second.tree_full->Fill();
      if (dofill && sel_cat==k2Drandomconesideband_template) if (isdata && !isstep2) sel->second.tree_extra->Fill();
    }

    else if (!(sel->second.is2d)){

      for (int i=0; i<passing.size(); i++){
	if (debug) cout << "pho" << endl;
      ResetVars();
      FillLead(passing.at(i),passing_jets);
      bool dofill = true;

      if (sel_cat==k1Drandomcone_template) {

	if (do_recalc_isolation){
	isolations_struct rcone_isos;
	rcone_isos = RandomConeIsolation(fTR,passing.at(i),"");
	pholead_phi = rcone_isos.newphi;
	pholead_SCphi = rcone_isos.newphi;
	pholead_PhoSCRemovalPFIsoPhoton = rcone_isos.photon;
	pholead_PhoSCRemovalPFIsoNeutral = rcone_isos.neutral;
	pholead_PhoSCRemovalPFIsoCharged = rcone_isos.charged;
	pholead_PhoSCRemovalPFIsoCombined = rcone_isos.photon+rcone_isos.neutral+rcone_isos.charged;
	pholead_Npfcandphotonincone = rcone_isos.nphotoncand;
	pholead_Npfcandchargedincone = rcone_isos.nchargedcand;
	pholead_Npfcandneutralincone = rcone_isos.nneutralcand;
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = rcone_isos.photoncandenergies.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = rcone_isos.chargedcandenergies.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = rcone_isos.neutralcandenergies.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = rcone_isos.photoncandets.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = rcone_isos.chargedcandets.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = rcone_isos.neutralcandets.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = rcone_isos.photoncanddetas.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = rcone_isos.chargedcanddetas.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = rcone_isos.neutralcanddetas.at(i);
	for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = rcone_isos.photoncanddphis.at(i);
	for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = rcone_isos.chargedcanddphis.at(i);
	for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = rcone_isos.neutralcanddphis.at(i);
	if (rcone_isos.photon==-999 || rcone_isos.neutral==-999 || rcone_isos.charged==-999) dofill=false;
	}
	else{
	pholead_phi = fTR->PhoSCRemovalRConePhi[passing.at(i)];
	pholead_SCphi = fTR->PhoSCRemovalRConePhi[passing.at(i)];
	pholead_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedPrimVtxRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoCombined = pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
	if (pholead_PhoSCRemovalPFIsoCharged==999 || pholead_PhoSCRemovalPFIsoNeutral==999 || pholead_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	}

	if (debug) cout << "here" << endl;
      }
      
//      if (sel_cat==k1Dsignal_template || sel_cat==k1Dbackground_template) for (int k=0; k<50; k++) {
//	pholead_test_rotatedphotoniso[k]=PFIsolation(passing.at(i),0.025*k,"photon",NULL,NULL,NULL,NULL,NULL,NULL);
//	if (!FindCloseJetsAndPhotons(fTR,0.025*k,passing.at(i),"")) pholead_test_rotatedwithcheckphotoniso[k]=PFIsolation(passing.at(i),0.025*k,"photon",NULL,NULL,NULL,NULL,NULL,NULL);
//      }


      if (sel_cat==k1Drandomcone_template || sel_cat==k1Dsideband_template) {
	std::set<int> removals = (do_recalc_isolation) ? GetPFCandInsideFootprint(fTR,passing.at(i),0,"photon") : GetPrecalculatedFootprintPhoEl(passing.at(i));
	int index=0;
	for (int k=0; k<fTR->NPfCand; k++){
	  //	  if (index==global_maxN_photonpfcandidates) {std::cout << "Too many pfcandidates" << std::endl; dofill=false; break;}
	  if (fTR->PfCandPdgId[k]!=22) continue;
	  //	  if (debug) cout << fTR->PfCandPt[k] << " " << fTR->PfCandEta[k] << " " << fTR->PfCandPhi[k] << endl;
	  float eta = fabs(fTR->PfCandEta[k]);
	  if (eta>1.4442 && eta<1.566) continue;
	  if (eta>2.5) continue;
	  if (fTR->PhoMatchedPFPhotonCand.size()>0) if (fTR->PhoMatchedPFPhotonCand[passing.at(i)]==k) continue;	
	  if (fTR->PhoMatchedPFElectronCand.size()>0) if (fTR->PhoMatchedPFElectronCand[passing.at(i)]==k) continue;
	  if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) if (fTR->PhoMatchedPFPhotonOrElectronCand[passing.at(i)]==k) continue;
	  bool removed = false;
	  for (set<int>::iterator j=removals.begin(); j!=removals.end(); j++) if (k==*j) removed=true;
	  if (removed) continue;
	  allphotonpfcand_pt[index] = fTR->PfCandPt[k];
	  allphotonpfcand_eta[index] = fTR->PfCandEta[k];
	  allphotonpfcand_phi[index] = fTR->PfCandPhi[k];
	  allphotonpfcand_vx[index] = fTR->PfCandVx[k];
	  allphotonpfcand_vy[index] = fTR->PfCandVy[k];
	  allphotonpfcand_vz[index] = fTR->PfCandVz[k];
	  //	  cout << "taken" << endl;
	  index++;
	}
	if (debug) cout << "done loop" << endl;
	allphotonpfcand_count = index;

	FillVetoObjects(fTR,passing.at(i),(sel_cat==k1Drandomcone_template) ? TString("") : TString("exclude_object_itself"));

      }

      if (dofill) sel->second.tree_full->Fill();
      if (dofill && (sel_cat==k1Drandomcone_template || sel_cat==k1Dsideband_template)) if (isdata && !isstep2) sel->second.tree_extra->Fill();
      }

    }

  }

  if (debug) cout << "start light tree loop" << endl;
  
  for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++){ // lightweight tree for efficiency and unfolding studies

    if (!(sel->second.islightsel)) continue;
    if (isdata && !(sel->second.runondata)) continue;
    if (!isdata && !(sel->second.runonmc)) continue;

    EnumSel thissel = sel->second.id;
    if (debug) cout << thissel << endl;

    // revert in any case to default corrections
    CorrectPhotonEnergy(kCorrectedNoShift,kCorrectedNoShift,ordering,unscaled_r9,unscaled_energy);
    JECJERCorrection(kCorrectedNoShift,kCorrectedNoShift,ordering_jets);
    StatusScaleUpScaleDown_Photons_EnergyScale==kCorrectedNoShift;
    StatusScaleUpScaleDown_Photons_EnergySmear==kCorrectedNoShift;
    StatusScaleUpScaleDown_Jets_JEC=kCorrectedNoShift;
    StatusScaleUpScaleDown_Jets_JER=kCorrectedNoShift;

    if (thissel==kLightDefault || thissel==kLightDiElectron);
    else if (thissel==kLightJECup){
      JECJERCorrection(kShiftUp,kCorrectedNoShift,ordering_jets);
      StatusScaleUpScaleDown_Jets_JEC=kShiftUp;
      StatusScaleUpScaleDown_Jets_JER=kCorrectedNoShift;
    }
    else if (thissel==kLightJECdown){
      JECJERCorrection(kShiftDown,kCorrectedNoShift,ordering_jets);
      StatusScaleUpScaleDown_Jets_JEC=kShiftDown;
      StatusScaleUpScaleDown_Jets_JER=kCorrectedNoShift;
      }
    else if (thissel==kLightJERup){
      JECJERCorrection(kCorrectedNoShift,kShiftUp,ordering_jets);
      StatusScaleUpScaleDown_Jets_JEC=kCorrectedNoShift;
      StatusScaleUpScaleDown_Jets_JER=kShiftUp;
      }
    else if (thissel==kLightJERdown){
      JECJERCorrection(kCorrectedNoShift,kShiftDown,ordering_jets);
      StatusScaleUpScaleDown_Jets_JEC=kCorrectedNoShift;
      StatusScaleUpScaleDown_Jets_JER=kShiftDown;
      }
    else if (thissel==kLightESCALEup){
      CorrectPhotonEnergy(kShiftUp,kCorrectedNoShift,ordering,unscaled_r9,unscaled_energy);
      StatusScaleUpScaleDown_Photons_EnergyScale==kShiftUp;
      StatusScaleUpScaleDown_Photons_EnergySmear==kCorrectedNoShift;
      }
    else if (thissel==kLightESCALEdown){
      CorrectPhotonEnergy(kShiftDown,kCorrectedNoShift,ordering,unscaled_r9,unscaled_energy);
      StatusScaleUpScaleDown_Photons_EnergyScale==kShiftDown;
      StatusScaleUpScaleDown_Photons_EnergySmear==kCorrectedNoShift;
      }
    else if (thissel==kLightESMEARup){
      CorrectPhotonEnergy(kCorrectedNoShift,kShiftUp,ordering,unscaled_r9,unscaled_energy);
      StatusScaleUpScaleDown_Photons_EnergyScale==kCorrectedNoShift;
      StatusScaleUpScaleDown_Photons_EnergySmear==kShiftUp;
      }
    else if (thissel==kLightESMEARdown){
      CorrectPhotonEnergy(kCorrectedNoShift,kShiftDown,ordering,unscaled_r9,unscaled_energy);
      StatusScaleUpScaleDown_Photons_EnergyScale==kCorrectedNoShift;
      StatusScaleUpScaleDown_Photons_EnergySmear==kShiftDown;
      }
    else assert(false);


    ResetVars(); 

    std::vector<int> passing;
    for (int i=0; i<fTR->NPhotons; i++){
      passing.push_back(ordering.at(i).first);
      if (debug) cout << sel->second.name.Data() << " pho " << ordering.at(i).second << " " << fTR->PhoPt[ordering.at(i).first] << endl;
    }
    for (int i=0; i<(int)(passing.size())-1; i++){
      assert(fTR->PhoPt[passing.at(i)]>=fTR->PhoPt[passing.at(i+1)]);
    }
    passing = PhotonPreSelection(fTR,passing);
    passing = PhotonSelection(fTR,passing);
    
    std::vector<int> passing_jets;
    for (int i=0; i<(int)(ordering_jets.size()); i++) {
      passing_jets.push_back(ordering_jets.at(i).first);
      if (debug) cout << sel->second.name.Data() << " jet " << ordering_jets.at(i).second << " " << fTR->JPt[ordering_jets.at(i).first] << endl;
    }
    for (int i=0; i<(int)(passing_jets.size())-1; i++){
      assert(fTR->JPt[passing_jets.at(i)]>=fTR->JPt[passing_jets.at(i+1)]);
    }
    JetSelection(passing_jets);

    tree_reco_in_acc = passtrigger && StandardEventSelection(fTR,passing,passing_jets);
    if (!tree_reco_in_acc) {passing.clear(); passing_jets.clear();}

    std::vector<int> passing_gen;
    std::vector<int> passing_gen_jets;
    int NGenObjects = 0;
    std::vector<float> GenObjectPt;
    std::vector<float> GenObjectEta;
    std::vector<float> GenObjectPhi;

    if (isdata){
      tree_matched=false;
      tree_gen_in_acc=false;
    }
    else if (thissel!=kLightDiElectron || dataset_id==dy_dataset_id){

    bool isdy = (dataset_id==dy_dataset_id && thissel==kLightDiElectron);
    
    if (!isdy || tree_reco_in_acc){
      
      if (!isdy){
	for (int i=0; i<fTR->NGenPhotons; i++){
	  GenObjectPt.push_back(fTR->GenPhotonPt[i]);
	  GenObjectEta.push_back(fTR->GenPhotonEta[i]);
	  GenObjectPhi.push_back(fTR->GenPhotonPhi[i]);
	  NGenObjects++;
	}
      }
    }
    if (isdy){
      for (int i=0; i<fTR->NGenLeptons; i++){
	if (abs(fTR->GenLeptonID[i])!=11) continue;
	GenObjectPt.push_back(fTR->GenLeptonPt[i]);
	GenObjectEta.push_back(fTR->GenLeptonEta[i]);
	GenObjectPhi.push_back(fTR->GenLeptonPhi[i]);
	NGenObjects++;
      }
    }
    

    {
      std::vector<std::pair<int,float> > ordering_gen;
      for (int i=0; i<NGenObjects; i++){
	ordering_gen.push_back(std::pair<int,float>(i,GenObjectPt[i]));
      }
      std::sort(ordering_gen.begin(),ordering_gen.end(),indexComparator);
      for (size_t i=0; i<ordering_gen.size(); i++){
	passing_gen.push_back(ordering_gen.at(i).first);
      }
    }
    for (int i=0; i<fTR->NGenJets; i++){
      if (i<fTR->NGenJets-1) assert(fTR->GenJetPt[i]>=fTR->GenJetPt[i+1]);
      passing_gen_jets.push_back(i);
    }
    GenJetSelection(passing_gen_jets);

    for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); ){
      bool pass=1;
      if (!isdy) if (fTR->GenPhotonIsoDR04[*it]>5) pass=0;
      if (!pass) it=passing_gen.erase(it); else it++;
    }
    for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); ){
      bool pass=0;
      if (!isdy){
	int mother = fTR->GenPhotonMotherID[*it];
	int mcode=-999;
	if (mother>=-6 && mother<=6) mcode=1;
	if (mother==21) mcode=1;
	if (mother==22 && fTR->GenPhotonMotherStatus[*it]==3) mcode=2;
	if (determine_matchingstatus(mcode,fTR->GenPhotonIsoDR04[*it])==kSignal) pass=1;
      }
      if (isdy) pass=1;
      if (!pass) it=passing_gen.erase(it); else it++;
    }
    for (size_t i=0; i<passing_gen.size(); i++){
      for (vector<int>::iterator it = passing_gen_jets.begin(); it != passing_gen_jets.end(); ){
	bool match = false;
	if (Util::GetDeltaR(GenObjectEta[passing_gen.at(i)],fTR->GenJetEta[*it],GenObjectPhi[passing_gen.at(i)],fTR->GenJetPhi[*it])<0.3) match=true;
	if (match) it=passing_gen_jets.erase(it); else it++;
      }
    }

    tree_matched = false;
    int m0 = -999;
    int m1 = -999;
    if (tree_reco_in_acc){
      bool match0;
      bool match1;
      if (!isdy){
      match0 = (check_matching_status(passing.at(0))==kSignal);
      match1 = (check_matching_status(passing.at(1))==kSignal);
      m0 = (match0) ? fTR->PhoMCmatchindex[passing.at(0)] : -999;
      m1 = (match1) ? fTR->PhoMCmatchindex[passing.at(1)] : -999;
      }
      if (isdy){
	m0 = -999;
	m1 = -999;
	for (size_t i=0; i<passing_gen.size(); i++){
	  if (Util::GetDeltaR(GenObjectEta[passing_gen.at(i)],fTR->PhoEta[passing.at(0)],GenObjectPhi[passing_gen.at(i)],fTR->PhoPhi[passing.at(0)])<0.2) {match0=true; m0=i; break;}
	}
	for (size_t i=0; i<passing_gen.size(); i++){
	  if (Util::GetDeltaR(GenObjectEta[passing_gen.at(i)],fTR->PhoEta[passing.at(1)],GenObjectPhi[passing_gen.at(i)],fTR->PhoPhi[passing.at(1)])<0.2) {match1=true; m1=i; break;}
	}
      }
      if (m0!=m1 && find(passing_gen.begin(),passing_gen.end(),m0)!=passing_gen.end() && find(passing_gen.begin(),passing_gen.end(),m1)!=passing_gen.end()) {
	tree_matched = true;
	for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); ){
	  if (*it==m0 || *it==m1) it++; else it = passing_gen.erase(it);
	}
      }
    }

    if (passing_gen.size()>=2){
      bool bad = false;
      if (GenObjectPt[passing_gen.at(0)]<=40) bad = true;
      if (GenObjectPt[passing_gen.at(1)]<=25) bad = true;
      if (!((fabs(GenObjectEta[passing_gen.at(0)])<1.4442) || (fabs(GenObjectEta[passing_gen.at(0)])>1.566 && fabs(GenObjectEta[passing_gen.at(0)])<2.5))) bad = true;
      if (!((fabs(GenObjectEta[passing_gen.at(1)])<1.4442) || (fabs(GenObjectEta[passing_gen.at(1)])>1.566 && fabs(GenObjectEta[passing_gen.at(1)])<2.5))) bad = true;
      if (bad) passing_gen.clear();
  }

    if (passing_gen.size()>=2 && Util::GetDeltaR(GenObjectEta[passing_gen.at(0)],GenObjectEta[passing_gen.at(1)],GenObjectPhi[passing_gen.at(0)],GenObjectPhi[passing_gen.at(1)])<global_dR_cut_acceptance) passing_gen.clear();

    tree_gen_in_acc = (passing_gen.size()>=2);
    if (!tree_gen_in_acc) {passing_gen.clear(); passing_gen_jets.clear();}

    }

    if (tree_reco_in_acc) {
      FillLead(passing.at(0),passing_jets);
      FillTrail(passing.at(1),passing_jets);
      FillJetsInfo(passing,passing_jets);
      TLorentzVector pho1(fTR->PhoPx[passing.at(0)],fTR->PhoPy[passing.at(0)],fTR->PhoPz[passing.at(0)],fTR->PhoEnergy[passing.at(0)]);
      TLorentzVector pho2(fTR->PhoPx[passing.at(1)],fTR->PhoPy[passing.at(1)],fTR->PhoPz[passing.at(1)],fTR->PhoEnergy[passing.at(1)]);
      float invmass = (pho1+pho2).M();
      dipho_mgg_photon = invmass;
    }

    if (tree_gen_in_acc){
      int m0 = passing_gen.at(0);
      int m1 = passing_gen.at(1);
      if (tree_matched && m0!=passing_gen.at(0)) {int temp = m0; m0=m1; m1=temp;}
      TLorentzVector genphotons[2];
      genphotons[0].SetPtEtaPhiM(GenObjectPt[m0],GenObjectEta[m0],GenObjectPhi[m0],0);
      genphotons[1].SetPtEtaPhiM(GenObjectPt[m1],GenObjectEta[m1],GenObjectPhi[m1],0);
      pholead_GEN_pt =      genphotons[0].Pt();
      photrail_GEN_pt =     genphotons[1].Pt();
      pholead_GEN_eta =     genphotons[0].Eta();
      photrail_GEN_eta =    genphotons[1].Eta();
      pholead_GEN_phi =     genphotons[0].Phi();
      photrail_GEN_phi =    genphotons[1].Phi();
      FillGenJetsInfo(passing_gen,passing_gen_jets);
    }

    if (tree_reco_in_acc || tree_gen_in_acc) sel->second.tree_light->Fill();

  }

};

void DiPhotonMiniTree::End(){
  if (isdata && !isstep2){
    fOutputExtraFile->cd();
    for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++) if (sel->second.tree_extra) sel->second.tree_extra->Write();
  }
  fOutputFile->cd();
  for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++) if (sel->second.tree_full) sel->second.tree_full->Write();
  for (map<EnumSel,Selection>::const_iterator sel=sels.begin(); sel!=sels.end(); sel++) if (sel->second.tree_light) sel->second.tree_light->Write();
  fHNumPU->Write();
  fHNumPU_noweight->Write();
  fHNumPUTrue->Write();
  fHNumPUTrue_noweight->Write();
  fHNumVtx->Write();
  fPNumEvents->Write();
  fPSumWeights->Write();
  fPSumWeights2->Write();

  fOutputFile->Close();
  if (isdata && !isstep2) fOutputExtraFile->Close();

  cout << "FINISHED: processed correctly uuid " << uuid << endl;

}


void DiPhotonMiniTree::FillPhoIso_NewTemplates(TreeReader *fTR, Int_t *n1_arr, Int_t *n2_arr, std::vector<int> passing, SigBkgMode mode, ChoiceMixingTemplates mixing){

  //  if (mixing!=k1Event || mode!=kSigBkg) return;

  //  cout << "EVENT " << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]] << " " << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]] << endl;

  int m1 = (mode==kSigSig || mode==kSigBkg) ? 0 : 1;
  int m2 = (mode==kSigSig || mode==kBkgSig) ? 0 : 1;
  
  int found = 0;

	for (int l=0; l<nclosest_inputmatching; l++){
	  if (found==nclosest) break;
	  Int_t n1 = n1_arr[l];
	  Int_t n2 = n2_arr[l];
	  if (m1==m2 && n1==n2) {cout << "Same event for axis 1 and 2, skipping" << endl; continue;}
	  bool skip_EBEE = false;

	  pfcandidates_struct pfcands;
	  pfcandidates_struct pfcands_1;
	  pfcandidates_struct pfcands_2;

	  pfcandidates_struct *pfcands1=&pfcands;
	  pfcandidates_struct *pfcands2=&pfcands;

	  if (mixing==k1Event && mode==kBkgBkg) {
	    pfcands1=&pfcands_1;
	    pfcands2=&pfcands_2;
	  }

	  // {eta,pt,rho,sigma}
	  float rewinfo_1[6]={-999,-999,-999,-999,-999,-999}; 
	  float rewinfo_2[6]={-999,-999,-999,-999,-999,-999}; 

//	  cout << " m1_" << m1 << " m2_" << m2 << " l_" << l << " mix_" << mixing << " mode_" << mode << endl;
//	  cout << "pho0 " << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]] << " " << fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]] << endl;
//	  cout << "pho1 " << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]] << " " << fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]] << endl;
//
	  if (n1>=0){
	    if (InputTree[m1]->GetEntry(n1)<=0) {cerr << "ERROR: EVENT NOT FOUND IN MATCHINGTREE" << endl; continue;}
	    if (m1==1 && input_event_pass12whoissiglike==0) {
	      float input_photemp_SCeta = input_pholead_SCeta; float input_photemp_SCphi = input_pholead_SCphi; float input_photemp_pt = input_pholead_pt;
	      input_pholead_SCeta = input_photrail_SCeta; input_pholead_SCphi = input_photrail_SCphi; input_pholead_pt = input_photrail_pt;
	      input_photrail_SCeta = input_photemp_SCeta; input_photrail_SCphi = input_photemp_SCphi; input_photrail_pt = input_photemp_pt;
	    }
	    float myrewinfo_1[6]={input_pholead_SCeta,input_photrail_SCeta,input_pholead_pt,input_photrail_pt,input_event_rho,input_event_sigma};
	    memcpy(rewinfo_1,myrewinfo_1,sizeof(myrewinfo_1));
	    for (int k=0; k<input_allphotonpfcand_count; k++){
	      pfcands1->PfCandPt.push_back(input_allphotonpfcand_pt[k]);
	      pfcands1->PfCandEta.push_back(input_allphotonpfcand_eta[k]-input_pholead_SCeta+fTR->SCEta[fTR->PhotSCindex[passing.at(0)]]);
	      pfcands1->PfCandPhi.push_back(DeltaPhiSigned(input_allphotonpfcand_phi[k],input_pholead_SCphi)+fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]]);
	      pfcands1->PfCandVx.push_back(input_allphotonpfcand_vx[k]);
	      pfcands1->PfCandVy.push_back(input_allphotonpfcand_vy[k]);
	      pfcands1->PfCandVz.push_back(input_allphotonpfcand_vz[k]);
	    }
	    std::vector<std::pair<float,float> > obj0;
	    for (int k=0; k<input_vetoobjects_count; k++) obj0.push_back(std::make_pair<float,float>(float(input_vetoobjects_eta[k]),float(input_vetoobjects_phi[k])));

//	      cout << "1_match " << input_pholead_SCeta << " " << input_pholead_SCphi << endl;
//	      cout << "vetoing around " << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]]-fTR->SCEta[fTR->PhotSCindex[passing.at(0)]]+input_pholead_SCeta << " " << fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]]-fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]]+input_pholead_SCphi << endl;

	    if (FindCloseJetsAndPhotons(obj0,fTR->SCEta[fTR->PhotSCindex[passing.at(1)]]-fTR->SCEta[fTR->PhotSCindex[passing.at(0)]]+input_pholead_SCeta,fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]]-fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]]+input_pholead_SCphi)) continue; // veto around the OTHER photon

//	    cout << "veto1 pass" << endl;

	    if ((fabs(fTR->SCEta[fTR->PhotSCindex[passing.at(0)]])<1.4442) && (fabs(input_pholead_SCeta) > 1.4442)) skip_EBEE=true;
	    if ((fabs(fTR->SCEta[fTR->PhotSCindex[passing.at(0)]])>1.4442) && (fabs(input_pholead_SCeta) < 1.4442)) skip_EBEE=true;
	    
//	    if (skip_EBEE) {
//	      cout << "EB/EE migration, skipping1 " << m1 << " " << m2 << " " << l << endl; 
//	      cout << fTR->Run << " " << fTR->LumiSection << " " << fTR->Event << endl;
//	      cout << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]] << " " << input_pholead_SCeta << endl;
//	      continue;
//	    }
	  }
	  
	  skip_EBEE = false;

	  if (n2>=0){
	    if (InputTree[m2]->GetEntry(n2)<=0) {cerr << "ERROR: EVENT NOT FOUND IN MATCHINGTREE" << endl; continue;}
	    if (m2==1 && input_event_pass12whoissiglike==0) {
	      float input_photemp_SCeta = input_pholead_SCeta; float input_photemp_SCphi = input_pholead_SCphi; float input_photemp_pt = input_pholead_pt;
	      input_pholead_SCeta = input_photrail_SCeta; input_pholead_SCphi = input_photrail_SCphi; input_pholead_pt = input_photrail_pt;
	      input_photrail_SCeta = input_photemp_SCeta; input_photrail_SCphi = input_photemp_SCphi; input_photrail_pt = input_photemp_pt;
	    }
	    float myrewinfo_2[6]={input_pholead_SCeta,input_photrail_SCeta,input_pholead_pt,input_photrail_pt,input_event_rho,input_event_sigma};
	    memcpy(rewinfo_2,myrewinfo_2,sizeof(myrewinfo_2));
	    for (int k=0; k<input_allphotonpfcand_count; k++){
	      pfcands2->PfCandPt.push_back(input_allphotonpfcand_pt[k]);
	      pfcands2->PfCandEta.push_back(input_allphotonpfcand_eta[k]-input_pholead_SCeta+fTR->SCEta[fTR->PhotSCindex[passing.at(1)]]);
	      pfcands2->PfCandPhi.push_back(DeltaPhiSigned(input_allphotonpfcand_phi[k],input_pholead_SCphi)+fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]]);
	      pfcands2->PfCandVx.push_back(input_allphotonpfcand_vx[k]);
	      pfcands2->PfCandVy.push_back(input_allphotonpfcand_vy[k]);
	      pfcands2->PfCandVz.push_back(input_allphotonpfcand_vz[k]);
	    }
	    std::vector<std::pair<float,float> > obj1;
	    for (int k=0; k<input_vetoobjects_count; k++) obj1.push_back(std::make_pair<float,float>(float(input_vetoobjects_eta[k]),float(input_vetoobjects_phi[k])));

//	      cout << "2_match " << input_pholead_SCeta << " " << input_pholead_SCphi << endl;
//	      cout << "vetoing around " << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]]-fTR->SCEta[fTR->PhotSCindex[passing.at(1)]]+input_pholead_SCeta << " " << fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]]-fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]]+input_pholead_SCphi << endl;


	    if (FindCloseJetsAndPhotons(obj1,fTR->SCEta[fTR->PhotSCindex[passing.at(0)]]-fTR->SCEta[fTR->PhotSCindex[passing.at(1)]]+input_pholead_SCeta,fTR->SCPhi[fTR->PhotSCindex[passing.at(0)]]-fTR->SCPhi[fTR->PhotSCindex[passing.at(1)]]+input_pholead_SCphi)) continue; // veto around the OTHER photon

//	    cout << "veto2 pass" << endl;

	    if ((fabs(fTR->SCEta[fTR->PhotSCindex[passing.at(1)]])<1.4442) && (fabs(input_pholead_SCeta) > 1.4442)) skip_EBEE=true;
	    if ((fabs(fTR->SCEta[fTR->PhotSCindex[passing.at(1)]])>1.4442) && (fabs(input_pholead_SCeta) < 1.4442)) skip_EBEE=true;
	    
//	    if (skip_EBEE) {
//	      cout << "EB/EE migration, skipping2 " << m1 << " " << m2 << " " << l << endl; 
//	      cout << fTR->Run << " " << fTR->LumiSection << " " << fTR->Event << endl;
//	      cout << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]] << " " << input_pholead_SCeta << endl;
//	      continue;
//	    }
	  }

//	  cout << "filled " << found << endl;

	  if (mixing==k1Event){
	    if (mode!=kBkgBkg){
	      if (mode==kSigSig){
		memcpy(rewinfo_template_1event_sigsig_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
		memcpy(rewinfo_template_1event_sigsig_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
		std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,true,true,rewinfo_1[0],rewinfo_1[0]-pholead_SCeta+photrail_SCeta);
		phoiso_template_1event_sigsig_1[found] = isos.first;
		phoiso_template_1event_sigsig_2[found] = isos.second;
	      }
	      else if (mode==kSigBkg){
		memcpy(rewinfo_template_1event_sigbkg_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
		memcpy(rewinfo_template_1event_sigbkg_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
		std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,true,false,rewinfo_2[0]-photrail_SCeta+pholead_SCeta,rewinfo_2[0]);
		phoiso_template_1event_sigbkg_1[found] = isos.first;
		phoiso_template_1event_sigbkg_2[found] = isos.second;
	      }
	      else if (mode==kBkgSig){
		memcpy(rewinfo_template_1event_bkgsig_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
		memcpy(rewinfo_template_1event_bkgsig_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
		std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,false,true,rewinfo_1[0],rewinfo_1[0]-pholead_SCeta+photrail_SCeta);
		phoiso_template_1event_bkgsig_1[found] = isos.first;
		phoiso_template_1event_bkgsig_2[found] = isos.second;
	      }
	    }
	    else {
	      memcpy(rewinfo_template_1event_bkgbkg_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
	      memcpy(rewinfo_template_1event_bkgbkg_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
	      phoiso_template_1event_bkgbkg_1[found] = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),pfcands1,false,false,rewinfo_1[0]).first;
	      phoiso_template_1event_bkgbkg_2[found] = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),pfcands2,false,false,rewinfo_2[0]).second;
	    }
	  }

	  else if (mixing==k2Events){
	    if (mode==kSigSig){
	      memcpy(rewinfo_template_2events_sigsig_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
	      memcpy(rewinfo_template_2events_sigsig_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
	      std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,true,true,rewinfo_1[0],rewinfo_2[0]);
	      phoiso_template_2events_sigsig_1[found] = isos.first;
	      phoiso_template_2events_sigsig_2[found] = isos.second;
	    }
	    else if (mode==kSigBkg){
	      memcpy(rewinfo_template_2events_sigbkg_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
	      memcpy(rewinfo_template_2events_sigbkg_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
	      std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,true,false,rewinfo_1[0],rewinfo_2[0]);
	      phoiso_template_2events_sigbkg_1[found] = isos.first;
	      phoiso_template_2events_sigbkg_2[found] = isos.second;
	    }
	    else if (mode==kBkgSig){
	      memcpy(rewinfo_template_2events_bkgsig_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
	      memcpy(rewinfo_template_2events_bkgsig_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
	      std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,false,true,rewinfo_1[0],rewinfo_2[0]);
	      phoiso_template_2events_bkgsig_1[found] = isos.first;
	      phoiso_template_2events_bkgsig_2[found] = isos.second;
	    }
	    else if (mode==kBkgBkg){
	      memcpy(rewinfo_template_2events_bkgbkg_1+(found*6),rewinfo_1,sizeof(rewinfo_1));
	      memcpy(rewinfo_template_2events_bkgbkg_2+(found*6),rewinfo_2,sizeof(rewinfo_2));
	      std::pair<float,float> isos = PFPhotonIsolationFromMinitree(passing.at(0),passing.at(1),&pfcands,false,false,rewinfo_1[0],rewinfo_2[0]);
	      phoiso_template_2events_bkgbkg_1[found] = isos.first;
	      phoiso_template_2events_bkgbkg_2[found] = isos.second;
	    }
	  }
	  
	  found++;

	}

	if (found<nclosest) {
	  std::cout << "Found only " << found << " matches in mode " << mode << " mixing " << mixing << std::endl;
	  //	  for (int l=0; l<nclosest; l++) cout << phoiso_template_1event_sigsig_1[l] << " " << phoiso_template_1event_sigsig_2[l] << endl;
	}
	//	else cout << "found ok" << endl;

};


void DiPhotonMiniTree::CorrectPhotonEnergy(StatusScaleUpScaleDown status_escale, StatusScaleUpScaleDown status_esmear, std::vector<pair<int,float> > &neworder, std::vector<float> const &unscaled_r9, std::vector<float> const &unscaled_energy){
  neworder.clear();
  for (int i=0; i<fTR->NPhotons; i++){
    if (debug) cout << "phoenergy " << fTR->PhoPt[i] << " -> ";
    CorrPhoton(fTR,i,unscaled_r9,unscaled_energy,status_escale,status_esmear);
    if (debug) cout << fTR->PhoPt[i] << endl;
    neworder.push_back(std::pair<int,float>(i,fTR->PhoPt[i]));
  }
  sort(neworder.begin(),neworder.end(),indexComparator);
};
StatusScaleUpScaleDown DiPhotonMiniTree::CorrectPhotonR9(StatusScaleUpScaleDown status, std::vector<float> &unscaled_r9){
  assert (unscaled_r9.size()==0 && status==kCorrectedNoShift);
  for (int i=0; i<fTR->NPhotons; i++){
    unscaled_r9.push_back(fTR->PhoR9[i]);
    if (status==kCorrectedNoShift) fTR->PhoR9[i] = R9Rescale(fTR->PhoR9[i],(bool)(fabs(fTR->PhoEta[i])<1.5));
  }
  return status;
};
StatusScaleUpScaleDown DiPhotonMiniTree::CorrectPhotonSieie(StatusScaleUpScaleDown status, std::vector<float> &unscaled_sieie){
  assert (unscaled_sieie.size()==0 && status==kCorrectedNoShift);
  for (int i=0; i<fTR->NPhotons; i++){
    unscaled_sieie.push_back(fTR->PhoSigmaIetaIeta[i]);
    if (status==kCorrectedNoShift) fTR->PhoSigmaIetaIeta[i] = SieieRescale(fTR->PhoSigmaIetaIeta[i],(bool)(fabs(fTR->PhoEta[i])<1.5));
  }
  return status;
};

void DiPhotonMiniTree::CorrPhoton(TreeReader *fTR, int i, std::vector<float> const &unscaled_r9, std::vector<float> const &unscaled_energy, StatusScaleUpScaleDown status_escale, StatusScaleUpScaleDown status_esmear){

  float invpresentcorr = unscaled_energy[i]/fTR->PhoEnergy[i];
  fTR->PhoPt[i]*=invpresentcorr;
  fTR->PhoPx[i]*=invpresentcorr;
  fTR->PhoPy[i]*=invpresentcorr;
  fTR->PhoPz[i]*=invpresentcorr;
  fTR->PhoEnergy[i]*=invpresentcorr;

  float corr = fTR->PhoRegrEnergy[i]/fTR->PhoEnergy[i];
  
  if (isdata && status_escale!=kUncorrected) {
    float sigmas=0;
    if (status_escale==kShiftUp) sigmas=1;
    if (status_escale==kShiftDown) sigmas=-1;
    corr*=EnergyScaleOffset(fTR->PhoSCEta[i],fTR->PhoR9[i],fTR->Run,sigmas);
  }
  if (!isdata && status_esmear!=kUncorrected) {
    float r9forsmearing = (year==2012) ? unscaled_r9.at(i) : fTR->PhoR9[i];
    float sigmas=0;
    if (status_esmear==kShiftUp) sigmas=1;
    if (status_esmear==kShiftDown) sigmas=-1;
    float width = EnergySmearingCorrection(fTR->PhoSCEta[i],r9forsmearing,fTR->Run,sigmas)/100;
    // implements deterministic smearing
    UInt_t seedBase = (UInt_t) fTR->Event + (UInt_t) fTR->Run + (UInt_t) fTR->LumiSection;
    UInt_t seed1    = seedBase + 100000*(UInt_t) (TMath::Abs(100.*fTR->PhoPhi[i])) + 1000*(UInt_t) (TMath::Abs(100.*fTR->PhoEta[i]));
    float mycorr = 1;
    randomgen_forEsmearing->SetSeed(seed1);
    if (width>0) mycorr=randomgen_forEsmearing->Gaus(1.,width);
    corr*=mycorr;
    if (debug) cout << "esmear " << mycorr << " width " << width << endl;
  }

  if (debug) cout << corr << " " << status_escale << " " << status_esmear << " -> ";

  fTR->PhoPt[i]*=corr;
  fTR->PhoPx[i]*=corr;
  fTR->PhoPy[i]*=corr;
  fTR->PhoPz[i]*=corr;
  fTR->PhoEnergy[i]*=corr;

};

std::vector<int> DiPhotonMiniTree::ApplyPixelVeto(TreeReader *fTR, vector<int> passing, bool forelectron){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Pixel veto (conversion safe)
    int wantpixelseed;
    if (forelectron) wantpixelseed=1; else wantpixelseed=0; 
    if (fTR->PhoHasPixSeed[*it]!=wantpixelseed) it=passing.erase(it); else it++;
    //    if (fTR->PhoPassConversionVeto[*it]==wantpixelseed) it=passing.erase(it); else it++;
  }

  return passing;

};

float DiPhotonMiniTree::getEtaCorrectionBarrel(float eta){
  int iEta = (int)(TMath::Abs(eta)*(5/0.087));
  if ( iEta < 40.2198 ) return 1;
  return 1.0/(1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198));
};

std::vector<int> DiPhotonMiniTree::PhotonPreSelection(TreeReader *fTR, std::vector<int> passing){

  // FILTERS
  // scraping veto done at ntuplizer level
  if (!PassPrimaryVertexFilter()) passing.clear();
  if (!fTR->HBHENoiseFilterResult) passing.clear();
  if (!fTR->CSCTightHaloID) passing.clear();

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // SC matching
    if (fTR->PhotSCindex[*it]<0) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // fiducial region
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if ((fabs(eta)>1.4442 && fabs(eta)<1.566) || (fabs(eta)>2.5)) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    float eta=fTR->PhoEta[*it];
    if (fabs(eta)>2.5) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Pt cut on RawEnCetaCorr
    float energy=fTR->SCRaw[fTR->PhotSCindex[*it]];
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    if (fabs(eta)<1.4442) energy*=getEtaCorrectionBarrel(eta);
    if (fabs(eta)>1.566) energy+=fTR->SCPre[fTR->PhotSCindex[*it]];
    if (energy/cosh(eta)<20) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // Remove eta/phi cracks
//    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//    float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
//    if ((phocorr->isInPhiCracks(phi,eta)) || (phocorr->isInEBEtaCracks(eta))) it=passing.erase(it); else it++;
//  }

  // MVA presel from Hgg

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // HoverE cut
    float r9=fTR->PhoR9[*it];
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float hoe=fTR->PhoHoverE[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && r9<0.9 && hoe<0.075) pass=1;
    if (fabs(eta)<1.4442 && r9>0.9 && hoe<0.082) pass=1;
    if (fabs(eta)>1.566 && r9<0.9 && hoe<0.075) pass=1;
    if (fabs(eta)>1.566 && r9>0.9 && hoe<0.075) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // sieie cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float sieie=fTR->PhoSigmaIetaIeta[*it];
    float sipip=sqrt(fTR->PhoSigmaIphiIphi[*it]); // to be fixed in the producer
    bool pass=0;
    if (fabs(eta)<1.4442 && sieie<0.014 && sieie>0.001 && sipip>0.001) pass=1;
    if (fabs(eta)>1.566 && sieie<0.034) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (trigger)
    float r9=fTR->PhoR9[*it];
    bool pass=0;
    float etcorrecaliso=fTR->PhoIso03Ecal[*it]-0.012*fTR->PhoPt[*it];
    float etcorrhcaliso=fTR->PhoIso03Hcal[*it]-0.005*fTR->PhoPt[*it];
    float etcorrtrkiso=fTR->PhoIso03TrkHollow[*it]-0.002*fTR->PhoPt[*it];
    if (r9<0.9 && etcorrecaliso<4 && etcorrhcaliso<4 && etcorrtrkiso<4) pass=1;
    if (r9>0.9 && etcorrecaliso<50 && etcorrhcaliso<50 && etcorrtrkiso<50) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (filter)
    bool pass=0;
    if(fTR->PhoCiCPFIsoChargedDR03[*it]<9.) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::MuonSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    float eta=fTR->MuEta[*it];
    if ((fabs(eta)>1.4442 && fabs(eta)<1.566) || (fabs(eta)>2.5)) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    if (fTR->MuPt[*it]<10) it=passing.erase(it); else it++;
  }

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    if (!(fTR->MuIsGlobalMuon[*it] || fTR->MuIsTrackerMuon[*it])) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    if (fTR->MuRelIso03[*it]>0.15) it=passing.erase(it); else it++;
//  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::GenLevelIsolationCut(TreeReader *fTR, std::vector<int> passing){

  cout << "DO NOT USE, SUPERSEDED" << endl;
  assert(false);

  // Genlevel isolation cut for MC
  // selects only photons matched to gen level, with gen-level iso < 5 GeV

  if (isdata) {
    std::cout << "Calling gen level isolation cut on data!!!" << std::endl;
    return passing;
  }
  
  return std::vector<int>();

};

std::vector<int> DiPhotonMiniTree::PhotonSelection(TreeReader *fTR, std::vector<int> passing, TString mode){

  if (mode!="" && mode!="invert_sieie_cut" && mode!="revert_pixel_veto"){
    std::cout << "Error" << std::endl;
    return std::vector<int>();
  }

  if (mode=="revert_pixel_veto") passing = ApplyPixelVeto(fTR,passing,1); // for DY pixel veto reverse
  else passing = ApplyPixelVeto(fTR,passing,0);


  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // HoverE cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float hoe=fTR->PhoHoverE[*it];
    bool pass=0;
    if (fabs(eta)<1.4442 && hoe<0.05) pass=1;
    if (fabs(eta)>1.566 && hoe<0.05) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // sieie cut
    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
    float sieie=fTR->PhoSigmaIetaIeta[*it];
    bool pass=0;
    if (mode=="invert_sieie_cut"){ // sieie sideband
      if (fabs(eta)<1.4442 && sieie>0.011 && sieie<0.014) pass=1;
      if (fabs(eta)>1.566 && sieie>0.030 && sieie<0.034) pass=1;
    }
    else {
    if (fabs(eta)<1.4442 && sieie<0.011) pass=1;
    if (fabs(eta)>1.566 && sieie<0.030) pass=1;
    }
    if (!pass) it=passing.erase(it); else it++;
  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    bool pass=1;
//    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//    float phi=fTR->SCPhi[fTR->PhotSCindex[*it]];
//    for (int i=0; i<fTR->NMus; i++){
//      float mueta = fTR->MuEta[i];
//      float muphi = fTR->MuPhi[i];
//      if (Util::GetDeltaR(mueta,eta,muphi,phi)<0.4) pass=0;
//    }
//    if (!pass) it=passing.erase(it); else it++;
//  }

//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    bool pass=0;
//    float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//    float combiso = PFIsolation(*it,0,"combined");
//    if (mode=="no_combiso_cut") pass=1; // pass in any case
//    else if (mode=="cut_combiso_sideband"){ // selection for sideband
//      //      if (combiso> && combiso<) pass=1;
//    }
//    else if (combiso<999) pass=1;
//    if (!pass) it=passing.erase(it); else it++;
//  }


//  if (mode=="invert_sieie_cut"){ // veto close objects
//    for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // HoverE cut
//      bool pass=1;
//      if (FindCloseJetsAndPhotons(fTR,0,*it,"exclude_object_itself")) pass=0;
//      if (!pass) it=passing.erase(it); else it++;
//    }
//  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::SignalSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (check_matching_status(*it)==kSignal) pass=1;
    if (pass) {
      float dR = Util::GetDeltaR(fTR->GenPhotonEta[fTR->PhoMCmatchindex[*it]],fTR->PhoEta[*it],fTR->GenPhotonPhi[fTR->PhoMCmatchindex[*it]],fTR->PhoPhi[*it]);
      if (dR>0.1) cout << "PATHOLOGICAL MATCHING DR " << dR << endl;
    }
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::BackgroundSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (check_matching_status(*it)==kBackground) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  return passing;

};

std::vector<int> DiPhotonMiniTree::DiPhotonInvariantMassCutSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoPt[*it]>25) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  if (passing.size()<2) return std::vector<int>();

  passing.resize(2);

  TLorentzVector pho1(fTR->PhoPx[passing.at(0)],fTR->PhoPy[passing.at(0)],fTR->PhoPz[passing.at(0)],fTR->PhoEnergy[passing.at(0)]);
  TLorentzVector pho2(fTR->PhoPx[passing.at(1)],fTR->PhoPy[passing.at(1)],fTR->PhoPz[passing.at(1)],fTR->PhoEnergy[passing.at(1)]);
  float invmass0 = (pho1+pho2).M();
  if (fabs(invmass0-91.2)>10) return std::vector<int>();

  return passing;

};

bool DiPhotonMiniTree::SinglePhotonEventSelection(TreeReader *fTR, std::vector<int> &passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (fTR->PhoPt[*it]>25) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

  if (passing.size()<1) return false;

  return true;

};

bool DiPhotonMiniTree::StandardEventSelection(TreeReader *fTR, std::vector<int> &passing, std::vector<int> &passing_jets){

  if (passing.size()<2) return false;

  passing.resize(2); // keep only the first two

  TLorentzVector pho1(fTR->PhoPx[passing.at(0)],fTR->PhoPy[passing.at(0)],fTR->PhoPz[passing.at(0)],fTR->PhoEnergy[passing.at(0)]);
  TLorentzVector pho2(fTR->PhoPx[passing.at(1)],fTR->PhoPy[passing.at(1)],fTR->PhoPz[passing.at(1)],fTR->PhoEnergy[passing.at(1)]);
  float invmass0 = (pho1+pho2).M();
  float deta=fTR->PhoEta[passing.at(0)]-fTR->PhoEta[passing.at(1)];
  float dphi=Util::DeltaPhi(fTR->PhoPhi[passing.at(0)],fTR->PhoPhi[passing.at(1)]);
  double dR=sqrt(dphi*dphi+deta*deta);

  if (fTR->PhoPt[passing.at(0)]<40) return false;
  if (fTR->PhoPt[passing.at(1)]<25) return false;
  //  if (invmass0<80) return false;
  if (dR<global_dR_cut_acceptance) return false;

  if (VetoJetPhotonOverlap(passing,passing_jets)) return false;

  return true;
};

bool DiPhotonMiniTree::VetoJetPhotonOverlap(std::vector<int> &passing, std::vector<int> &passing_jets){

  // remove jets that coincide with selected photons
  for (size_t i=0; i<passing.size(); i++){
    int m = photon_jet_matching.at(passing.at(i)).m_jet;
    for (vector<int>::iterator it = passing_jets.begin(); it != passing_jets.end(); ){
      if ((*it)==m) it=passing_jets.erase(it); else it++;
    }
  }

  bool out=false;

//  // minimum dR between photon and jet
//  for (size_t i=0; i<passing.size(); i++){
//    for (vector<int>::iterator it = passing_jets.begin(); it != passing_jets.end(); ){
//      float dR = Util::GetDeltaR(fTR->PhoEta[passing.at(i)],fTR->JEta[*it],fTR->PhoPhi[passing.at(i)],fTR->JPhi[*it]);
//      if (dR<global_mindR_photon_jet) {out=true; it=passing_jets.erase(it);} else it++;
//    }
//  }

  return out;

};

bool DiPhotonMiniTree::DiMuonFromZSelection(TreeReader *fTR, std::vector<int> &passing){

  if (passing.size()<2) return false;

  passing.resize(2); // keep only the first two

  TLorentzVector mu1(fTR->MuPx[passing.at(0)],fTR->MuPy[passing.at(0)],fTR->MuPz[passing.at(0)],fTR->MuE[passing.at(0)]);
  TLorentzVector mu2(fTR->MuPx[passing.at(1)],fTR->MuPy[passing.at(1)],fTR->MuPz[passing.at(1)],fTR->MuE[passing.at(1)]);

  float invmass0 = (mu1+mu2).M();

  if (fabs(invmass0-91.2)>10) return false;

  return true;
};

bool DiPhotonMiniTree::TriggerSelection(){
  if (!isdata) return true; // trigger sel off in MC
#include "DiPhotonTriggerSelection.cc"
};

int DiPhotonMiniTree::Count_part_isrfsr_gamma(TreeReader *fTR){

  if (isdata){
    std::cout << "calling count of gammas for k-factor on data!!!" << std::endl;
    return 0;
  }

  int res=0;

  for (int i=0; i<fTR->NGenPhotons; i++){
    bool pass=0;
    int mother = fTR->GenPhotonMotherID[i];
    if (mother>=-6 && mother<=6) pass=1;
    if (mother==21) pass=1;
    if (mother==22 && fTR->GenPhotonMotherStatus[i]==3) pass=1;
    if (pass) res++;
  }

  if (res>2) res=2;
  return res;

};


double DiPhotonMiniTree::phiNorm(float phi) {

  const float pi = 3.1415927;
  const float twopi = 2.0*pi;

  if(phi >  pi) {phi = phi - twopi;}
  if(phi < -pi) {phi = phi + twopi;}

  return phi;
}

bool DiPhotonMiniTree::FindCloseJetsAndPhotons(TreeReader *fTR, float rotation_phi, int phoqi, TString mod){

  if (mod!="" && mod!="exclude_object_itself") {std::cout << "error" << std::endl; return true;}

  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);
  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    photon_position *= r;
  }
  double eta = photon_position.Eta();
  double phi = photon_position.Phi();
  
  //  if (debug) std::cout << "calling FindCloseJetsAndPhotons eta=" << eta << " phi=" << phi << std::endl;

  const float mindR = 0.8;
  bool found=false;

  for (int i=0; i<fTR->NJets; i++){
    if (fTR->JPt[i]<20) continue;
    if (!(fTR->JEcorr.at(i)>0)) continue;
    //    cout << "vj " << fTR->JPt[i] << " " << fTR->JEta[i] << " " << fTR->JPhi[i] << endl;
    //    if (!(fTR->JPassPileupIDT0[fTR->JVrtxListStart[i]+0])) continue;
    float dR = Util::GetDeltaR(eta,fTR->JEta[i],phi,fTR->JPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    if (dR<mindR) found=true;
    //    if (debug) if (dR<mindR) std::cout << "Found jet eta=" << fTR->JEta[i] << " phi=" << fTR->JPhi[i] << std::endl;
  }

  for (int i=0; i<fTR->NPhotons; i++){
    if (fTR->PhoPt[i]<10) continue;
    //    cout << "vg " << fTR->PhoPt[i] << " " << fTR->PhoEta[i] << " " << fTR->PhoPhi[i] << endl;
    float dR = Util::GetDeltaR(eta,fTR->PhoEta[i],phi,fTR->PhoPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    if (dR<mindR) found=true;
    //    if (debug) if (dR<mindR) std::cout << "Found phot eta=" << fTR->PhoEta[i] << " phi=" << fTR->PhoPhi[i] << std::endl;
  }

  //  if (debug) std::cout << "returning " << found << std::endl;
  return found;

};

bool DiPhotonMiniTree::FindCloseJetsAndPhotons(std::vector<std::pair<float,float> > obj, float eta, float phi){

  const float mindR = 0.8;
  bool found=false;

  //  cout << "looking for close obj at " << eta << " " << phi << endl;

  for (int i=0; i<obj.size(); i++){
    float dR = Util::GetDeltaR(eta,obj.at(i).first,phi,obj.at(i).second);
    if (dR<mindR) {
      found=true;
      //      cout << "found close " << obj.at(i).first << " " << obj.at(i).second << " dR=" << dR << endl;
    }
  }

  return found;

};

void DiPhotonMiniTree::FillVetoObjects(TreeReader *fTR, int phoqi, TString mod){

  //  if (debug) cout << "start FillVetoObjects" << endl;

  std::vector<std::pair<TVector3,int>> obj;

  if (mod!="" && mod!="exclude_object_itself") {std::cout << "error" << std::endl;}

  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);

  double eta = photon_position.Eta();
  double phi = photon_position.Phi();
  
  for (int i=0; i<fTR->NJets; i++){
    if (fTR->JPt[i]<20) continue;
    if (!(fTR->JEcorr.at(i)>0)) continue;
    float dR = Util::GetDeltaR(eta,fTR->JEta[i],phi,fTR->JPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    TVector3 a;
    a.SetPtEtaPhi(fTR->JPt[i],fTR->JEta[i],fTR->JPhi[i]);
    obj.push_back(std::pair<TVector3,int>(a,0));
  }

  //  if (debug) cout << "bla1" <<endl;

  for (int i=0; i<fTR->NPhotons; i++){
    if (fTR->PhoPt[i]<10) continue;
    float dR = Util::GetDeltaR(eta,fTR->PhoEta[i],phi,fTR->PhoPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    TVector3 a;
    a.SetPtEtaPhi(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i]);
    obj.push_back(std::pair<TVector3,int>(a,1));
  }

  //  if (debug) cout << "bla2" <<endl;

  if (obj.size()>global_maxN_vetoobjects) {std::cout << "MaxN vetoobjects reached" << std::endl; obj.resize(global_maxN_vetoobjects);}
  for (int i=0; i<obj.size(); i++){
    vetoobjects_pt[i]=obj.at(i).first.Pt();
    vetoobjects_eta[i]=obj.at(i).first.Eta();
    vetoobjects_phi[i]=obj.at(i).first.Phi();
    vetoobjects_type[i]=obj.at(i).second;
  }
  vetoobjects_count = obj.size();

  //  if (debug) cout << "end FillVetoObjects" << endl;

};

std::set<int> DiPhotonMiniTree::GetPFCandIDedRemovals(TreeReader *fTR, int phoqi){
  std::set<int> out;
  if (fTR->PhoMatchedPFPhotonCand.size()>0) out.insert(fTR->PhoMatchedPFPhotonCand[phoqi]);
  if (fTR->PhoMatchedPFElectronCand.size()>0) out.insert(fTR->PhoMatchedPFElectronCand[phoqi]);
  if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) out.insert(fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]);
  return out;
};

std::set<int> DiPhotonMiniTree::GetPFCandInsideFootprint(TreeReader *fTR, int phoqi, float rotation_phi, TString component){
  std::set<int> removals = GetPFCandWithFootprintRemoval(fTR,phoqi,rotation_phi,false,component);
  return removals;
};

std::set<int> DiPhotonMiniTree::GetPFCandInsideFootprint(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, TString component){
  return GetPFCandWithFootprintRemoval(fTR,pfcands,phoqi,rotation_phi,false,component);
};

std::set<int> DiPhotonMiniTree::GetPrecalculatedFootprintPhoEl(int phoqi){
  int stop = (phoqi<fTR->NPhotons-1) ? fTR->PhoFootprintPfCandsListStart[phoqi+1] : fTR->PhoFootprintPfCands.size();
  std::set<int> out;
  for (int i = fTR->PhoFootprintPfCandsListStart[phoqi]; i<stop; i++){
    int pfcand = fTR->PhoFootprintPfCands.at(i);
    int type = FindPFCandType(fTR->PfCandPdgId[pfcand]);
    if (type!=2 && type!=3) continue;
    out.insert(pfcand);
  }
  return out;
};

std::set<int> DiPhotonMiniTree::GetPFCandWithFootprintRemoval(TreeReader *fTR, int phoqi, float rotation_phi, bool outoffootprint, TString component){

  if (component!="photon"){
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return std::set<int>();
  }

  if (component!="neutral" && component!="charged" && component!="photon" && component!="combined") {
    std::cout << "Wrong choice for component" << std::endl;
    return std::set<int>();
  }

  int scindex = fTR->PhotSCindex[phoqi];
  
  if (scindex<0) {
    std::cout << "Error in GetPFCandOverlappingSC" << std::endl;
    std::cout << scindex << " " << phoqi << " " << fTR->PhoPt[phoqi] << " " << fTR->PhoEta[phoqi] << " " << fTR->PhoPhi[phoqi] << std::endl;
    for (int i=0; i<fTR->NSuperClusters; i++){
      cout << i << " " << fTR->SCEnergy[i]/TMath::CosH(fTR->SCEta[phoqi]) << " " << fTR->SCEta[phoqi] << " " << fTR->SCRaw[i]/fTR->SCEnergy[i] << endl;
    }
    return std::set<int>();
  }

  bool isbarrel = fTR->PhoisEB[phoqi];
  int nxtals = fTR->SCNXtals[scindex];

  std::set<int> result;

  for (int i=0; i<fTR->NPfCand; i++){

    int type = FindPFCandType(fTR->PfCandPdgId[i]);
    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    TVector3 sc_position = TVector3(fTR->SCX[scindex],fTR->SCY[scindex],fTR->SCZ[scindex]);

    bool inside=false;

    if (fTR->PhoMatchedPFPhotonCand.size()>0)           if (fTR->PhoMatchedPFPhotonCand[phoqi]==i) continue;
    if (fTR->PhoMatchedPFElectronCand.size()>0)         if (fTR->PhoMatchedPFElectronCand[phoqi]==i) continue;
    if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) if (fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]==i) continue;

    for (int j=0; j<nxtals; j++){
      
      TVector3 xtal_position = TVector3(fTR->XtalX[fTR->SCXtalListStart[scindex]+j],fTR->XtalY[fTR->SCXtalListStart[scindex]+j],fTR->XtalZ[fTR->SCXtalListStart[scindex]+j]);

      if (rotation_phi!=0) {
	TRotation r; r.RotateZ(rotation_phi);
	xtal_position *= r;
      }

      TVector3 ecalpfhit = PropagatePFCandToEcal(fTR,i,isbarrel ? xtal_position.Perp() : xtal_position.z(), isbarrel); // this would be the most correct

      if (isbarrel){
	float xtalEtaWidth = fTR->XtalEtaWidth[fTR->SCXtalListStart[scindex]+j]*(1+global_linkbyrechit_enlargement);
	float xtalPhiWidth = fTR->XtalPhiWidth[fTR->SCXtalListStart[scindex]+j]*(1+global_linkbyrechit_enlargement);
	if (fabs(ecalpfhit.Eta()-xtal_position.Eta())<xtalEtaWidth/2 && Util::DeltaPhi(ecalpfhit.Phi(),xtal_position.Phi())<xtalPhiWidth/2) inside=true;
      }
      else { // EE
	if (ecalpfhit.z()*xtal_position.z()>0){
	  TVector3 xtal_corners[4];
	  xtal_corners[0]=TVector3(fTR->XtalFront1X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront1Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront1Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[1]=TVector3(fTR->XtalFront2X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront2Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront2Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[2]=TVector3(fTR->XtalFront3X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront3Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront3Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[3]=TVector3(fTR->XtalFront4X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront4Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront4Z[fTR->SCXtalListStart[scindex]+j]);
	  if (rotation_phi!=0) {
	    TRotation r; r.RotateZ(rotation_phi);
	    for (int k=0; k<4; k++) xtal_corners[k] *= r;
	  }
	  float hitx = ecalpfhit.x();
	  float hity = ecalpfhit.y();
	  float polx[5];
	  float poly[5];
	  for (int k=0; k<4; k++) polx[k] = xtal_corners[k].x();
	  for (int k=0; k<4; k++) poly[k] = xtal_corners[k].y();
	  polx[4]=polx[0]; poly[4]=poly[0]; // closed polygon
	  float centerx = (polx[0]+polx[1]+polx[2]+polx[3])/4;
	  float centery = (poly[0]+poly[1]+poly[2]+poly[3])/4;
	  hitx = centerx + (hitx-centerx)/(1.0+global_linkbyrechit_enlargement);
	  hity = centery + (hity-centery)/(1.0+global_linkbyrechit_enlargement);
	  if (TMath::IsInside(hitx,hity,5,polx,poly)) inside=true;
	}
      }

    }

    if (outoffootprint) inside=!inside;
    if (inside) result.insert(i);

  }

  return result;

};

std::set<int> DiPhotonMiniTree::GetPFCandWithFootprintRemoval(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, bool outoffootprint, TString component){

  assert(rotation_phi==0);

  if (component!="photon"){
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return std::set<int>();
  }

  if (phoqi<0) return std::set<int>();

  int scindex = fTR->PhotSCindex[phoqi];
  
  if (scindex<0) {
    std::cout << "Error in GetPFCandOverlappingSC" << std::endl;
    std::cout << scindex << " " << phoqi << std::endl;
    return std::set<int>();
  }

  bool isbarrel = fTR->PhoisEB[phoqi];
  int nxtals = fTR->SCNXtals[scindex];

  std::set<int> result;

  for (int i=0; i<pfcands->PfCandPt.size(); i++){

    TVector3 sc_position = TVector3(fTR->SCX[scindex],fTR->SCY[scindex],fTR->SCZ[scindex]);

    bool inside=false;

    for (int j=0; j<nxtals; j++){
      
      TVector3 xtal_position = TVector3(fTR->XtalX[fTR->SCXtalListStart[scindex]+j],fTR->XtalY[fTR->SCXtalListStart[scindex]+j],fTR->XtalZ[fTR->SCXtalListStart[scindex]+j]);

      TVector3 ecalpfhit = PropagatePFCandToEcal(pfcands,i,isbarrel ? xtal_position.Perp() : xtal_position.z(), isbarrel); // this would be the most correct

      if (isbarrel){
	float xtalEtaWidth = fTR->XtalEtaWidth[fTR->SCXtalListStart[scindex]+j]*(1+global_linkbyrechit_enlargement);
	float xtalPhiWidth = fTR->XtalPhiWidth[fTR->SCXtalListStart[scindex]+j]*(1+global_linkbyrechit_enlargement);
	if (fabs(ecalpfhit.Eta()-xtal_position.Eta())<xtalEtaWidth/2 && Util::DeltaPhi(ecalpfhit.Phi(),xtal_position.Phi())<xtalPhiWidth/2) inside=true;
      }
      else { // EE
	if (ecalpfhit.z()*xtal_position.z()>0){
	  TVector3 xtal_corners[4];
	  xtal_corners[0]=TVector3(fTR->XtalFront1X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront1Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront1Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[1]=TVector3(fTR->XtalFront2X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront2Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront2Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[2]=TVector3(fTR->XtalFront3X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront3Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront3Z[fTR->SCXtalListStart[scindex]+j]);
	  xtal_corners[3]=TVector3(fTR->XtalFront4X[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront4Y[fTR->SCXtalListStart[scindex]+j],fTR->XtalFront4Z[fTR->SCXtalListStart[scindex]+j]);
	  if (rotation_phi!=0) {
	    TRotation r; r.RotateZ(rotation_phi);
	    for (int k=0; k<4; k++) xtal_corners[k] *= r;
	  }
	  float hitx = ecalpfhit.x();
	  float hity = ecalpfhit.y();
	  float polx[5];
	  float poly[5];
	  for (int k=0; k<4; k++) polx[k] = xtal_corners[k].x();
	  for (int k=0; k<4; k++) poly[k] = xtal_corners[k].y();
	  polx[4]=polx[0]; poly[4]=poly[0]; // closed polygon
	  float centerx = (polx[0]+polx[1]+polx[2]+polx[3])/4;
	  float centery = (poly[0]+poly[1]+poly[2]+poly[3])/4;
	  hitx = centerx + (hitx-centerx)/(1.0+global_linkbyrechit_enlargement);
	  hity = centery + (hity-centery)/(1.0+global_linkbyrechit_enlargement);
	  if (TMath::IsInside(hitx,hity,5,polx,poly)) inside=true;
	}
      }

    }

    if (outoffootprint) inside=!inside;
    if (inside) result.insert(i);

  }

  return result;

};

TVector3 DiPhotonMiniTree::PropagatePFCandToEcal(TreeReader *fTR, int pfcandindex, float position, bool isbarrel){
  // WARNING: this propagates until EE+ or EE- at the given TMath::Abs(position.z()) for isbarrel=0, depending on where the candidate is pointing.

  int i = pfcandindex;

  if (FindPFCandType(fTR->PfCandPdgId[i])!=2) {
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return TVector3(0,0,0);
  }

  TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
  TVector3 pfmomentum;
  pfmomentum.SetPtEtaPhi(fTR->PfCandPt[i],fTR->PfCandEta[i],fTR->PfCandPhi[i]);
  pfmomentum = pfmomentum.Unit();
  TVector3 ecalpfhit(0,0,0);
  if (isbarrel){
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    double dist = TGeoTube::DistFromInsideS(p,d,0,position,1e+10);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }
  else { // EE
    double dim[6]={1e+10,1e+10,fabs(position),0,0,0};
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    eegeom.SetDimensions(dim);
    double dist = eegeom.DistFromInside(p,d);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }

  return ecalpfhit;

};

TVector3 DiPhotonMiniTree::PropagatePFCandToEcal(pfcandidates_struct *pfcands, int pfcandindex, float position, bool isbarrel){
  // WARNING: this propagates until EE+ or EE- at the given TMath::Abs(position.z()) for isbarrel=0, depending on where the candidate is pointing.

  int i = pfcandindex;

  TVector3 pfvertex(pfcands->PfCandVx[i],pfcands->PfCandVy[i],pfcands->PfCandVz[i]);
  TVector3 pfmomentum; pfmomentum.SetPtEtaPhi(pfcands->PfCandPt[i],pfcands->PfCandEta[i],pfcands->PfCandPhi[i]);
  pfmomentum = pfmomentum.Unit();
  TVector3 ecalpfhit(0,0,0);
  if (isbarrel){
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    double dist = TGeoTube::DistFromInsideS(p,d,0,position,1e+10);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }
  else { // EE
    double dim[6]={1e+10,1e+10,fabs(position),0,0,0};
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    eegeom.SetDimensions(dim);
    double dist = eegeom.DistFromInside(p,d);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }

  return ecalpfhit;

};

isolations_struct DiPhotonMiniTree::RandomConeIsolation(TreeReader *fTR, int phoqi, TString mod){

  float result=0;
  const double pi = TMath::Pi();

  double rotation_phi = pi/2;

  //  double rotation_phi = fTR->PhoSCRemovalRConePhi[phoqi]-TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]).Phi();

//  if (fTR->PhoSCRemovalPFIsoPhoton[phoqi]>900){
//    isolations_struct out;
//    out.nphotoncand=0; out.nchargedcand=0; out.nneutralcand=0;
//    out.photon = -999;
//    out.charged = -999;
//    out.neutral = -999;
//    out.newphi=-999;
//    return out;
//  }

  bool isok = !(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
  if (!isok) {
    
    //    cout << "WRONG: MISMATCH IN OUTPUT OF CHECK " << fTR->PhoSCEta[phoqi] << " " << fTR->PhoPt[phoqi] << endl;
    
    rotation_phi = -pi/2;
    isok=!(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
  }
  
  int count=0;
  while (!isok && count<20) {
    rotation_phi = randomgen->Uniform(0.8,2*pi-0.8);
    isok=!(FindCloseJetsAndPhotons(fTR,rotation_phi,phoqi,mod));
    count++;
  }

  isolations_struct out;
  out.nphotoncand=0; out.nchargedcand=0; out.nneutralcand=0;


  if (count==20){
    std::cout << "Error in random cone generation!!!"  << std::endl;
    out.photon = -999;
    out.charged = -999;
    out.neutral = -999;
    out.newphi=-999;
    return out;
  };


  out.photon = PFIsolation(phoqi,rotation_phi,"photon",&(out.nphotoncand),&(out.photoncandenergies),&(out.photoncandets),&(out.photoncanddetas),&(out.photoncanddphis),&(out.newphi));
  out.charged = PFIsolation(phoqi,rotation_phi,"charged",&(out.nchargedcand),&(out.chargedcandenergies),&(out.chargedcandets),&(out.chargedcanddetas),&(out.chargedcanddphis),&(out.newphi));
  out.neutral = PFIsolation(phoqi,rotation_phi,"neutral",&(out.nneutralcand),&(out.neutralcandenergies),&(out.neutralcandets),&(out.neutralcanddetas),&(out.neutralcanddphis),&(out.newphi));

//  if (fabs(out.photon/fTR->PhoSCRemovalPFIsoPhotonRCone[phoqi]-1)>1e-3)  {
//    cout << "(" << out.photon << "," << fTR->PhoSCRemovalPFIsoPhotonRCone[phoqi] << ") ";
//    cout << "(" << fTR->PhoSCRemovalRConePhi[phoqi] << "," << out.newphi << ") ";
//    cout << "(" << fTR->PhoSCRemovalRConeEta[phoqi] << "," << fTR->PhoSCEta[phoqi] << ") ";
//    cout << "(" << pholead_PhoSCRemovalPFIsoPhoton << "," << fTR->PhoSCRemovalPFIsoPhoton[phoqi] << ") ";
//    cout << endl;
//  }

  return out;

};

isolations_struct DiPhotonMiniTree::PFConeIsolation(TreeReader *fTR, int phoqi){
  isolations_struct out;
  out.nphotoncand=0; out.nchargedcand=0; out.nneutralcand=0; out.newphi=-999;
  out.photon = PFIsolation(phoqi,0,"photon",&(out.nphotoncand),&(out.photoncandenergies),&(out.photoncandets),&(out.photoncanddetas),&(out.photoncanddphis));
  out.charged = PFIsolation(phoqi,0,"charged",&(out.nchargedcand),&(out.chargedcandenergies),&(out.chargedcandets),&(out.chargedcanddetas),&(out.chargedcanddphis));
  out.neutral = PFIsolation(phoqi,0,"neutral",&(out.nneutralcand),&(out.neutralcandenergies),&(out.neutralcandets),&(out.neutralcanddetas),&(out.neutralcanddphis));
  return out;
};

std::pair<float,float> DiPhotonMiniTree::PFPhotonIsolationFromMinitree(int phoqi1, int phoqi2, pfcandidates_struct *pfcands, bool doremoval1, bool doremoval2, float matched_eta1, float matched_eta2){

  float result1=0;
  float result2=0;
  
  std::set<int> footprint;
  if (doremoval1) {std::set<int> footprint1 = GetPFCandInsideFootprint(fTR,pfcands,phoqi1,0,"photon"); footprint.insert(footprint1.begin(),footprint1.end());} // phoqi<0 ritorna vuoto
  if (doremoval2) {std::set<int> footprint2 = GetPFCandInsideFootprint(fTR,pfcands,phoqi2,0,"photon"); footprint.insert(footprint2.begin(),footprint2.end());} // phoqi<0 ritorna vuoto

  for (int i=0; i<pfcands->PfCandPt.size(); i++){

    bool removed=false;
    for (set<int>::iterator j=footprint.begin(); j!=footprint.end(); j++) {
      if (i==*j) removed=true;
    }
    if (removed) continue;

    if (phoqi1>=0) if (GetPFCandDeltaRFromSC(fTR,pfcands,phoqi1,i,matched_eta1).dR<0.4) result1+=pfcands->PfCandPt[i];
    if (phoqi2>=0) if (GetPFCandDeltaRFromSC(fTR,pfcands,phoqi2,i,matched_eta2).dR<0.4) result2+=pfcands->PfCandPt[i];
    

  } // end pf cand loop


  return std::make_pair<float,float>(float(result1),float(result2));

};

float DiPhotonMiniTree::PFIsolation(int phoqi, float rotation_phi, TString component, int *counter, std::vector<float> *energies, std::vector<float> *ets, std::vector<float> *detas, std::vector<float> *dphis, float* newphi, std::set<int> removals){

  if (component=="combined") return PFIsolation(phoqi,rotation_phi,"photon",counter,energies,ets,detas,dphis,newphi,removals) \
    + PFIsolation(phoqi,rotation_phi,"charged",counter,energies,ets,detas,dphis,newphi,removals) \
    + PFIsolation(phoqi,rotation_phi,"neutral",counter,energies,ets,detas,dphis,newphi,removals);

  if (component!="neutral" && component!="charged" && component!="photon") {std::cout << "wrong" << std::endl; return -999;} 

  float minimal_pfphotoncand_threshold_EB = (component=="photon") ? global_minthrpfphotoncandEB : 0.0;
  float minimal_pfphotoncand_threshold_EE = (component=="photon") ? global_minthrpfphotoncandEE : 0.0;

  // footprint removal only for photons, charged use Poter's veto cones, nothing for neutrals

  float result=0;
  float scaleresult=1;
  
  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);

  bool isbarrel = (fabs(photon_position.Eta())<1.4442);

  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    photon_position *= r;
  } 

  if (newphi) *newphi = photon_position.Phi();

  if (component=="photon"){
    std::set<int> footprint = GetPFCandInsideFootprint(fTR,phoqi,rotation_phi,"photon");
    for (set<int>::iterator i=footprint.begin(); i!=footprint.end(); i++) removals.insert(*i);
    //    scaleresult = scareaSF[fTR->PhotSCindex[phoqi]];
  }

//  if (rotation_phi!=0) {
//    cout << "rotated removals ";
//    for (int i=0; i<removals.size(); i++) cout << fTR->PfCandPt[removals.at(i)] << " ";
//    cout << endl;
//  }

  for (int i=0; i<fTR->NPfCand; i++){

    float pfeta = fTR->PfCandEta[i];
    if (pfeta>1.4442 && pfeta<1.566) continue;
    if (pfeta>2.5) continue;
    if (isbarrel && fabs(pfeta)>1.4442) continue;
    if (!isbarrel && fabs(pfeta)<1.566) continue;

    if (fTR->PhoMatchedPFPhotonCand.size()>0)           if (fTR->PhoMatchedPFPhotonCand[phoqi]==i) continue;
    if (fTR->PhoMatchedPFElectronCand.size()>0)         if (fTR->PhoMatchedPFElectronCand[phoqi]==i) continue;
    if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) if (fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]==i) continue;

    int type = FindPFCandType(fTR->PfCandPdgId[i]);

    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    bool removed = false;
    for (set<int>::iterator j=removals.begin(); j!=removals.end(); j++) {
      if (i==*j) removed=true;
    }
    if (removed) continue;

    double pt = fTR->PfCandPt[i];
    double energy = fTR->PfCandEnergy[i];

    double dR = 999;
    double dEta = 999;
    double dPhi = 999;

    if (type==2){
      angular_distances_struct angles = GetPFCandDeltaRFromSC(fTR,phoqi,i,rotation_phi);
      dR = angles.dR;
      dEta = angles.dEta;
      dPhi = angles.dPhi;
      if (dR>0.4) continue;
    }
    else if (type==0 || type==1) {

      TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
      TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
      TVector3 photon_direction = photon_position-pfvertex;
      float sceta = photon_direction.Eta();
      float scphi = photon_direction.Phi();
      double dxy;
      double dz;

//      if (type==1){
//	
//	TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);
//
//	if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
//	  std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
//	  vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
//	}
//      
//	dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
//	dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
//	dxy=fabs(dxy);
//	dz=fabs(dz);
//
//      }
//

      dEta = fTR->PfCandEta[i] - sceta;
      dPhi = DeltaPhiSigned(fTR->PfCandPhi[i],scphi);
      dR = sqrt(dEta*dEta+dPhi*dPhi);
      
      if (dR>0.4) continue;



//      if (type==2){ 
//	if (fTR->PhoisEB[phoqi]){
//	  if (fabs(dEta)<0.015) continue;
//	}
//	else if (fTR->PhoisEE[phoqi]){
//	  float limit_dR = 0.00864*fabs(sinh(sceta))*4;
//	  if (dR<limit_dR) continue;
//	}
//	else {
//	  std::cout << "Something wrong" << std::endl;
//	  return -999;
//	}
//      }

      if (type==1){
	if (dz>0.2) continue;
	if (dxy>0.1) continue;
	if (dR<0.02) continue;
	//	scaleresult = 1+2.5e-3;
      }

    }
    else {
      std::cout << "something wrong" << std::endl;
      continue;
    }

    // pfcandidate threshold
    if (type==2){ 
      if (fTR->PhoisEB[phoqi]){
	if (pt<minimal_pfphotoncand_threshold_EB) continue;
      }
      else if (fTR->PhoisEE[phoqi]){
	if (pt<minimal_pfphotoncand_threshold_EE) continue;
      }
      else {
	std::cout << "Something wrong" << std::endl;
	return -999;
      }
    }

    if (counter) (*counter)++;
    if (energies) energies->push_back(energy);
    if (ets) ets->push_back(pt);
    if (detas) detas->push_back(dEta);
    if (dphis) dphis->push_back(dPhi);
    result+=pt;

  } // end pf cand loop

  return result*scaleresult-GetPUEnergy(fTR,component,fTR->PhoisEB[phoqi]);

};

float DiPhotonMiniTree::PFPhotonIsolationAroundMuon(int muqi, int *counter, std::vector<float> *energies, std::vector<float> *ets, std::vector<float> *detas, std::vector<float> *dphis){
  
  float minimal_pfphotoncand_threshold_EB = global_minthrpfphotoncandEB;
  float minimal_pfphotoncand_threshold_EE = global_minthrpfphotoncandEE;

  float result=0;
  
  TLorentzVector mu(fTR->MuPx[muqi],fTR->MuPy[muqi],fTR->MuPz[muqi],fTR->MuE[muqi]);
  bool isbarrel = (fabs(mu.Eta())<1.4442);

  for (int i=0; i<fTR->NPfCand; i++){
    
    int type = FindPFCandType(fTR->PfCandPdgId[i]);
    if (type!=2) continue;
    
    double pt = fTR->PfCandPt[i];
    double energy = fTR->PfCandEnergy[i];
    double dEta = fTR->PfCandEta[i] - mu.Eta();
    double dPhi = DeltaPhiSigned(fTR->PfCandPhi[i],mu.Phi());
    double dR = sqrt(dEta*dEta+dPhi*dPhi);
  
    if (dR>0.4) continue;
    if (dR<0.1) continue;
  
  
    // pfcandidate threshold
    if (isbarrel){
      if (pt<minimal_pfphotoncand_threshold_EB) continue;
    }
    else {
      if (pt<minimal_pfphotoncand_threshold_EE) continue;
    }
    
    if (counter) (*counter)++;
    if (energies) energies->push_back(energy);
    if (ets) ets->push_back(pt);
    if (detas) detas->push_back(dEta);
    if (dphis) dphis->push_back(dPhi);

    result+=pt;
    
  } // end pf cand loop

  //  const float scaleresult = 1.066;      
  return result; 
  
};


angular_distances_struct DiPhotonMiniTree::GetPFCandDeltaRFromSC(TreeReader *fTR, int phoqi, int pfindex, float rotation_phi){

  int i = pfindex;
  int scindex = fTR->PhotSCindex[phoqi];

  if (scindex<0) {
    std::cout << "Error in GetPFCandDeltaRFromSC" << std::endl;
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

  bool isbarrel = fTR->PhoisEB[phoqi];

  TVector3 sc_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);

  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    sc_position *= r;
  }

  TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
  TVector3 pfmomentum;
  pfmomentum.SetPtEtaPhi(fTR->PfCandPt[i],fTR->PfCandEta[i],fTR->PfCandPhi[i]);
  pfmomentum = pfmomentum.Unit();

  TVector3 ecalpfhit = PropagatePFCandToEcal(fTR,i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel);

  if ((isbarrel && fabs(fTR->PfCandEta[i])>1.4442) || (!isbarrel && fabs(fTR->PfCandEta[i])<1.566) || (fabs(fTR->PfCandEta[i])>1.4442 && fabs(fTR->PfCandEta[i])<1.566) || (fabs(fTR->PfCandEta[i])>2.5)) {
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

  angular_distances_struct out;
  out.dR = Util::GetDeltaR(sc_position.Eta(),ecalpfhit.Eta(),sc_position.Phi(),ecalpfhit.Phi());
  out.dEta = ecalpfhit.Eta()-sc_position.Eta();
  out.dPhi = DeltaPhiSigned(ecalpfhit.Phi(),sc_position.Phi());

  return out;

}

angular_distances_struct DiPhotonMiniTree::GetPFCandDeltaRFromSC(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, int pfindex, float matched_eta){

  int i = pfindex;
  int scindex = fTR->PhotSCindex[phoqi];

  if (scindex<0) {
    std::cout << "Error in GetPFCandDeltaRFromSC" << std::endl;
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

  bool isbarrel = fTR->PhoisEB[phoqi];

  TVector3 sc_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);

  TVector3 pfvertex(pfcands->PfCandVx[i],pfcands->PfCandVy[i],pfcands->PfCandVz[i]);
  TVector3 pfmomentum; pfmomentum.SetPtEtaPhi(pfcands->PfCandPt[i],pfcands->PfCandEta[i],pfcands->PfCandPhi[i]);
  pfmomentum = pfmomentum.Unit();

  TVector3 ecalpfhit = PropagatePFCandToEcal(pfcands,i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel);

  if (fabs(matched_eta)>2.5) matched_eta=sc_position.Eta(); // necessary protection: if crazy number, don't do any shift correction

  bool bad=false;
  if (isbarrel && (pfcands->PfCandEta[i]>1.4442+sc_position.Eta()-matched_eta || pfcands->PfCandEta[i]<-1.4442+sc_position.Eta()-matched_eta)) bad=true;
  if (!isbarrel && (pfcands->PfCandEta[i]>-1.4442+sc_position.Eta()-matched_eta && pfcands->PfCandEta[i]<1.4442+sc_position.Eta()-matched_eta)) bad=true; 
  if (bad) {
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

  angular_distances_struct out;
  out.dR = Util::GetDeltaR(sc_position.Eta(),ecalpfhit.Eta(),sc_position.Phi(),ecalpfhit.Phi());
  out.dEta = ecalpfhit.Eta()-sc_position.Eta();
  out.dPhi = DeltaPhiSigned(ecalpfhit.Phi(),sc_position.Phi());

  return out;

}

int DiPhotonMiniTree::FindPFCandType(int id){

  int type = -1;

  if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
  if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
  if (id==22) type=2; //photons
  if (fabs(id)==11) type=3; //electrons
  if (fabs(id)==13) type=4; //muons

  return type;
}


void DiPhotonMiniTree::FillLead(int index, std::vector<int> passing_jets){

  pholead_eta = fTR->PhoEta[index];
  pholead_phi = fTR->PhoPhi[index];
  pholead_pt = fTR->PhoPt[index];
  pholead_energy = fTR->PhoEnergy[index];
  pholead_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  pholead_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  pholead_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  pholead_r9 = fTR->PhoR9[index];
  pholead_sieie = fTR->PhoSigmaIetaIeta[index];
  pholead_hoe = fTR->PhoHoverE[index];
  pholead_PhoIso03Ecal = fTR->PhoIso03Ecal[index];
  pholead_PhoIso03Hcal = fTR->PhoIso03Hcal[index];
  pholead_PhoIso03TrkHollow = fTR->PhoIso03TrkHollow[index];

  if (do_recalc_isolation){
    isolations_struct isos = PFConeIsolation(fTR,index);
    pholead_PhoSCRemovalPFIsoCharged=isos.charged;
    pholead_PhoSCRemovalPFIsoNeutral=isos.neutral;
    pholead_PhoSCRemovalPFIsoPhoton=isos.photon;
    if (isos.photon != fTR->PhoSCRemovalPFIsoPhoton[index]){
      cout << "WRONG " << pholead_eta << " " << pholead_SCeta << " " << isos.photon << " " << fTR->PhoSCRemovalPFIsoPhoton[index] << endl;
    }
    pholead_PhoSCRemovalPFIsoCombined=pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
  pholead_Npfcandphotonincone = isos.nphotoncand;
  pholead_Npfcandchargedincone = isos.nchargedcand;
  pholead_Npfcandneutralincone = isos.nneutralcand;
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = isos.photoncandenergies.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i] = isos.chargedcandenergies.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i] = isos.neutralcandenergies.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = isos.photoncandets.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i] = isos.chargedcandets.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i] = isos.neutralcandets.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = isos.photoncanddetas.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i] = isos.chargedcanddetas.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i] = isos.neutralcanddetas.at(i);
  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = isos.photoncanddphis.at(i);
  for (int i=0; i<pholead_Npfcandchargedincone && i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i] = isos.chargedcanddphis.at(i);
  for (int i=0; i<pholead_Npfcandneutralincone && i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i] = isos.neutralcanddphis.at(i);
  }
  else{
    pholead_PhoSCRemovalPFIsoCharged=fTR->PhoSCRemovalPFIsoChargedPrimVtx[index];
    pholead_PhoSCRemovalPFIsoNeutral=fTR->PhoSCRemovalPFIsoNeutral[index];
    pholead_PhoSCRemovalPFIsoPhoton=fTR->PhoSCRemovalPFIsoPhoton[index];
    pholead_PhoSCRemovalPFIsoCombined=pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
  }
  pholead_PhoPassConversionVeto=fTR->PhoPassConversionVeto[index];
  if (!isdata){
    if (fTR->PhoMCmatchindex[index]>=0) pholead_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
    pholead_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
  }
//  pholead_scarea = scarea[fTR->PhotSCindex[index]];
//  pholead_scareaSF = scareaSF[fTR->PhotSCindex[index]];

  {
    jetmatching_struct m = photon_jet_matching.at(index);
    pholead_m_jet_ptcorr = (m.m_jet>=0) ? fTR->JPt[m.m_jet] : -999;
    pholead_m_jet_dR = (m.m_jet>=0) ? Util::GetDeltaR(pholead_eta,fTR->JEta[m.m_jet],pholead_phi,fTR->JPhi[m.m_jet]) : 999;
    pholead_phopt_footprint_total = m.phopt_footprint_total;
    pholead_phopt_footprint_m_frac = m.phopt_footprint_m_frac;
    pholead_jetpt_pf = m.jetpt_pf;
    pholead_jetpt_m_frac = m.jetpt_m_frac;
    pholead_jetpt_m_frac_PhoComp = m.jetpt_m_frac_PhoComp;
    pholead_dR_closestjet = 999;
    pholead_pt_closestjet = 999;
    for (size_t i=0; i<passing_jets.size(); i++){
      if (passing_jets.at(i)==m.m_jet) continue;
      float dR = Util::GetDeltaR(pholead_eta,fTR->JEta[passing_jets.at(i)],pholead_phi,fTR->JPhi[passing_jets.at(i)]);
      if (dR<pholead_dR_closestjet) {pholead_dR_closestjet = dR; pholead_pt_closestjet = fTR->JPt[passing_jets.at(i)];}
    }
  }


//  if (fabs(pholead_phopt_footprint_total/pholead_pt-1)>0.5){
//    int phoqi = index;
//    cout << event_run << " " << event_lumi << " " << event_number << endl;
//    cout << pholead_pt << " " << pholead_SCeta << " " << pholead_SCphi << endl;
//    cout << "qui mancano circa " << pholead_pt - pholead_phopt_footprint_total << endl;
//    cout << pholead_pt << " " << pholead_phopt_footprint_total << endl;
//    cout << pholead_PhoSCRemovalPFIsoCharged << " " << pholead_PhoSCRemovalPFIsoNeutral << " " << pholead_PhoSCRemovalPFIsoPhoton << endl;
//    std::set<int> pfcands = GetPFCandInsideFootprint(fTR,phoqi,0,"photon");
//    if (fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]>=0) {
//      int m = fTR->PhoMatchedPFPhotonOrElectronCand[phoqi];
//      for (set<int>::iterator i=pfcands.begin(); i!=pfcands.end(); i++) if (*i!=m) continue;
//      pfcands.insert(m);
//    }
//    cout << "PF cands in footprint" << endl;
//    float sum=0;
//    for (set<int>::iterator i=pfcands.begin(); i!=pfcands.end(); i++) cout << *i << " " << fTR->PfCandPt[*i] << " " << fTR->PfCandEta[*i] << " " << fTR->PfCandPhi[*i] << " " << fTR->PfCandPdgId[*i] << endl;
//    for (set<int>::iterator i=pfcands.begin(); i!=pfcands.end(); i++) sum+=fTR->PfCandPt[*i];
//    cout << sum << " " << pholead_pt << endl;
//    if (fabs(sum/pholead_pt-1)<0.5) cout << "RECOVERED" << endl; else cout << "WROOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOONG" << endl;
//    cout << "other pfcands nearby" << endl;
//    for (int i=0; i<fTR->NPfCand; i++){
//      bool good = true;
//      for (set<int>::iterator j=pfcands.begin(); j!=pfcands.end(); j++) if (i==*j) good=false;
//      if (!good) continue;
//      if (Util::GetDeltaR(fTR->PfCandEta[i],fTR->PhoEta[phoqi],fTR->PfCandPhi[i],fTR->PhoPhi[phoqi])>0.6) continue;
//      cout << i << " " << fTR->PfCandPt[i] << " " << fTR->PfCandEta[i] << " " << fTR->PfCandPhi[i] << " " << fTR->PfCandPdgId[i] << endl;
//    }
//  }

};

float DiPhotonMiniTree::SieieRescale(float sieie, bool isbarrel){
  if (isdata) return sieie; // rescale sieie only in MC
  return sieie;
  //  SIEIE NOT SCALED NEITHER IN CIC NOR MVA
};

float DiPhotonMiniTree::R9Rescale(float r9, bool isbarrel){
  if (isdata) return r9; // rescale r9 only in MC
  if (global_is2011) return isbarrel ? 1.00153*r9+0.0008543 : 1.0005*r9+0.001231; // 2011 legacy, from globe
  if (global_is2012) return isbarrel ? 1.00793*r9+(-0.00532538) : 1.00017*r9+(-0.0016474); // 2012 legacy, from globe
};

void DiPhotonMiniTree::FillTrail(int index, std::vector<int> passing_jets){

  photrail_eta = fTR->PhoEta[index];
  photrail_phi = fTR->PhoPhi[index];
  photrail_pt = fTR->PhoPt[index];
  photrail_energy = fTR->PhoEnergy[index];
  photrail_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  photrail_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  photrail_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  photrail_r9 = fTR->PhoR9[index];
  photrail_sieie = fTR->PhoSigmaIetaIeta[index];
  photrail_hoe = fTR->PhoHoverE[index];
  photrail_PhoIso03Ecal = fTR->PhoIso03Ecal[index];
  photrail_PhoIso03Hcal = fTR->PhoIso03Hcal[index];
  photrail_PhoIso03TrkHollow = fTR->PhoIso03TrkHollow[index];

  if (do_recalc_isolation){
    isolations_struct isos = PFConeIsolation(fTR,index);
    photrail_PhoSCRemovalPFIsoCharged=isos.charged;
    photrail_PhoSCRemovalPFIsoNeutral=isos.neutral;
    photrail_PhoSCRemovalPFIsoPhoton=isos.photon;
    photrail_PhoSCRemovalPFIsoCombined=photrail_PhoSCRemovalPFIsoCharged+photrail_PhoSCRemovalPFIsoNeutral+photrail_PhoSCRemovalPFIsoPhoton;
  photrail_Npfcandphotonincone = isos.nphotoncand;
  photrail_Npfcandchargedincone = isos.nchargedcand;
  photrail_Npfcandneutralincone = isos.nneutralcand;
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i] = isos.photoncandenergies.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i] = isos.chargedcandenergies.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i] = isos.neutralcandenergies.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i] = isos.photoncandets.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i] = isos.chargedcandets.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i] = isos.neutralcandets.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i] = isos.photoncanddetas.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i] = isos.chargedcanddetas.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i] = isos.neutralcanddetas.at(i);
  for (int i=0; i<photrail_Npfcandphotonincone && i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i] = isos.photoncanddphis.at(i);
  for (int i=0; i<photrail_Npfcandchargedincone && i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i] = isos.chargedcanddphis.at(i);
  for (int i=0; i<photrail_Npfcandneutralincone && i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i] = isos.neutralcanddphis.at(i);
  }
  else{
    photrail_PhoSCRemovalPFIsoCharged=fTR->PhoSCRemovalPFIsoChargedPrimVtx[index];
    photrail_PhoSCRemovalPFIsoNeutral=fTR->PhoSCRemovalPFIsoNeutral[index];
    photrail_PhoSCRemovalPFIsoPhoton=fTR->PhoSCRemovalPFIsoPhoton[index];
    photrail_PhoSCRemovalPFIsoCombined=photrail_PhoSCRemovalPFIsoCharged+photrail_PhoSCRemovalPFIsoNeutral+photrail_PhoSCRemovalPFIsoPhoton;
  }

  photrail_PhoPassConversionVeto=fTR->PhoPassConversionVeto[index];
  if (!isdata){
    if (fTR->PhoMCmatchindex[index]>=0) photrail_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
    photrail_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
  }
//  photrail_scarea = scarea[fTR->PhotSCindex[index]];
//  photrail_scareaSF = scareaSF[fTR->PhotSCindex[index]];

  {
    jetmatching_struct m = photon_jet_matching.at(index);
    photrail_m_jet_ptcorr = (m.m_jet>=0) ? fTR->JPt[m.m_jet] : -999;
    photrail_m_jet_dR = (m.m_jet>=0) ? Util::GetDeltaR(photrail_eta,fTR->JEta[m.m_jet],photrail_phi,fTR->JPhi[m.m_jet]) : 999;
    photrail_phopt_footprint_total = m.phopt_footprint_total;
    photrail_phopt_footprint_m_frac = m.phopt_footprint_m_frac;
    photrail_jetpt_pf = m.jetpt_pf;
    photrail_jetpt_m_frac = m.jetpt_m_frac;
    photrail_jetpt_m_frac_PhoComp = m.jetpt_m_frac_PhoComp;
    photrail_dR_closestjet = 999;
    photrail_pt_closestjet = 999;
    for (size_t i=0; i<passing_jets.size(); i++){
      if (passing_jets.at(i)==m.m_jet) continue;
      float dR = Util::GetDeltaR(photrail_eta,fTR->JEta[passing_jets.at(i)],photrail_phi,fTR->JPhi[passing_jets.at(i)]);
      if (dR<photrail_dR_closestjet) {photrail_dR_closestjet = dR; photrail_pt_closestjet = fTR->JPt[passing_jets.at(i)];}
    }
  }

};

void DiPhotonMiniTree::FillJetsInfo(std::vector<int> passing, std::vector<int> passing_jets){

  n_jets = passing_jets.size();
  for (size_t i=0; i<passing_jets.size(); i++){
    jet_pt[i] =  fTR->JPt[passing_jets.at(i)];
    jet_eta[i] = fTR->JEta[passing_jets.at(i)];
    jet_phi[i] = fTR->JPhi[passing_jets.at(i)];
    jet_energy[i] = fTR->JE[passing_jets.at(i)];
  }

};

void DiPhotonMiniTree::FillGenJetsInfo(std::vector<int> passing_gen, std::vector<int> passing_gen_jets){

  n_GEN_jets = passing_gen_jets.size();
  for (size_t i=0; i<passing_gen_jets.size(); i++){
    jet_GEN_pt[i] =  fTR->GenJetPt[passing_gen_jets.at(i)];
    jet_GEN_eta[i] = fTR->GenJetEta[passing_gen_jets.at(i)];
    jet_GEN_phi[i] = fTR->GenJetPhi[passing_gen_jets.at(i)];
    jet_GEN_energy[i] = fTR->GenJetE[passing_gen_jets.at(i)];
  }

};

void DiPhotonMiniTree::JetSelection(std::vector<int> &passing_jets){

  for (vector<int>::iterator it = passing_jets.begin(); it != passing_jets.end(); ){
    float pt = fTR->JPt[*it];
    float eta = fTR->JEta[*it];
    bool pass=1;
    if (pt<min_jet_pt_alljets) pass=0;
    if (fabs(eta)>max_eta_jets) pass=0;
    if (!(fTR->JEcorr[*it]>0)) pass=0;
    if (!(fTR->JPassPileupIDM0[fTR->JVrtxListStart[*it]+0])) pass=0; // pileup ID medium cut-based
    if (!pass) it=passing_jets.erase(it); else it++;
  }

}

void DiPhotonMiniTree::GenJetSelection(std::vector<int> &passing_gen_jets){

  for (vector<int>::iterator it = passing_gen_jets.begin(); it != passing_gen_jets.end(); ){
    float pt = fTR->GenJetPt[*it];
    float eta = fTR->GenJetEta[*it];
    bool pass=1;
    if (pt<min_jet_pt_alljets) pass=0;
    if (fabs(eta)>max_eta_jets) pass=0;
    if (!pass) it=passing_gen_jets.erase(it); else it++;
  }

}

void DiPhotonMiniTree::FillMuonInfo(int index){
//
//  pholead_eta = fTR->MuEta[index];
//  pholead_phi = fTR->MuPhi[index];
//  pholead_px = fTR->MuPx[index];
//  pholead_py = fTR->MuPy[index];
//  pholead_pt = fTR->MuPt[index];
//  pholead_pz = fTR->MuPz[index];
//  pholead_energy = fTR->MuE[index];
//  pholead_SCeta = fTR->MuEta[index];
//  pholead_SCphi = fTR->MuPhi[index];
//  pholead_pho_Cone03PFCombinedIso = fTR->MuRelIso03[index];
//  pholead_Npfcandphotonincone = 0;
//  std::vector<float> energies;
//  std::vector<float> ets;
//  std::vector<float> detas;
//  std::vector<float> dphis;
//  pholead_PhoSCRemovalPFIsoCharged = PFPhotonIsolationAroundMuon(index,&pholead_Npfcandphotonincone,&energies,&ets,&detas,&dphis);
//  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i] = energies.at(i);
//  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i] = ets.at(i);
//  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i] = detas.at(i);
//  for (int i=0; i<pholead_Npfcandphotonincone && i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i] = dphis.at(i);
};

void DiPhotonMiniTree::ResetVars(){

  event_pass12whoissiglike = -999;
  dipho_mgg_photon = -999;
  pholead_eta = -999;
  photrail_eta = -999;
  pholead_phi = -999;
  photrail_phi = -999;
  pholead_pt = -999;
  photrail_pt = -999;
  pholead_energy = -999;
  photrail_energy = -999;
  pholead_SCeta = -999;
  photrail_SCeta = -999;
  pholead_SCphi = -999;
  photrail_SCphi = -999;
  pholead_PhoHasPixSeed = -999;
  photrail_PhoHasPixSeed = -999;
  pholead_r9 = -999;
  photrail_r9 = -999;
  pholead_sieie = -999;
  photrail_sieie = -999;
  pholead_hoe = -999;
  photrail_hoe = -999;
  pholead_PhoIso03Ecal = -999;
  pholead_PhoIso03Hcal = -999;
  pholead_PhoIso03TrkHollow = -999;
  photrail_PhoIso03Ecal = -999;
  photrail_PhoIso03Hcal = -999;
  photrail_PhoIso03TrkHollow = -999;
  pholead_PhoSCRemovalPFIsoCharged = -999;
  photrail_PhoSCRemovalPFIsoCharged = -999;
  pholead_PhoSCRemovalPFIsoNeutral = -999;
  photrail_PhoSCRemovalPFIsoNeutral = -999;
  pholead_PhoSCRemovalPFIsoPhoton = -999;
  photrail_PhoSCRemovalPFIsoPhoton = -999;
  pholead_PhoSCRemovalPFIsoCombined = -999;
  photrail_PhoSCRemovalPFIsoCombined = -999;
  pholead_PhoPassConversionVeto = -999;
  photrail_PhoPassConversionVeto = -999;
  pholead_GenPhotonIsoDR04 = -999;
  photrail_GenPhotonIsoDR04 = -999;
  pholead_PhoMCmatchexitcode = -999;
  photrail_PhoMCmatchexitcode = -999;
//  pholead_scarea = -999;
//  photrail_scarea = -999;
//  pholead_scareaSF = -999;
//  photrail_scareaSF = -999;
  pholead_m_jet_ptcorr = -999;
  pholead_m_jet_dR = -999;
  pholead_phopt_footprint_total = -999;
  pholead_phopt_footprint_m_frac = -999;
  pholead_jetpt_pf = -999;
  pholead_jetpt_m_frac = -999;
  pholead_jetpt_m_frac_PhoComp = -999;
  photrail_m_jet_ptcorr = -999;
  photrail_m_jet_dR = -999;
  photrail_phopt_footprint_total = -999;
  photrail_phopt_footprint_m_frac = -999;
  photrail_jetpt_pf = -999;
  photrail_jetpt_m_frac = -999;
  photrail_jetpt_m_frac_PhoComp = -999;
  pholead_pt_closestjet = -999;
  pholead_dR_closestjet = -999;
  photrail_pt_closestjet = -999;
  photrail_dR_closestjet = -999;
  if (do_recalc_isolation){
  pholead_Npfcandphotonincone = -999;
  pholead_Npfcandchargedincone = -999;
  pholead_Npfcandneutralincone = -999;
  photrail_Npfcandphotonincone = -999;
  photrail_Npfcandchargedincone = -999;
  photrail_Npfcandneutralincone = -999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_photonpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_chargedpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) pholead_neutralpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcandenergies[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcandets[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcanddetas[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_photonpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_chargedpfcanddphis[i]=-999;
  for (int i=0; i<global_size_pfcandarrays; i++) photrail_neutralpfcanddphis[i]=-999;
  }
//
//  for (int i=0; i<50; i++) pholead_test_rotatedphotoniso[i]=-999;
//  for (int i=0; i<50; i++) pholead_test_rotatedwithcheckphotoniso[i]=-999;
//
  pholead_GEN_eta = -999;
  photrail_GEN_eta = -999;
  pholead_GEN_phi = -999;
  photrail_GEN_phi = -999;
  pholead_GEN_pt = -999;
  photrail_GEN_pt = -999;

  tree_gen_in_acc = false;
  tree_reco_in_acc = false;
  tree_matched = false;

  allphotonpfcand_count = 0;
  for (int i=0; i<global_maxN_photonpfcandidates; i++){
    allphotonpfcand_pt[i]=-999;
    allphotonpfcand_eta[i]=-999;
    allphotonpfcand_phi[i]=-999;
    allphotonpfcand_vx[i]=-999;
    allphotonpfcand_vy[i]=-999;
    allphotonpfcand_vz[i]=-999;
  }

  for (int i=0; i<nclosest; i++){
    phoiso_template_1event_sigsig_1[i]=-999;
    phoiso_template_1event_sigsig_2[i]=-999;
    phoiso_template_1event_sigbkg_1[i]=-999;
    phoiso_template_1event_sigbkg_2[i]=-999;
    phoiso_template_1event_bkgsig_1[i]=-999;
    phoiso_template_1event_bkgsig_2[i]=-999;
    phoiso_template_1event_bkgbkg_1[i]=-999;
    phoiso_template_1event_bkgbkg_2[i]=-999;
    phoiso_template_2events_sigsig_1[i]=-999;
    phoiso_template_2events_sigsig_2[i]=-999;
    phoiso_template_2events_sigbkg_1[i]=-999;
    phoiso_template_2events_sigbkg_2[i]=-999;
    phoiso_template_2events_bkgsig_1[i]=-999;
    phoiso_template_2events_bkgsig_2[i]=-999;
    phoiso_template_2events_bkgbkg_1[i]=-999;
    phoiso_template_2events_bkgbkg_2[i]=-999;
    for (int k=0; k<6; k++){
      rewinfo_template_1event_sigsig_1[i*6+k]=-999;
      rewinfo_template_1event_sigsig_2[i*6+k]=-999;
      rewinfo_template_1event_sigbkg_1[i*6+k]=-999;
      rewinfo_template_1event_sigbkg_2[i*6+k]=-999;
      rewinfo_template_1event_bkgsig_1[i*6+k]=-999;
      rewinfo_template_1event_bkgsig_2[i*6+k]=-999;
      rewinfo_template_1event_bkgbkg_1[i*6+k]=-999;
      rewinfo_template_1event_bkgbkg_2[i*6+k]=-999;
      rewinfo_template_2events_sigsig_1[i*6+k]=-999;
      rewinfo_template_2events_sigsig_2[i*6+k]=-999;
      rewinfo_template_2events_sigbkg_1[i*6+k]=-999;
      rewinfo_template_2events_sigbkg_2[i*6+k]=-999;
      rewinfo_template_2events_bkgsig_1[i*6+k]=-999;
      rewinfo_template_2events_bkgsig_2[i*6+k]=-999;
      rewinfo_template_2events_bkgbkg_1[i*6+k]=-999;
      rewinfo_template_2events_bkgbkg_2[i*6+k]=-999;
    }
  }

  vetoobjects_count=0;
  for (int i=0; i<global_maxN_vetoobjects; i++){
    vetoobjects_pt[i]=-999;
    vetoobjects_eta[i]=-999;
    vetoobjects_phi[i]=-999;
    vetoobjects_type[i]=-999;
  }

  n_jets = -999;
  for (int i=0; i<global_maxN_jets; i++){
    jet_pt[i] = -999;
    jet_eta[i] = -999;
    jet_phi[i] = -999;
    jet_energy[i] = -999;
  }

  n_GEN_jets = -999;
  for (int i=0; i<global_maxN_jets; i++){
    jet_GEN_pt[i] = -999;
    jet_GEN_eta[i] = -999;
    jet_GEN_phi[i] = -999;
    jet_GEN_energy[i] = -999;
  }

};

float DiPhotonMiniTree::CalculateSCArea(TreeReader *fTR, int scindex){
  //  std::cout << "call scarea " << scindex << " " << fTR->SCNXtals[scindex] << std::endl; // DEBUUUUG
  if (scindex>=fTR->NSuperClusters) return -999;
  float area=0;
  for (int i=0; i<fTR->SCNXtals[scindex]; i++) area+=fTR->XtalEtaWidth[fTR->SCXtalListStart[scindex]+i]*fTR->XtalPhiWidth[fTR->SCXtalListStart[scindex]+i];
  return area;
};

float DiPhotonMiniTree::GetPUEnergy(TreeReader *fTR, TString mode, float eta){

  return 0;

  cout << "HAVE YOU INITIALIZED THE EA ARRAYS?" << endl;

  eta=fabs(eta);
  bool isbarrel = (eta<1.4442);

  int bin = Choose_bin_eta(fabs(eta),isbarrel ? 0 : 1);

  float eff_area = 0;

  std::vector<float> *eff_areas_EB_photon = isdata ?  &eff_areas_EB_photon_data :  &eff_areas_EB_photon_MC;
  std::vector<float> *eff_areas_EB_charged = isdata ? &eff_areas_EB_charged_data : &eff_areas_EB_charged_MC;
  std::vector<float> *eff_areas_EB_neutral = isdata ? &eff_areas_EB_neutral_data : &eff_areas_EB_neutral_MC;
  std::vector<float> *eff_areas_EE_photon = isdata ?  &eff_areas_EB_photon_data :  &eff_areas_EB_photon_MC;
  std::vector<float> *eff_areas_EE_charged = isdata ? &eff_areas_EE_charged_data : &eff_areas_EE_charged_MC;
  std::vector<float> *eff_areas_EE_neutral = isdata ? &eff_areas_EE_neutral_data : &eff_areas_EE_neutral_MC;

  if (mode=="photon") eff_area = isbarrel ?  (*eff_areas_EB_photon)[bin] :  (*eff_areas_EE_photon)[bin];
  if (mode=="charged") eff_area = isbarrel ? (*eff_areas_EB_charged)[bin] : (*eff_areas_EE_charged)[bin];
  if (mode=="neutral") eff_area = isbarrel ? (*eff_areas_EB_neutral)[bin] : (*eff_areas_EE_neutral)[bin];

  if (eff_area==0) std::cout << "Warning: problem in eta-dependent EA" << std::endl;

  return TMath::Pi()*0.4*0.4*eff_area*fTR->Rho;

};

Int_t DiPhotonMiniTree::Choose_bin_eta(float eta, int region){

  cout << "HAVE YOU INITIALIZED THE ETA ARRAYS?" <<endl;

  eta=fabs(eta);

  int index;

  std::vector<float> *cuts;

  if (region==0) {cuts=&binsdef_single_gamma_EB_eta; index=n_templates_EB;}
  if (region==1) {cuts=&binsdef_single_gamma_EE_eta; index=n_templates_EE;}

  assert (index!=0);

  (*cuts)[index]=9999;

  if (eta<(*cuts)[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << eta << " cuts[0]=" << (*cuts)[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((eta>=(*cuts)[i]) && (eta<(*cuts)[i+1])) return i;

  std::cout << "WARNING: called bin choice for out-of-range value " << eta << std::endl;
  return -999;

};

std::vector<int> DiPhotonMiniTree::ImpingingTrackSelection(TreeReader *fTR, std::vector<int> passing, bool invert){

  return std::vector<int>();

//
//  for (int i=0; i<100; i++) impinging_track_pfcand[i]=-999;
//
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    
//    bool found=0;
//    int phoqi=*it;
//    
//    found = FindImpingingTrack(fTR,phoqi,impinging_track_pfcand[phoqi]);
//
//    if (invert) { // selection = 0 impinging tracks
//      if (found) it=passing.erase(it); else it++;
//    }
//
//    else { // selection = 1 impinging track
//
//      if (found){    
//	std::vector<int> remove;
//	remove.push_back(impinging_track_pfcand[phoqi]);
//	float eta=fTR->SCEta[fTR->PhotSCindex[*it]];
//	const float eff_area = (fabs(eta)<1.4442) ? eff_area_EB : eff_area_EE;
//	const float dR=0.4;
//	float puenergy =3.14*dR*dR*eff_area*fTR->Rho;
//	if (PFIsolation(phoqi,0,"combined",remove)-puenergy>5) found=0;
//      }
//      
//      if (!found) it=passing.erase(it); else it++;
//    
//    }
//
//  }
//
//  return passing;
//  
};


bool DiPhotonMiniTree::FindImpingingTrack(TreeReader *fTR, int phoqi, int &reference_index_found, bool dofootprintremoval, std::set<int> removals){

  return false;

  // WARNING THE FOLLOWING IS NOT UP TO DATE!!!
//  bool found = false;
//
//  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);
//
//  if (dofootprintremoval){
//    std::vector<int> footprint = GetPFCandInsideFootprint(fTR,phoqi,0,"charged");
//    for (int i=0; i<footprint.size(); i++) removals.push_back(footprint.at(i));
//  }
//
//  for (int i=0; i<fTR->NPfCand; i++){
//
//    if (fTR->PhoisPFPhoton[phoqi] && fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]==i) continue;
//
//    int type = FindPFCandType(fTR->PfCandPdgId[i]);
//    if (type!=1) continue;
//
//    bool removed = false;
//    for (int j=0; j<removals.size(); j++) {
//      if (i==removals.at(j)) removed=true;
//    }
//    if (removed) continue;
//
//    double pt = fTR->PfCandPt[i];
//    if (pt<1.5) continue;
//
//    if (!(fTR->PfCandHasHitInFirstPixelLayer[i])) continue;
//
//    if (dofootprintremoval){
//      
//      std::cout << "Wrong! No propagator for charged stuff" << std::endl;
//      continue;
//
//      //      double dR = GetPFCandDeltaRFromSC(fTR,phoqi,i,0);
//      //      if (dR>0.4) continue;
//    }
//
//    if (!dofootprintremoval){ // veto cones
//
//      TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
//      TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
//      TVector3 photon_direction = photon_position-pfvertex;
//      float sceta = photon_direction.Eta();
//      float scphi = photon_direction.Phi();
//      double dxy;
//      double dz;
//	
//      TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);
//
//      if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
//	std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
//	vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
//      }
//      
//      dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
//      dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
//      dxy=fabs(dxy);
//      dz=fabs(dz);
//
//
//      double dEta = fTR->PfCandEta[i] - sceta;
//      double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);
//      double dR = sqrt(dEta*dEta+dPhi*dPhi);
//
//      if (dR>0.4) continue;
//
//      // additional veto sieie-orthogonal
//      if (fTR->PhoisEB[phoqi]){
//	if (fabs(dEta)<0.05) continue;
//	if (fabs(dPhi)<0.05) continue;
//      }
//
//      if (dz>0.2) continue;
//      if (dxy>0.1) continue;
//      if (dR<0.02) continue;
//
//    }
//
//    found=1;
//    reference_index_found = i;
//    break;
//
//  } // end pf cand loop
//
//  return found;
//
};

int DiPhotonMiniTree::CountChargedHadronsInCone(TreeReader *fTR, int phoqi, std::set<int> removals, bool skipvetocones){

  return 0;

//  int found = 0;
//  
//  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);
//  TVector3 phovtx(fTR->PhoVx[phoqi],fTR->PhoVy[phoqi],fTR->PhoVz[phoqi]);
//    
//  for (int i=0; i<fTR->NPfCand; i++){
//    
//    bool removed = false;
//    for (int j=0; j<removals.size(); j++) {
//      if (i==removals.at(j)) removed=true;
//    }
//    if (removed) continue;
//    
//    float id = fTR->PfCandPdgId[i];
//    if (!(fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212)) continue;
//    
//    TVector3 pfvertex(fTR->PfCandVx[i],fTR->PfCandVy[i],fTR->PfCandVz[i]);
//    TVector3 photon_direction = photon_position-pfvertex;
//      
//    double sceta = photon_direction.Eta();
//    double scphi = photon_direction.Phi();
//
//    double dEta = fTR->PfCandEta[i] - sceta;
//    double dPhi = Util::DeltaPhi(fTR->PfCandPhi[i],scphi);
//    double dR = sqrt(dEta*dEta+dPhi*dPhi);
//
//    double pt = fTR->PfCandPt[i];
//
//    TVector3 vtxmom(fTR->PfCandTrackRefPx[i],fTR->PfCandTrackRefPy[i],fTR->PfCandTrackRefPz[i]);
//
//    if (vtxmom.x()==-999 || vtxmom.y()==-999 || vtxmom.z()==-999) {
//      std::cout << "Something wrong with vtxmom from trackref, fallback" << std::endl;
//      vtxmom = TVector3(fTR->PfCandPx[i],fTR->PfCandPy[i],fTR->PfCandPz[i]);
//    }
//      
//    double dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
//    double dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
//    dxy=fabs(dxy);
//    dz=fabs(dz);
//
//    if (dR>0.4) continue;
//
//    if (!skipvetocones){
//      // additional veto sieie-orthogonal
//      if (fabs(sceta)<1.4442){
//	if (fabs(dEta)<0.05) continue;
//	if (fabs(dPhi)<0.05) continue;
//      }
//      
//      if (dz>0.2) continue;
//      if (dxy>0.1) continue;
//      if (dR<0.02) continue;
//    }
//
//    if (!(fTR->PfCandHasHitInFirstPixelLayer[i])) continue;
//
//    found++;
//      
//  } // end pf cand loop
//
//  return found;
//
};

std::set<int> DiPhotonMiniTree::NChargedHadronsInConeSelection(TreeReader *fTR, std::vector<int> passing, int minimum, int maximum){

  return std::set<int>();

//  if (maximum<minimum) {
//    std::cout << "Wrong values for NChargedHadronsInCone selection!!!" << std::endl;
//    return std::vector<int>();
//  }
//
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
//    bool pass=0;
//    std::vector<int> remove = GetPFCandRemovals(fTR,*it);
//    std::vector<int> footprint = GetPFCandInFootprint(fTR,*it,-999);
//    for (int i=0; i<footprint.size(); i++) remove.push_back(footprint.at(i));
//    int n = CountChargedHadronsInCone(fTR,*it,remove,trueblabla);
//    if (n>=minimum && n<=maximum) pass=1; else pass=0;
//    if (!pass) it=passing.erase(it); else it++;
//  }
//
//  return passing;
//
};

void DiPhotonMiniTree::Fillhist_PFPhotonDepositAroundImpingingTrack(int phoqi, int trkindex){
  return;
};

float DiPhotonMiniTree::DeltaPhiSigned(float phi1, float phi2){
  //copy-paste from reco::deltaPhi
  double result = phi1 - phi2;
  const float pi = TMath::Pi();
  while (result > pi) result -= 2*pi;
  while (result <= -pi) result += 2*pi;
  return result;
};

bool DiPhotonMiniTree::PassPrimaryVertexFilter(){
  for (int i=0; i<fTR->NVrtx; i++){
    if (fTR->VrtxNdof[i]>4 && fabs(fTR->VrtxZ[i])<=24 && TVector3(fTR->VrtxX[i],fTR->VrtxY[i],fTR->VrtxZ[i]).Perp()<=2) return true;
  }
  return false;
};


void DiPhotonMiniTree::InitInputTree(){
  if (isstep2){
    TString title = Form("%s/matchingtree_%u.root",input_filename.Data(),uuid);
    cout << "opening " << title.Data() << endl;
    f_input = TFile::Open(title.Data(),"read");
    if (f_input->IsZombie()) cerr << "WRONG: impossible to open file " << title.Data() << endl;
    f_input->GetObject("Tree_1Drandomcone_template_EXTRA",InputTree[0]);
    f_input->GetObject("Tree_2Drandomconesideband_template_EXTRA",InputTree[1]);
    for (int m=0; m<2; m++){
      InputTree[m]->SetBranchAddress("event_number", &input_event_number, &b_input_event_number);
      InputTree[m]->SetBranchAddress("allphotonpfcand_count", &input_allphotonpfcand_count, &b_input_allphotonpfcand_count);
      InputTree[m]->SetBranchAddress("allphotonpfcand_pt", input_allphotonpfcand_pt   , &b_input_allphotonpfcand_pt   );
      InputTree[m]->SetBranchAddress("allphotonpfcand_eta", input_allphotonpfcand_eta  , &b_input_allphotonpfcand_eta  );
      InputTree[m]->SetBranchAddress("allphotonpfcand_phi", input_allphotonpfcand_phi  , &b_input_allphotonpfcand_phi  );
      InputTree[m]->SetBranchAddress("allphotonpfcand_vx", input_allphotonpfcand_vx   , &b_input_allphotonpfcand_vx   );
      InputTree[m]->SetBranchAddress("allphotonpfcand_vy", input_allphotonpfcand_vy   , &b_input_allphotonpfcand_vy   );
      InputTree[m]->SetBranchAddress("allphotonpfcand_vz", input_allphotonpfcand_vz   , &b_input_allphotonpfcand_vz   );
      InputTree[m]->SetBranchAddress("pholead_SCeta", &input_pholead_SCeta, &b_input_pholead_SCeta);
      InputTree[m]->SetBranchAddress("pholead_SCphi", &input_pholead_SCphi, &b_input_pholead_SCphi);
      InputTree[m]->SetBranchAddress("pholead_pt", &input_pholead_pt, &b_input_pholead_pt);
      InputTree[m]->SetBranchAddress("photrail_SCeta", &input_photrail_SCeta, &b_input_photrail_SCeta);
      InputTree[m]->SetBranchAddress("photrail_SCphi", &input_photrail_SCphi, &b_input_photrail_SCphi);
      InputTree[m]->SetBranchAddress("photrail_pt", &input_photrail_pt, &b_input_photrail_pt);
      InputTree[m]->SetBranchAddress("event_rho", &input_event_rho, &b_input_event_rho);
      InputTree[m]->SetBranchAddress("event_sigma", &input_event_sigma, &b_input_event_sigma);
      InputTree[m]->SetBranchAddress("event_pass12whoissiglike", &input_event_pass12whoissiglike, &b_input_event_pass12whoissiglike);
      InputTree[m]->SetBranchAddress("vetoobjects_count",&input_vetoobjects_count,&b_input_vetoobjects_count);
      InputTree[m]->SetBranchAddress("vetoobjects_eta", input_vetoobjects_eta, &b_input_vetoobjects_eta);
      InputTree[m]->SetBranchAddress("vetoobjects_phi", input_vetoobjects_phi, &b_input_vetoobjects_phi);
    }

    f_input->GetObject("matchingtree",matchingtree);
    matchingtree->SetBranchAddress("matchingtree_event_run",&matchingtree_event_run,&b_matchingtree_event_run);
    matchingtree->SetBranchAddress("matchingtree_event_lumi",&matchingtree_event_lumi,&b_matchingtree_event_lumi);
    matchingtree->SetBranchAddress("matchingtree_event_number",&matchingtree_event_number,&b_matchingtree_event_number);
    matchingtree->SetBranchAddress("matchingtree_index_1event_sigsig_1",matchingtree_index_1event_sigsig_1,&b_matchingtree_index_1event_sigsig_1);
    matchingtree->SetBranchAddress("matchingtree_index_1event_sigsig_2",matchingtree_index_1event_sigsig_2,&b_matchingtree_index_1event_sigsig_2);
    matchingtree->SetBranchAddress("matchingtree_index_1event_sigbkg_1",matchingtree_index_1event_sigbkg_1,&b_matchingtree_index_1event_sigbkg_1);
    matchingtree->SetBranchAddress("matchingtree_index_1event_sigbkg_2",matchingtree_index_1event_sigbkg_2,&b_matchingtree_index_1event_sigbkg_2);
    matchingtree->SetBranchAddress("matchingtree_index_1event_bkgsig_1",matchingtree_index_1event_bkgsig_1,&b_matchingtree_index_1event_bkgsig_1);
    matchingtree->SetBranchAddress("matchingtree_index_1event_bkgsig_2",matchingtree_index_1event_bkgsig_2,&b_matchingtree_index_1event_bkgsig_2);
    matchingtree->SetBranchAddress("matchingtree_index_1event_bkgbkg_1",matchingtree_index_1event_bkgbkg_1,&b_matchingtree_index_1event_bkgbkg_1);
    matchingtree->SetBranchAddress("matchingtree_index_1event_bkgbkg_2",matchingtree_index_1event_bkgbkg_2,&b_matchingtree_index_1event_bkgbkg_2);
    matchingtree->SetBranchAddress("matchingtree_index_2events_sigsig_1",matchingtree_index_2events_sigsig_1,&b_matchingtree_index_2events_sigsig_1);
    matchingtree->SetBranchAddress("matchingtree_index_2events_sigsig_2",matchingtree_index_2events_sigsig_2,&b_matchingtree_index_2events_sigsig_2);
    matchingtree->SetBranchAddress("matchingtree_index_2events_sigbkg_1",matchingtree_index_2events_sigbkg_1,&b_matchingtree_index_2events_sigbkg_1);
    matchingtree->SetBranchAddress("matchingtree_index_2events_sigbkg_2",matchingtree_index_2events_sigbkg_2,&b_matchingtree_index_2events_sigbkg_2);
    matchingtree->SetBranchAddress("matchingtree_index_2events_bkgsig_1",matchingtree_index_2events_bkgsig_1,&b_matchingtree_index_2events_bkgsig_1);
    matchingtree->SetBranchAddress("matchingtree_index_2events_bkgsig_2",matchingtree_index_2events_bkgsig_2,&b_matchingtree_index_2events_bkgsig_2);
    matchingtree->SetBranchAddress("matchingtree_index_2events_bkgbkg_1",matchingtree_index_2events_bkgbkg_1,&b_matchingtree_index_2events_bkgbkg_1);
    matchingtree->SetBranchAddress("matchingtree_index_2events_bkgbkg_2",matchingtree_index_2events_bkgbkg_2,&b_matchingtree_index_2events_bkgbkg_2);

    matchingtree->BuildIndex("matchingtree_event_number%2","matchingtree_event_number>>1");

    inputtree_isinitialized = true;

  }

};


jetmatching_struct DiPhotonMiniTree::PFMatchPhotonToJet(int phoqi){ // returns (jet_index,fraction)

  jetmatching_struct out;
  out.m_jet=-999;
  out.phopt_footprint_total=-999;
  out.phopt_footprint_m_frac=-999;
  out.jetpt_pf=-999;
  out.jetpt_m_frac=-999;

  if (fTR->PhotSCindex[phoqi]<0) return out;
  if (fTR->PhoPt[phoqi]<25) return out;

  // prepare list of pfcands to represent the photon deposit
  std::set<int> pfcands = GetPrecalculatedFootprintPhoEl(phoqi);
  if (fTR->PhoMatchedPFPhotonCand.size()>0)           if (fTR->PhoMatchedPFPhotonCand[phoqi]>=0) pfcands.insert(fTR->PhoMatchedPFPhotonCand[phoqi]);
  if (fTR->PhoMatchedPFElectronCand.size()>0)         if (fTR->PhoMatchedPFElectronCand[phoqi]>=0) pfcands.insert(fTR->PhoMatchedPFElectronCand[phoqi]);
  if (fTR->PhoMatchedPFPhotonOrElectronCand.size()>0) if (fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]>=0) pfcands.insert(fTR->PhoMatchedPFPhotonOrElectronCand[phoqi]);

  // init ranking
  std::vector<std::pair<int,float> > ranking;
  for (int i=0; i<fTR->NJets; i++) ranking.push_back(std::pair<int,float>(i,0));
  if (ranking.size()==0) {
    //    cout << "PFMatchPhotonToJet: no jets in the event! Returning error state" << endl;
    return out;
  }

  // loop on candidates
  out.phopt_footprint_total=0;
  for (set<int>::iterator i=pfcands.begin(); i!=pfcands.end(); i++){
    float ecalfrac = fTR->PfCandEcalEnergy.at(*i)/fTR->PfCandEnergy.at(*i);
    float pt = fTR->PfCandPt.at(*i)*ecalfrac;
    out.phopt_footprint_total+=pt;
    int j = fTR->PfCandBelongsToJet.at(*i);
    if (j<0) continue; // unclustered or belonging to jet not present in ntuple
    ranking.at(j).second+=pt;
  }

  // sort sharings
  std::sort(ranking.begin(),ranking.end(),indexComparator); // sort done w.r.t. absolute pt sharing
  std::pair<int,float> chosen = ranking.at(0);

  if (!(chosen.second>0)) return out;

  // output
  out.m_jet = chosen.first;
  out.phopt_footprint_m_frac = chosen.second/out.phopt_footprint_total;
  out.jetpt_pf = fTR->JPt.at(chosen.first)/fTR->JEcorr.at(chosen.first);
  out.jetpt_m_frac = chosen.second/out.jetpt_pf;
  out.jetpt_m_frac_PhoComp = out.jetpt_m_frac/fTR->JPhoFrac.at(chosen.first);

  return out;

};

float DiPhotonMiniTree::EnergyScaleOffset(float eta, float r9, int run, float sigmas){
  assert(isdata);
  assert(energyScaleDatabase.size()>0);
  std::vector<struct_escale_item>::const_iterator it;
  for (it = energyScaleDatabase.begin(); it!=energyScaleDatabase.end(); it++){
    if (run<it->mrun) continue;
    if (run>it->Mrun) continue;
    if (fabs(eta)<it->meta) continue;
    if (fabs(eta)>=it->Meta) continue;
    if (r9<it->mr9) continue;
    if (r9>=it->Mr9) continue;
    break;
  }
  if (it==energyScaleDatabase.end()) {cout << "WRONG: energy scale item not found " << run << " " << fabs(eta) << " " << r9 << endl; return 1;}
  return it->val+sigmas*it->err;
};

float DiPhotonMiniTree::EnergySmearingCorrection(float eta, float r9, int run, float sigmas){
  assert(!isdata);
  assert(energySmearingDatabase.size()>0);
  std::vector<struct_escale_item>::const_iterator it;
  for (it = energySmearingDatabase.begin(); it!=energySmearingDatabase.end(); it++){
    if (run<it->mrun) continue;
    if (run>it->Mrun) continue;
    if (fabs(eta)<it->meta) continue;
    if (fabs(eta)>=it->Meta) continue;
    if (r9<it->mr9) continue;
    if (r9>=it->Mr9) continue;
    break;
  }
  if (it==energySmearingDatabase.end()) {cout << "WRONG: energy smearing item not found " << run << " " << fabs(eta) << " " << r9 << endl; return 1;}
  return it->val+sigmas*it->err;
};

void DiPhotonMiniTree::InitEnergyScalesAndSmearingsDatabase(){

  assert (energyScaleDatabase.size()==0);
  assert (energySmearingDatabase.size()==0);

  if (year==2011){
  // https://twiki.cern.ch/twiki/pub/CMS/ECALELF/21Jun2012_7TeV-step2-invMass_SC_regrCorrSemiPar7TeVtrainV8_pho-loose-Et_25-noPF-HggRunEtaR9.dat
  InsertEnergyScaleItem(0,1,-999,0.94,160431,165547,0.9961,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,165548,167042,0.9971,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,167043,172400,0.9974,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,172401,173663,0.9976,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,173664,176840,0.9964,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,176841,177775,0.9971,0.0001);
  InsertEnergyScaleItem(0,1,-999,0.94,177776,178723,0.9962,0.0002);
  InsertEnergyScaleItem(0,1,-999,0.94,178724,180252,0.9961,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,160431,165547,0.9929,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,165548,167042,0.9939,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,167043,172400,0.9942,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,172401,173663,0.9944,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,173664,176840,0.9932,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,176841,177775,0.9939,0.0001);
  InsertEnergyScaleItem(0,1,0.94,999,177776,178723,0.9929,0.0002);
  InsertEnergyScaleItem(0,1,0.94,999,178724,180252,0.9929,0.0002);
  InsertEnergyScaleItem(1,1.5,-999,0.94,160431,165547,0.9981,0.0008);
  InsertEnergyScaleItem(1,1.5,-999,0.94,165548,167042,0.9973,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,167043,172400,0.9981,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,172401,173663,0.9980,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,173664,176840,0.9952,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,176841,177775,0.9965,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,177776,178723,0.9964,0.0006);
  InsertEnergyScaleItem(1,1.5,-999,0.94,178724,180252,0.9965,0.0006);
  InsertEnergyScaleItem(1,1.5,0.94,999,160431,165547,0.9883,0.0009);
  InsertEnergyScaleItem(1,1.5,0.94,999,165548,167042,0.9876,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,167043,172400,0.9883,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,172401,173663,0.9883,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,173664,176840,0.9854,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,176841,177775,0.9868,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,177776,178723,0.9867,0.0007);
  InsertEnergyScaleItem(1,1.5,0.94,999,178724,180252,0.9868,0.0007);
  InsertEnergyScaleItem(1.5,2,-999,0.94,160431,165547,0.9961,0.0011);
  InsertEnergyScaleItem(1.5,2,-999,0.94,165548,167042,0.9973,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,167043,172400,0.9973,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,172401,173663,0.9973,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,173664,176840,0.9985,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,176841,177775,0.9975,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,177776,178723,0.9969,0.0008);
  InsertEnergyScaleItem(1.5,2,-999,0.94,178724,180252,0.9962,0.0008);
  InsertEnergyScaleItem(1.5,2,0.94,999,160431,165547,0.9944,0.0012);
  InsertEnergyScaleItem(1.5,2,0.94,999,165548,167042,0.9956,0.0009);
  InsertEnergyScaleItem(1.5,2,0.94,999,167043,172400,0.9956,0.0009);
  InsertEnergyScaleItem(1.5,2,0.94,999,172401,173663,0.9956,0.0009);
  InsertEnergyScaleItem(1.5,2,0.94,999,173664,176840,0.9968,0.0009);
  InsertEnergyScaleItem(1.5,2,0.94,999,176841,177775,0.9958,0.0008);
  InsertEnergyScaleItem(1.5,2,0.94,999,177776,178723,0.9952,0.0009);
  InsertEnergyScaleItem(1.5,2,0.94,999,178724,180252,0.9945,0.0009);
  InsertEnergyScaleItem(2,3,-999,0.94,160431,165547,0.9987,0.0011);
  InsertEnergyScaleItem(2,3,-999,0.94,165548,167042,0.9982,0.0008);
  InsertEnergyScaleItem(2,3,-999,0.94,167043,172400,0.9993,0.0007);
  InsertEnergyScaleItem(2,3,-999,0.94,172401,173663,0.9993,0.0008);
  InsertEnergyScaleItem(2,3,-999,0.94,173664,176840,0.9983,0.0007);
  InsertEnergyScaleItem(2,3,-999,0.94,176841,177775,0.9994,0.0007);
  InsertEnergyScaleItem(2,3,-999,0.94,177776,178723,1.0006,0.0008);
  InsertEnergyScaleItem(2,3,-999,0.94,178724,180252,0.9995,0.0008);
  InsertEnergyScaleItem(2,3,0.94,999,160431,165547,0.9941,0.0011);
  InsertEnergyScaleItem(2,3,0.94,999,165548,167042,0.9937,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,167043,172400,0.9948,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,172401,173663,0.9948,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,173664,176840,0.9938,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,176841,177775,0.9949,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,177776,178723,0.9961,0.0007);
  InsertEnergyScaleItem(2,3,0.94,999,178724,180252,0.9950,0.0008);

  // https://twiki.cern.ch/twiki/pub/CMS/ECALELF/21Jun2012_7TeV-outProfile-scaleStep2smearing-Et_25-noPF-FitResult.config
  InsertEnergySmearingItem(0,1,-999,0.94,0,999999,0.96,0.03);
  InsertEnergySmearingItem(0,1,0.94,999,0,999999,0.68,0.04);
  InsertEnergySmearingItem(1,1.5,-999,0.94,0,999999,1.85,0.04);
  InsertEnergySmearingItem(1,1.5,0.94,999,0,999999,1.01,0.14);
  InsertEnergySmearingItem(1.5,2,-999,0.94,0,999999,1.85,0.07);
  InsertEnergySmearingItem(1.5,2,0.94,999,0,999999,1.58,0.18);
  InsertEnergySmearingItem(2,3,-999,0.94,0,999999,1.83,0.09);
  InsertEnergySmearingItem(2,3,0.94,999,0,999999,2.01,0.06);

  }
  if (year==2012){
    // https://twiki.cern.ch/twiki/pub/CMS/ECALELF/22Jan2012-runDepPowheg-noR9shift-step2-invMass_SC_regrCorrSemiParV5_pho-loose-Et_25-trigger-noPF-HggRunEtaR9.dat
    InsertEnergyScaleItem(0,1,-999,0.94,190645,190781,0.9922,0.0004);
    InsertEnergyScaleItem(0,1,-999,0.94,190782,191042,0.9989,0.0004);
    InsertEnergyScaleItem(0,1,-999,0.94,191043,191720,0.9931,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,191721,193833,0.9922,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,193834,194116,0.9929,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,194117,194427,0.9935,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,194428,194618,0.9929,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,194619,194789,0.9932,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,194790,195111,0.9938,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,195112,195377,0.9940,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,195378,195398,0.9931,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,195399,195657,0.9936,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,195658,195918,0.9942,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,195919,196198,0.9936,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,196199,196356,0.9943,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,196357,198115,0.9938,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,198116,198940,0.9934,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,198941,199317,0.9936,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,199318,199428,0.9933,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,199429,199697,0.9935,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,199698,199832,0.9938,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,199833,199960,0.9940,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,199961,200151,0.9942,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,200152,200490,0.9940,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,200491,200991,0.9947,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,200992,201201,0.9937,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,201202,201624,0.9943,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,201625,201707,0.9945,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,201708,202059,0.9944,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,202060,202204,0.9947,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,202205,202332,0.9951,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,202333,202972,0.9949,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,202973,203002,0.9944,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,203003,203852,0.9958,0.0006);
    InsertEnergyScaleItem(0,1,-999,0.94,203853,204099,0.9935,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,204100,204562,0.9939,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,204563,205085,0.9938,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,205086,205310,0.9939,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,205311,205617,0.9938,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,205618,205825,0.9942,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,205826,206207,0.9949,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,206208,206389,0.9946,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,206390,206483,0.9945,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,206484,206597,0.9945,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,206598,206896,0.9941,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,206897,207220,0.9952,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,207221,207315,0.9949,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,207316,207489,0.9950,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,207490,207919,0.9951,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,207920,208351,0.9947,0.0002);
    InsertEnergyScaleItem(0,1,-999,0.94,208352,208686,0.9952,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,190645,190781,0.9894,0.0004);
    InsertEnergyScaleItem(0,1,0.94,+999,190782,191042,0.9961,0.0004);
    InsertEnergyScaleItem(0,1,0.94,+999,191043,191720,0.9902,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,191721,193833,0.9894,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,193834,194116,0.9900,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,194117,194427,0.9907,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,194428,194618,0.9901,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,194619,194789,0.9904,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,194790,195111,0.9910,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,195112,195377,0.9912,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,195378,195398,0.9903,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,195399,195657,0.9908,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,195658,195918,0.9914,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,195919,196198,0.9908,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,196199,196356,0.9915,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,196357,198115,0.9910,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,198116,198940,0.9906,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,198941,199317,0.9907,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,199318,199428,0.9904,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,199429,199697,0.9907,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,199698,199832,0.9910,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,199833,199960,0.9912,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,199961,200151,0.9913,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,200152,200490,0.9912,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,200491,200991,0.9919,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,200992,201201,0.9909,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,201202,201624,0.9915,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,201625,201707,0.9917,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,201708,202059,0.9916,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,202060,202204,0.9919,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,202205,202332,0.9923,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,202333,202972,0.9921,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,202973,203002,0.9916,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,203003,203852,0.9930,0.0006);
    InsertEnergyScaleItem(0,1,0.94,+999,203853,204099,0.9906,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,204100,204562,0.9910,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,204563,205085,0.9910,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,205086,205310,0.9911,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,205311,205617,0.9909,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,205618,205825,0.9914,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,205826,206207,0.9921,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,206208,206389,0.9918,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,206390,206483,0.9917,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,206484,206597,0.9917,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,206598,206896,0.9913,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,206897,207220,0.9923,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,207221,207315,0.9921,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,207316,207489,0.9922,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,207490,207919,0.9922,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,207920,208351,0.9919,0.0002);
    InsertEnergyScaleItem(0,1,0.94,+999,208352,208686,0.9924,0.0002);
    InsertEnergyScaleItem(1,1.5,-999,0.94,190645,190781,0.9982,0.0016);
    InsertEnergyScaleItem(1,1.5,-999,0.94,190782,191042,1.0014,0.0017);
    InsertEnergyScaleItem(1,1.5,-999,0.94,191043,191720,0.9963,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,191721,193833,0.9982,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,193834,194116,0.9970,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,194117,194427,0.9975,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,194428,194618,0.9973,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,194619,194789,0.9979,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,194790,195111,0.9992,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,195112,195377,0.9976,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,195378,195398,0.9968,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,195399,195657,0.9993,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,195658,195918,0.9983,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,195919,196198,0.9980,0.0009);
    InsertEnergyScaleItem(1,1.5,-999,0.94,196199,196356,0.9983,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,196357,198115,0.9977,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,198116,198940,0.9970,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,198941,199317,0.9970,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,199318,199428,0.9973,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,199429,199697,0.9976,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,199698,199832,0.9985,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,199833,199960,0.9983,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,199961,200151,0.9970,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,200152,200490,0.9984,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,200491,200991,0.9971,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,200992,201201,0.9969,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,201202,201624,0.9992,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,201625,201707,0.9983,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,201708,202059,0.9974,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,202060,202204,0.9977,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,202205,202332,0.9989,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,202333,202972,0.9994,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,202973,203002,0.9962,0.0010);
    InsertEnergyScaleItem(1,1.5,-999,0.94,203003,203852,0.9936,0.0025);
    InsertEnergyScaleItem(1,1.5,-999,0.94,203853,204099,0.9991,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,204100,204562,0.9999,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,204563,205085,1.0006,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,205086,205310,0.9991,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,205311,205617,0.9995,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,205618,205825,0.9982,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,205826,206207,1.0012,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,206208,206389,0.9995,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,206390,206483,0.9995,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,206484,206597,0.9997,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,206598,206896,0.9998,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,206897,207220,0.9997,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,207221,207315,1.0006,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,207316,207489,1.0001,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,207490,207919,1.0006,0.0007);
    InsertEnergyScaleItem(1,1.5,-999,0.94,207920,208351,1.0000,0.0008);
    InsertEnergyScaleItem(1,1.5,-999,0.94,208352,208686,1.0011,0.0007);
    InsertEnergyScaleItem(1,1.5,0.94,+999,190645,190781,0.9876,0.0016);
    InsertEnergyScaleItem(1,1.5,0.94,+999,190782,191042,0.9909,0.0017);
    InsertEnergyScaleItem(1,1.5,0.94,+999,191043,191720,0.9858,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,191721,193833,0.9876,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,193834,194116,0.9865,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,194117,194427,0.9870,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,194428,194618,0.9867,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,194619,194789,0.9874,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,194790,195111,0.9887,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,195112,195377,0.9871,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,195378,195398,0.9863,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,195399,195657,0.9888,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,195658,195918,0.9878,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,195919,196198,0.9874,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,196199,196356,0.9877,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,196357,198115,0.9871,0.0007);
    InsertEnergyScaleItem(1,1.5,0.94,+999,198116,198940,0.9864,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,198941,199317,0.9865,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,199318,199428,0.9868,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,199429,199697,0.9870,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,199698,199832,0.9879,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,199833,199960,0.9877,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,199961,200151,0.9865,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,200152,200490,0.9879,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,200491,200991,0.9866,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,200992,201201,0.9864,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,201202,201624,0.9886,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,201625,201707,0.9878,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,201708,202059,0.9869,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,202060,202204,0.9871,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,202205,202332,0.9883,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,202333,202972,0.9888,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,202973,203002,0.9857,0.0010);
    InsertEnergyScaleItem(1,1.5,0.94,+999,203003,203852,0.9831,0.0025);
    InsertEnergyScaleItem(1,1.5,0.94,+999,203853,204099,0.9886,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,204100,204562,0.9894,0.0009);
    InsertEnergyScaleItem(1,1.5,0.94,+999,204563,205085,0.9901,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,205086,205310,0.9886,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,205311,205617,0.9890,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,205618,205825,0.9877,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,205826,206207,0.9908,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,206208,206389,0.9890,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,206390,206483,0.9890,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,206484,206597,0.9892,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,206598,206896,0.9893,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,206897,207220,0.9892,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,207221,207315,0.9901,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,207316,207489,0.9896,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,207490,207919,0.9901,0.0007);
    InsertEnergyScaleItem(1,1.5,0.94,+999,207920,208351,0.9895,0.0008);
    InsertEnergyScaleItem(1,1.5,0.94,+999,208352,208686,0.9906,0.0007);
    InsertEnergyScaleItem(1.5,2,-999,0.94,190645,190781,0.9919,0.0020);
    InsertEnergyScaleItem(1.5,2,-999,0.94,190782,191042,0.9917,0.0022);
    InsertEnergyScaleItem(1.5,2,-999,0.94,191043,191720,0.9944,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,191721,193833,0.9917,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,193834,194116,0.9906,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,194117,194427,0.9918,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,194428,194618,0.9905,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,194619,194789,0.9940,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,194790,195111,0.9963,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,195112,195377,0.9948,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,195378,195398,0.9930,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,195399,195657,0.9929,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,195658,195918,0.9962,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,195919,196198,0.9942,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,196199,196356,0.9944,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,196357,198115,0.9923,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,198116,198940,0.9940,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,198941,199317,0.9920,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,199318,199428,0.9921,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,199429,199697,0.9916,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,199698,199832,0.9909,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,199833,199960,0.9938,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,199961,200151,0.9924,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,200152,200490,0.9942,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,200491,200991,0.9924,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,200992,201201,0.9914,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,201202,201624,0.9920,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,201625,201707,0.9919,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,201708,202059,0.9929,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,202060,202204,0.9914,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,202205,202332,0.9950,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,202333,202972,0.9920,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,202973,203002,0.9924,0.0012);
    InsertEnergyScaleItem(1.5,2,-999,0.94,203003,203852,0.9769,0.0035);
    InsertEnergyScaleItem(1.5,2,-999,0.94,203853,204099,0.9914,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,204100,204562,0.9931,0.0011);
    InsertEnergyScaleItem(1.5,2,-999,0.94,204563,205085,0.9934,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,205086,205310,0.9932,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,205311,205617,0.9916,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,205618,205825,0.9930,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,205826,206207,0.9945,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,206208,206389,0.9906,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,206390,206483,0.9957,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,206484,206597,0.9949,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,206598,206896,0.9910,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,206897,207220,0.9935,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,207221,207315,0.9934,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,207316,207489,0.9944,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,207490,207919,0.9930,0.0009);
    InsertEnergyScaleItem(1.5,2,-999,0.94,207920,208351,0.9930,0.0010);
    InsertEnergyScaleItem(1.5,2,-999,0.94,208352,208686,0.9934,0.0008);
    InsertEnergyScaleItem(1.5,2,0.94,+999,190645,190781,0.9856,0.0020);
    InsertEnergyScaleItem(1.5,2,0.94,+999,190782,191042,0.9854,0.0022);
    InsertEnergyScaleItem(1.5,2,0.94,+999,191043,191720,0.9880,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,191721,193833,0.9853,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,193834,194116,0.9842,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,194117,194427,0.9855,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,194428,194618,0.9841,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,194619,194789,0.9876,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,194790,195111,0.9900,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,195112,195377,0.9885,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,195378,195398,0.9867,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,195399,195657,0.9866,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,195658,195918,0.9898,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,195919,196198,0.9878,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,196199,196356,0.9880,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,196357,198115,0.9859,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,198116,198940,0.9877,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,198941,199317,0.9856,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,199318,199428,0.9857,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,199429,199697,0.9853,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,199698,199832,0.9845,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,199833,199960,0.9875,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,199961,200151,0.9860,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,200152,200490,0.9878,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,200491,200991,0.9861,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,200992,201201,0.9850,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,201202,201624,0.9857,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,201625,201707,0.9856,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,201708,202059,0.9866,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,202060,202204,0.9850,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,202205,202332,0.9887,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,202333,202972,0.9857,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,202973,203002,0.9860,0.0013);
    InsertEnergyScaleItem(1.5,2,0.94,+999,203003,203852,0.9705,0.0035);
    InsertEnergyScaleItem(1.5,2,0.94,+999,203853,204099,0.9851,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,204100,204562,0.9868,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,204563,205085,0.9871,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,205086,205310,0.9868,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,205311,205617,0.9852,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,205618,205825,0.9866,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,205826,206207,0.9882,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,206208,206389,0.9842,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,206390,206483,0.9894,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,206484,206597,0.9886,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,206598,206896,0.9846,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,206897,207220,0.9872,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,207221,207315,0.9870,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,207316,207489,0.9881,0.0010);
    InsertEnergyScaleItem(1.5,2,0.94,+999,207490,207919,0.9866,0.0009);
    InsertEnergyScaleItem(1.5,2,0.94,+999,207920,208351,0.9867,0.0011);
    InsertEnergyScaleItem(1.5,2,0.94,+999,208352,208686,0.9871,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,190645,190781,0.9862,0.0019);
    InsertEnergyScaleItem(2,3,-999,0.94,190782,191042,0.9810,0.0019);
    InsertEnergyScaleItem(2,3,-999,0.94,191043,191720,0.9845,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,191721,193833,0.9841,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,193834,194116,0.9843,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,194117,194427,0.9855,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,194428,194618,0.9836,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,194619,194789,0.9838,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,194790,195111,0.9853,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,195112,195377,0.9866,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,195378,195398,0.9851,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,195399,195657,0.9856,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,195658,195918,0.9846,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,195919,196198,0.9846,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,196199,196356,0.9831,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,196357,198115,0.9841,0.0008);
    InsertEnergyScaleItem(2,3,-999,0.94,198116,198940,0.9868,0.0008);
    InsertEnergyScaleItem(2,3,-999,0.94,198941,199317,0.9855,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,199318,199428,0.9858,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,199429,199697,0.9853,0.0008);
    InsertEnergyScaleItem(2,3,-999,0.94,199698,199832,0.9853,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,199833,199960,0.9873,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,199961,200151,0.9867,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,200152,200490,0.9858,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,200491,200991,0.9864,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,200992,201201,0.9883,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,201202,201624,0.9876,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,201625,201707,0.9870,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,201708,202059,0.9883,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,202060,202204,0.9877,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,202205,202332,0.9874,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,202333,202972,0.9861,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,202973,203002,0.9865,0.0012);
    InsertEnergyScaleItem(2,3,-999,0.94,203003,203852,0.9778,0.0030);
    InsertEnergyScaleItem(2,3,-999,0.94,203853,204099,0.9827,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,204100,204562,0.9864,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,204563,205085,0.9839,0.0008);
    InsertEnergyScaleItem(2,3,-999,0.94,205086,205310,0.9832,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,205311,205617,0.9816,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,205618,205825,0.9842,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,205826,206207,0.9866,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,206208,206389,0.9844,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,206390,206483,0.9873,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,206484,206597,0.9841,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,206598,206896,0.9839,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,206897,207220,0.9863,0.0010);
    InsertEnergyScaleItem(2,3,-999,0.94,207221,207315,0.9848,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,207316,207489,0.9841,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,207490,207919,0.9819,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,207920,208351,0.9863,0.0009);
    InsertEnergyScaleItem(2,3,-999,0.94,208352,208686,0.9859,0.0008);
    InsertEnergyScaleItem(2,3,0.94,+999,190645,190781,0.9812,0.0019);
    InsertEnergyScaleItem(2,3,0.94,+999,190782,191042,0.9760,0.0019);
    InsertEnergyScaleItem(2,3,0.94,+999,191043,191720,0.9795,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,191721,193833,0.9791,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,193834,194116,0.9793,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,194117,194427,0.9805,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,194428,194618,0.9786,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,194619,194789,0.9788,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,194790,195111,0.9803,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,195112,195377,0.9817,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,195378,195398,0.9801,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,195399,195657,0.9806,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,195658,195918,0.9797,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,195919,196198,0.9796,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,196199,196356,0.9781,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,196357,198115,0.9791,0.0008);
    InsertEnergyScaleItem(2,3,0.94,+999,198116,198940,0.9818,0.0008);
    InsertEnergyScaleItem(2,3,0.94,+999,198941,199317,0.9805,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,199318,199428,0.9808,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,199429,199697,0.9803,0.0008);
    InsertEnergyScaleItem(2,3,0.94,+999,199698,199832,0.9803,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,199833,199960,0.9823,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,199961,200151,0.9817,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,200152,200490,0.9809,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,200491,200991,0.9814,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,200992,201201,0.9833,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,201202,201624,0.9827,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,201625,201707,0.9820,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,201708,202059,0.9833,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,202060,202204,0.9827,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,202205,202332,0.9825,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,202333,202972,0.9812,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,202973,203002,0.9816,0.0012);
    InsertEnergyScaleItem(2,3,0.94,+999,203003,203852,0.9728,0.0030);
    InsertEnergyScaleItem(2,3,0.94,+999,203853,204099,0.9777,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,204100,204562,0.9814,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,204563,205085,0.9789,0.0008);
    InsertEnergyScaleItem(2,3,0.94,+999,205086,205310,0.9782,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,205311,205617,0.9766,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,205618,205825,0.9792,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,205826,206207,0.9816,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,206208,206389,0.9794,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,206390,206483,0.9824,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,206484,206597,0.9791,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,206598,206896,0.9789,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,206897,207220,0.9813,0.0010);
    InsertEnergyScaleItem(2,3,0.94,+999,207221,207315,0.9798,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,207316,207489,0.9791,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,207490,207919,0.9769,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,207920,208351,0.9813,0.0009);
    InsertEnergyScaleItem(2,3,0.94,+999,208352,208686,0.9809,0.0008);

    // https://twiki.cern.ch/twiki/pub/CMS/ECALELF/22Jan2012-runDepPowheg-noR9shift-outProfile-scaleStep2smearing-Et_25-trigger-noPF-FitResult.config
    InsertEnergySmearingItem(0,1,-999,0.94,0,999999  ,0.86,0.0181);
    InsertEnergySmearingItem(0,1,0.94,999,0,999999   ,0.75,0.03  );
    InsertEnergySmearingItem(1,1.5,-999,0.94,0,999999,1.88,0.0188);
    InsertEnergySmearingItem(1,1.5,0.94,999,0,999999 ,1.22,0.0992);
    InsertEnergySmearingItem(1.5,2,-999,0.94,0,999999,1.98,0.0323);
    InsertEnergySmearingItem(1.5,2,0.94,999,0,999999 ,1.63,0.0878);
    InsertEnergySmearingItem(2,3,-999,0.94,0,999999  ,1.92,0.0444);
    InsertEnergySmearingItem(2,3,0.94,999,0,999999   ,1.86,0.0353);

  }

};

void DiPhotonMiniTree::InsertEnergyScaleItem(float meta, float Meta, float mr9, float Mr9, int mrun, int Mrun, float val, float err){
  
  struct_escale_item a;
  a.meta = meta;
  a.Meta = Meta;
  a.mr9 = mr9;
  a.Mr9 = Mr9;
  a.mrun = mrun;
  a.Mrun = Mrun;
  a.val = val;
  a.err = err;

  energyScaleDatabase.push_back(a);

};

void DiPhotonMiniTree::InsertEnergySmearingItem(float meta, float Meta, float mr9, float Mr9, int mrun, int Mrun, float val, float err){
  
  struct_escale_item a;
  a.meta = meta;
  a.Meta = Meta;
  a.mr9 = mr9;
  a.Mr9 = Mr9;
  a.mrun = mrun;
  a.Mrun = Mrun;
  a.val = val;
  a.err = err;

  energySmearingDatabase.push_back(a);

};

void DiPhotonMiniTree::JECJERCorrection(StatusScaleUpScaleDown status_escale, StatusScaleUpScaleDown status_esmear, vector<pair<int,float> > &ordering_jets){

  // ADD PROTECTIONS FOR ZERO, HIGH RHO ETC.

  // JEC

  if (!MyJetCorrector) {
    cout << "ERROR: NO JET ENERGY CORRECTOR INITIALIZED!!!" << endl;
    return;
  }

  ordering_jets.clear();
  
  for (int i=0; i<fTR->NJets; i++){

    if (!(fTR->JEcorr[i]>0)) continue;

    float invcorr = 1./fTR->JEcorr[i];
    fTR->JPt[i] *= invcorr;
    fTR->JE[i] *= invcorr;
    fTR->JEcorr[i] *= invcorr;
    
    if (status_escale!=kUncorrected){

      std::vector<float> corr = MyJetCorrector->getCorrPtECorr(fTR->JPt[i],fTR->JEta[i],fTR->JE[i],fTR->JEcorr[i],fTR->JArea[i],fTR->Rho);
      fTR->JPt[i] = corr.at(0);
      fTR->JE[i] = corr.at(1);
      fTR->JEcorr[i] = corr.at(2);
      // WARNING: the other jet quantities are not recorrected, do not use them!!!
      
      {
	float eta = fTR->JEta[i];
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	MyJetCorrector->fJetUncertainty->setJetPt(fTR->JPt[i]);
	MyJetCorrector->fJetUncertainty->setJetEta(eta);
	float jecuncup = MyJetCorrector->fJetUncertainty->getUncertainty(true);
	MyJetCorrector->fJetUncertainty->setJetPt(fTR->JPt[i]);
	MyJetCorrector->fJetUncertainty->setJetEta(eta);
	float jecuncdown = MyJetCorrector->fJetUncertainty->getUncertainty(false);
	float addcorr = 1;
	//	if (debug) cout << jecuncup << " " << jecuncdown << endl;
	if (status_escale==kShiftUp) addcorr=1+jecuncup;
	if (status_escale==kShiftDown) addcorr=1-jecuncdown;
	fTR->JPt[i] *= addcorr;
	fTR->JE[i] *= addcorr;
	fTR->JEcorr[i] *= addcorr;
      }

      
    }

    // JER, good for both 2011 and 2012

    {

      float thiseta = fabs(fTR->JEta[i]);
      float smear,smearup,smeardown;

      if (thiseta<0.5) smear=1.052;
      else if (thiseta<1.1) smear=1.057;
      else if (thiseta<1.7) smear=1.096;
      else if (thiseta<2.3) smear=1.134;
      else smear=1.288;

      if (thiseta<0.5) smearup=1.115;
      else if (thiseta<1.1) smearup=1.114;
      else if (thiseta<1.7) smearup=1.161;
      else if (thiseta<2.3) smearup=1.228;
      else smearup=1.488;

      if (thiseta<0.5) smeardown=0.990;
      else if (thiseta<1.1) smeardown=1.001;
      else if (thiseta<1.7) smeardown=1.032;
      else if (thiseta<2.3) smeardown=1.142;
      else smeardown=1.089;

      float pt = fTR->JPt[i];
      float ptgen = pt;
      if (!isdata && fTR->NGenJets>0 && fTR->JGenJetIndex[i]>=0) ptgen = fTR->GenJetPt[fTR->JGenJetIndex[i]]; // no smearing for jets not matched to genjets
      else if (!isdata && debug) {
	cout << "no match to genjet " << isdata << " " << fTR->NGenJets << " " << fTR->JGenJetIndex[i] << endl;
	for (int k=0; k<fTR->NGenJets; k++) cout << "genjet" << k << " " << fTR->GenJetPt[k] << " " << fTR->GenJetEta[k] << " " << fTR->GenJetPhi[k] << endl;
      }
      float jer = TMath::Max(float(0),float((ptgen+smear*(pt-ptgen))))/pt;
      float jerup = TMath::Max(float(0),float((ptgen+smearup*(pt-ptgen))))/pt;
      float jerdown = TMath::Max(float(0),float((ptgen+smeardown*(pt-ptgen))))/pt;
      float myjer = 1;

      if (status_esmear==kUncorrected) ;
      else if (status_esmear==kCorrectedNoShift) myjer=jer;
      else if (status_esmear==kShiftUp) myjer=jerup;
      else if (status_esmear==kShiftDown) myjer=jerdown;
      else assert(false);
    
      fTR->JPt[i] *= myjer;
      fTR->JE[i] *= myjer;
      fTR->JEcorr[i] *= myjer;

    }

    ordering_jets.push_back(std::pair<int,float>(i,fTR->JPt[i]));

  }

  sort(ordering_jets.begin(),ordering_jets.end(),indexComparator);


//  float OnTheFlyCorrections::getJECUncertainty(float pt, float eta){
//    if      (eta> 5.0) eta = 5.0;
//    else if (eta<-5.0) eta =-5.0;
//    fJetUncertainty->setJetPt(pt);
//    fJetUncertainty->setJetEta(eta);
//    float uncert= fJetUncertainty->getUncertainty(true);
//    return uncert;
//
//  }
//
//  std::vector< float > OnTheFlyCorrections::getCorrPtECorr(float jetpt, float jeteta, float jetenergy, float jecorr, float jetarea, float rho){
//    // give this function a corrected jet and it will return a vector with the 
//    // corrected jet-pt, jet-energy and the new correction according to the global tag that is set
//
//    float rawpt = jetpt/jecorr;
//    float rawe  = jetenergy/jecorr;
//    fJetCorrector->setJetPt(rawpt); // raw-pT here
//    fJetCorrector->setJetEta(jeteta);
//    fJetCorrector->setRho(rho);
//    fJetCorrector->setJetA(jetarea);
//
//    float corr = fJetCorrector->getCorrection();
//    std::vector< float > corrected;
//
//    // new jetpt = old jet-pt * new correction / old correction
//    corrected.push_back(corr*rawpt); // new corrected pt as 0th item
//    corrected.push_back(corr*rawe);  // new corrected energy as 1st item
//    corrected.push_back(corr);       // new correction as 2nd item
//    return corrected;
//  }

};

MatchingStatus DiPhotonMiniTree::check_matching_status(int phoindex){

  assert (!isdata);

  int status = fTR->PhoMCmatchexitcode[phoindex];
  float geniso = (fTR->PhoMCmatchindex[phoindex]>=0) ? fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[phoindex]] : 999;

  return determine_matchingstatus(status,geniso);

};

MatchingStatus DiPhotonMiniTree::determine_matchingstatus(int status, float geniso){

  assert (!isdata);

  if ((status==1 || status==2 || dataset_id==sherpa_dataset_id) && (geniso<5)) return kSignal;
  else if (status==4) return kElectron;
  else return kBackground;

};

void DiPhotonMiniTree::FixMatchingStatusElectrons(int phoindex){

  // Correct matching code to 4 if matched to electron

  if (fTR->PhoMCmatchexitcode[phoindex]!=0) return;
  
  for (int i=0; i<fTR->NGenLeptons; i++){

    if (abs(fTR->GenLeptonID[i])!=11) continue;

    float dr = Util::GetDeltaR(fTR->GenLeptonEta[i],fTR->PhoEta[phoindex],fTR->GenLeptonPhi[i],fTR->PhoPhi[phoindex]);
    if(dr > 0.2) continue;

    double ndpt = fabs(fTR->GenLeptonPt[i] - fTR->PhoPt[phoindex])/fTR->GenLeptonPt[i];
    if(ndpt > 2.) continue;

    fTR->PhoMCmatchexitcode[phoindex]=4;
    return;

  }

};

void DiPhotonMiniTree::FixMatchingStatusSherpa(int phoindex){

  assert (dataset_id==sherpa_dataset_id);

  // Correct matching code to 2 in any case, if matched to status 1 photon
  // irrespectively of the mother particle

  if (fTR->PhoMCmatchexitcode[phoindex]<=0) return;
  else if (fTR->PhoMCmatchexitcode[phoindex]==4) return;
  else fTR->PhoMCmatchexitcode[phoindex]=2;
  
};
