#include "helper/Utilities.hh"
#include "DiPhotonMiniTree.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

DiPhotonMiniTree::DiPhotonMiniTree(TreeReader *tr, std::string dataType, Float_t aw, Float_t* _kfac, Float_t _minthrpfphotoncandEB, Float_t _minthrpfphotoncandEE, bool _isstep2, TString _input_filename, UInt_t _uuid) : UserAnalysisBase(tr), fDataType_(dataType), AddWeight(aw), kfactors(_kfac), global_minthrpfphotoncandEB(_minthrpfphotoncandEB), global_minthrpfphotoncandEE(_minthrpfphotoncandEE), isstep2(_isstep2), input_filename(_input_filename), uuid(_uuid){
  Util::SetStyle();	
  if (fDataType_ == "mc") isdata=false;
  else if (fDataType_ == "data") isdata=true; 
  else {
    std::cout << "wrong data type" << std::endl;
    assert(1==0);
  }
  randomgen = new TRandom3(0);

  global_linkbyrechit_enlargement = 0.25; // xtal_size_eff = (1+global_linkbyrechit_enlargement)*xtal_size

  eegeom = TGeoPara(1,1,1,0,0,0);

//  // 020616 no cleaning numbers (superseded: should be updated with new ones after removing pf ch iso from presel) XXX
//  // EA SET TO ZERO HARDCODED BELOW, these values are ignored;
//  float _binsdef_single_gamma_EB_eta[n_templates_EB+1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442}; 
//  float _binsdef_single_gamma_EE_eta[n_templates_EE+1] = {1.566,1.653,1.8,2,2.2,2.5};
//  float _eff_areas_EB_photon_data[n_templates_EB] = {2.601118e-01,2.584915e-01,2.640072e-01,2.656851e-01,2.564615e-01,2.396511e-01,1.645776e-01};
//  float _eff_areas_EE_photon_data[n_templates_EE] = {5.783452e-02,8.321881e-02,1.177009e-01,1.422445e-01,1.139434e-01};
//  float _eff_areas_EB_charged_data[n_templates_EB] = {2.526058e-02,2.619860e-02,2.118318e-02,2.078713e-02,2.066629e-02,1.908721e-02,1.732135e-02};
//  float _eff_areas_EE_charged_data[n_templates_EE] = {1.538706e-02,1.737151e-02,1.410831e-02,1.558452e-02,1.056310e-02};
//  float _eff_areas_EB_neutral_data[n_templates_EB] = {9.534780e-02,8.826486e-02,8.199201e-02,8.729101e-02,1.018310e-01,1.244805e-01,1.542784e-01};
//  float _eff_areas_EE_neutral_data[n_templates_EE] = {1.547500e-01,1.480421e-01,1.391524e-01,1.546148e-01,2.487922e-01};
//  float _eff_areas_EB_photon_MC[n_templates_EB] = {2.728305e-01,2.747126e-01,2.742736e-01,2.689392e-01,2.642000e-01,2.430941e-01,1.701134e-01};
//  float _eff_areas_EE_photon_MC[n_templates_EE] = {5.066866e-02,8.511977e-02,1.342678e-01,1.693861e-01,1.448627e-01};
//  float _eff_areas_EB_charged_MC[n_templates_EB] = {2.333960e-02,2.682862e-02,2.747902e-02,2.494792e-02,2.206384e-02,1.864648e-02,2.470047e-02};
//  float _eff_areas_EE_charged_MC[n_templates_EE] = {2.264766e-02,1.934945e-02,1.834243e-02,2.223225e-02,8.191396e-03};
//  float _eff_areas_EB_neutral_MC[n_templates_EB] = {7.527103e-02,7.490141e-02,7.510347e-02,7.806471e-02,9.715760e-02,1.193150e-01,1.449191e-01};
//  float _eff_areas_EE_neutral_MC[n_templates_EE] = {1.507419e-01,1.383201e-01,1.489503e-01,1.460166e-01,2.548731e-01};
//
//
//
//  binsdef_single_gamma_EB_eta.assign(_binsdef_single_gamma_EB_eta,_binsdef_single_gamma_EB_eta+n_templates_EB+1);
//  binsdef_single_gamma_EE_eta.assign(_binsdef_single_gamma_EE_eta,_binsdef_single_gamma_EE_eta+n_templates_EE+1);
//  eff_areas_EB_photon_data.assign(_eff_areas_EB_photon_data ,_eff_areas_EB_photon_data +n_templates_EB);
//  eff_areas_EE_photon_data.assign(_eff_areas_EE_photon_data ,_eff_areas_EE_photon_data +n_templates_EE);
//  eff_areas_EB_charged_data.assign(_eff_areas_EB_charged_data ,_eff_areas_EB_charged_data +n_templates_EB);
//  eff_areas_EE_charged_data.assign(_eff_areas_EE_charged_data ,_eff_areas_EE_charged_data +n_templates_EE);
//  eff_areas_EB_neutral_data.assign(_eff_areas_EB_neutral_data ,_eff_areas_EB_neutral_data +n_templates_EB);
//  eff_areas_EE_neutral_data.assign(_eff_areas_EE_neutral_data ,_eff_areas_EE_neutral_data +n_templates_EE);
//  eff_areas_EB_photon_MC.assign(_eff_areas_EB_photon_MC ,_eff_areas_EB_photon_MC +n_templates_EB);
//  eff_areas_EE_photon_MC.assign(_eff_areas_EE_photon_MC ,_eff_areas_EE_photon_MC +n_templates_EE);
//  eff_areas_EB_charged_MC.assign(_eff_areas_EB_charged_MC ,_eff_areas_EB_charged_MC +n_templates_EB);
//  eff_areas_EE_charged_MC.assign(_eff_areas_EE_charged_MC ,_eff_areas_EE_charged_MC +n_templates_EE);
//  eff_areas_EB_neutral_MC.assign(_eff_areas_EB_neutral_MC ,_eff_areas_EB_neutral_MC +n_templates_EB);
//  eff_areas_EE_neutral_MC.assign(_eff_areas_EE_neutral_MC ,_eff_areas_EE_neutral_MC +n_templates_EE);
//

}

DiPhotonMiniTree::~DiPhotonMiniTree(){
  delete randomgen;
}

void DiPhotonMiniTree::Begin(){

  cout << "Begin" << endl;

  fOutputFile->cd();

  treename[0] = TString("Tree_2Dstandard_selection");
  treename[1] = TString("Tree_1Drandomcone_template");
  treename[2] = TString("Tree_1Dsideband_template");
  treename[3] = TString("Tree_2DZee_pixelvetoreversed_selection");
  treename[4] = TString("Tree_1Dpreselection");
  treename[5] = TString("Tree_1Dselection");
  treename[6] = TString("Tree_2Drandomcone_template");
  treename[7] = TString("Tree_2Drandomconesideband_template");
  treename[8] = TString("Tree_2Dsideband_template");
  treename[9] = TString("Tree_2Dstandard_preselection");
  treename[10] = TString("Tree_2DZmumu_selection");
  treename[11] = TString("Tree_1Dsignal_template");
  treename[12] = TString("Tree_1Dbackground_template");
  treename[13] = TString("Tree_2Dtruesigsig_template");
  treename[14] = TString("Tree_2Dtruesigbkg_template");
  treename[15] = TString("Tree_2Dtruebkgbkg_template");
  treename[16] = TString("Tree_2Drconeplusgenfake_template");
  treename[17] = TString("Tree_2Dgenpromptplussideband_template");

  for (int i=0; i<18; i++){
    OutputTree[i] = new TTree(treename[i].Data(),treename[i].Data());
    is2d[i] = (treename[i].First("2")>-1);
  }

  for (int i=0; i<18; i++){

  OutputTree[i]->Branch("event_fileuuid",&event_fileuuid,"event_fileuuid/i");

  OutputTree[i]->Branch("event_run",&event_run,"event_run/I");
  OutputTree[i]->Branch("event_lumi",&event_lumi,"event_lumi/I");
  OutputTree[i]->Branch("event_number",&event_number,"event_number/I");

  OutputTree[i]->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
  OutputTree[i]->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");

  OutputTree[i]->Branch("event_weight",&event_weight,"event_weight/F");
  OutputTree[i]->Branch("event_rho",&event_rho,"event_rho/F");
  OutputTree[i]->Branch("event_sigma",&event_sigma,"event_sigma/F");
  OutputTree[i]->Branch("event_nPU",&event_nPU,"event_nPU/I");
  OutputTree[i]->Branch("event_nRecVtx",&event_nRecVtx,"event_nRecVtx/I");
  OutputTree[i]->Branch("event_pass12whoissiglike",&event_pass12whoissiglike,"event_pass12whoissiglike/I");

  OutputTree[i]->Branch("dipho_mgg_photon",&dipho_mgg_photon,"dipho_mgg_photon/F");

  OutputTree[i]->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
  OutputTree[i]->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
  OutputTree[i]->Branch("pholead_phi",&pholead_phi,"pholead_phi/F");
  OutputTree[i]->Branch("photrail_phi",&photrail_phi,"photrail_phi/F");
  OutputTree[i]->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
  OutputTree[i]->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
  OutputTree[i]->Branch("pholead_energy",&pholead_energy,"pholead_energy/F");
  OutputTree[i]->Branch("photrail_energy",&photrail_energy,"photrail_energy/F");

  OutputTree[i]->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
  OutputTree[i]->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
  OutputTree[i]->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
  OutputTree[i]->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");
  
  OutputTree[i]->Branch("pholead_PhoHasPixSeed",&pholead_PhoHasPixSeed,"pholead_PhoHasPixSeed/I");
  OutputTree[i]->Branch("photrail_PhoHasPixSeed",&photrail_PhoHasPixSeed,"photrail_PhoHasPixSeed/I");

  OutputTree[i]->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
  OutputTree[i]->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");
  OutputTree[i]->Branch("pholead_sieie",&pholead_sieie,"pholead_sieie/F");
  OutputTree[i]->Branch("photrail_sieie",&photrail_sieie,"photrail_sieie/F");
  OutputTree[i]->Branch("pholead_hoe",&pholead_hoe,"pholead_hoe/F");
  OutputTree[i]->Branch("photrail_hoe",&photrail_hoe,"photrail_hoe/F");

  OutputTree[i]->Branch("pholead_PhoSCRemovalPFIsoCharged",&pholead_PhoSCRemovalPFIsoCharged,"pholead_PhoSCRemovalPFIsoCharged/F");
  OutputTree[i]->Branch("photrail_PhoSCRemovalPFIsoCharged",&photrail_PhoSCRemovalPFIsoCharged,"photrail_PhoSCRemovalPFIsoCharged/F");
  OutputTree[i]->Branch("pholead_PhoSCRemovalPFIsoNeutral",&pholead_PhoSCRemovalPFIsoNeutral,"pholead_PhoSCRemovalPFIsoNeutral/F");
  OutputTree[i]->Branch("photrail_PhoSCRemovalPFIsoNeutral",&photrail_PhoSCRemovalPFIsoNeutral,"photrail_PhoSCRemovalPFIsoNeutral/F");
  OutputTree[i]->Branch("pholead_PhoSCRemovalPFIsoPhoton",&pholead_PhoSCRemovalPFIsoPhoton,"pholead_PhoSCRemovalPFIsoPhoton/F");
  OutputTree[i]->Branch("photrail_PhoSCRemovalPFIsoPhoton",&photrail_PhoSCRemovalPFIsoPhoton,"photrail_PhoSCRemovalPFIsoPhoton/F");
  OutputTree[i]->Branch("pholead_PhoSCRemovalPFIsoCombined",&pholead_PhoSCRemovalPFIsoCombined,"pholead_PhoSCRemovalPFIsoCombined/F");
  OutputTree[i]->Branch("photrail_PhoSCRemovalPFIsoCombined",&photrail_PhoSCRemovalPFIsoCombined,"photrail_PhoSCRemovalPFIsoCombined/F");

  OutputTree[i]->Branch("pholead_PhoPassConversionVeto",&pholead_PhoPassConversionVeto,"pholead_PhoPassConversionVeto/I");
  OutputTree[i]->Branch("photrail_PhoPassConversionVeto",&photrail_PhoPassConversionVeto,"photrail_PhoPassConversionVeto/I");

  OutputTree[i]->Branch("pholead_GenPhotonIsoDR04",&pholead_GenPhotonIsoDR04,"pholead_GenPhotonIsoDR04/F");
  OutputTree[i]->Branch("photrail_GenPhotonIsoDR04",&photrail_GenPhotonIsoDR04,"photrail_GenPhotonIsoDR04/F");

  OutputTree[i]->Branch("pholead_PhoMCmatchexitcode",&pholead_PhoMCmatchexitcode,"pholead_PhoMCmatchexitcode/I");
  OutputTree[i]->Branch("photrail_PhoMCmatchexitcode",&photrail_PhoMCmatchexitcode,"photrail_PhoMCmatchexitcode/I");

//  OutputTree[i]->Branch("pholead_scarea",&pholead_scarea,"pholead_scarea/F");
//  OutputTree[i]->Branch("pholead_scareaSF",&pholead_scareaSF,"pholead_scareaSF/F");
//  OutputTree[i]->Branch("photrail_scarea",&photrail_scarea,"photrail_scarea/F");
//  OutputTree[i]->Branch("photrail_scareaSF",&photrail_scareaSF,"photrail_scareaSF/F");

  if (do_recalc_isolation){
  OutputTree[i]->Branch("pholead_Npfcandphotonincone",&pholead_Npfcandphotonincone,"pholead_Npfcandphotonincone/I");
  OutputTree[i]->Branch("pholead_Npfcandchargedincone",&pholead_Npfcandchargedincone,"pholead_Npfcandchargedincone/I");
  OutputTree[i]->Branch("pholead_Npfcandneutralincone",&pholead_Npfcandneutralincone,"pholead_Npfcandneutralincone/I");
  OutputTree[i]->Branch("photrail_Npfcandphotonincone",&photrail_Npfcandphotonincone,"photrail_Npfcandphotonincone/I");
  OutputTree[i]->Branch("photrail_Npfcandchargedincone",&photrail_Npfcandchargedincone,"photrail_Npfcandchargedincone/I");
  OutputTree[i]->Branch("photrail_Npfcandneutralincone",&photrail_Npfcandneutralincone,"photrail_Npfcandneutralincone/I");
  OutputTree[i]->Branch("pholead_photonpfcandenergies",&pholead_photonpfcandenergies,"pholead_photonpfcandenergies[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcandenergies",&pholead_chargedpfcandenergies,"pholead_chargedpfcandenergies[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcandenergies",&pholead_neutralpfcandenergies,"pholead_neutralpfcandenergies[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcandets",&pholead_photonpfcandets,"pholead_photonpfcandets[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcandets",&pholead_chargedpfcandets,"pholead_chargedpfcandets[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcandets",&pholead_neutralpfcandets,"pholead_neutralpfcandets[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcanddetas",&pholead_photonpfcanddetas,"pholead_photonpfcanddetas[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcanddetas",&pholead_chargedpfcanddetas,"pholead_chargedpfcanddetas[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcanddetas",&pholead_neutralpfcanddetas,"pholead_neutralpfcanddetas[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("pholead_photonpfcanddphis",&pholead_photonpfcanddphis,"pholead_photonpfcanddphis[pholead_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("pholead_chargedpfcanddphis",&pholead_chargedpfcanddphis,"pholead_chargedpfcanddphis[pholead_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("pholead_neutralpfcanddphis",&pholead_neutralpfcanddphis,"pholead_neutralpfcanddphis[pholead_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcandenergies",&photrail_photonpfcandenergies,"photrail_photonpfcandenergies[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcandenergies",&photrail_chargedpfcandenergies,"photrail_chargedpfcandenergies[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcandenergies",&photrail_neutralpfcandenergies,"photrail_neutralpfcandenergies[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcandets",&photrail_photonpfcandets,"photrail_photonpfcandets[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcandets",&photrail_chargedpfcandets,"photrail_chargedpfcandets[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcandets",&photrail_neutralpfcandets,"photrail_neutralpfcandets[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcanddetas",&photrail_photonpfcanddetas,"photrail_photonpfcanddetas[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcanddetas",&photrail_chargedpfcanddetas,"photrail_chargedpfcanddetas[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcanddetas",&photrail_neutralpfcanddetas,"photrail_neutralpfcanddetas[photrail_Npfcandneutralincone]/F");
  OutputTree[i]->Branch("photrail_photonpfcanddphis",&photrail_photonpfcanddphis,"photrail_photonpfcanddphis[photrail_Npfcandphotonincone]/F");
  OutputTree[i]->Branch("photrail_chargedpfcanddphis",&photrail_chargedpfcanddphis,"photrail_chargedpfcanddphis[photrail_Npfcandchargedincone]/F");
  OutputTree[i]->Branch("photrail_neutralpfcanddphis",&photrail_neutralpfcanddphis,"photrail_neutralpfcanddphis[photrail_Npfcandneutralincone]/F");
  }

//  // IF SCAN ROTATED CONES
//  OutputTree[i]->Branch("pholead_test_rotatedphotoniso",&pholead_test_rotatedphotoniso,"pholead_test_rotatedphotoniso[50]/F");
//  OutputTree[i]->Branch("pholead_test_rotatedwithcheckphotoniso",&pholead_test_rotatedwithcheckphotoniso,"pholead_test_rotatedwithcheckphotoniso[50]/F");

  OutputTree[i]->Branch("allphotonpfcand_count",&allphotonpfcand_count,"allphotonpfcand_count/I");
  OutputTree[i]->Branch("allphotonpfcand_pt",&allphotonpfcand_pt,"allphotonpfcand_pt[allphotonpfcand_count]/F");
  OutputTree[i]->Branch("allphotonpfcand_eta",&allphotonpfcand_eta,"allphotonpfcand_eta[allphotonpfcand_count]/F");
  OutputTree[i]->Branch("allphotonpfcand_phi",&allphotonpfcand_phi,"allphotonpfcand_phi[allphotonpfcand_count]/F");
  OutputTree[i]->Branch("allphotonpfcand_vx",&allphotonpfcand_vx,"allphotonpfcand_vx[allphotonpfcand_count]/F");
  OutputTree[i]->Branch("allphotonpfcand_vy",&allphotonpfcand_vy,"allphotonpfcand_vy[allphotonpfcand_count]/F");
  OutputTree[i]->Branch("allphotonpfcand_vz",&allphotonpfcand_vz,"allphotonpfcand_vz[allphotonpfcand_count]/F");

  OutputTree[i]->Branch("phoiso_template_1event_sigsig_1",&phoiso_template_1event_sigsig_1,Form("phoiso_template_1event_sigsig_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_sigsig_2",&phoiso_template_1event_sigsig_2,Form("phoiso_template_1event_sigsig_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_sigbkg_1",&phoiso_template_1event_sigbkg_1,Form("phoiso_template_1event_sigbkg_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_sigbkg_2",&phoiso_template_1event_sigbkg_2,Form("phoiso_template_1event_sigbkg_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_bkgsig_1",&phoiso_template_1event_bkgsig_1,Form("phoiso_template_1event_bkgsig_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_bkgsig_2",&phoiso_template_1event_bkgsig_2,Form("phoiso_template_1event_bkgsig_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_bkgbkg_1",&phoiso_template_1event_bkgbkg_1,Form("phoiso_template_1event_bkgbkg_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_1event_bkgbkg_2",&phoiso_template_1event_bkgbkg_2,Form("phoiso_template_1event_bkgbkg_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_sigsig_1",&phoiso_template_2events_sigsig_1,Form("phoiso_template_2events_sigsig_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_sigsig_2",&phoiso_template_2events_sigsig_2,Form("phoiso_template_2events_sigsig_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_sigbkg_1",&phoiso_template_2events_sigbkg_1,Form("phoiso_template_2events_sigbkg_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_sigbkg_2",&phoiso_template_2events_sigbkg_2,Form("phoiso_template_2events_sigbkg_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_bkgsig_1",&phoiso_template_2events_bkgsig_1,Form("phoiso_template_2events_bkgsig_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_bkgsig_2",&phoiso_template_2events_bkgsig_2,Form("phoiso_template_2events_bkgsig_2[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_bkgbkg_1",&phoiso_template_2events_bkgbkg_1,Form("phoiso_template_2events_bkgbkg_1[%d]/F",nclosest));
  OutputTree[i]->Branch("phoiso_template_2events_bkgbkg_2",&phoiso_template_2events_bkgbkg_2,Form("phoiso_template_2events_bkgbkg_2[%d]/F",nclosest));

  // rewinfo = {eta1, eta2, pt1, pt2, rho, sigma}
  OutputTree[i]->Branch("rewinfo_template_1event_sigsig_1",&rewinfo_template_1event_sigsig_1,Form("rewinfo_template_1event_sigsig_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_sigsig_2",&rewinfo_template_1event_sigsig_2,Form("rewinfo_template_1event_sigsig_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_sigbkg_1",&rewinfo_template_1event_sigbkg_1,Form("rewinfo_template_1event_sigbkg_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_sigbkg_2",&rewinfo_template_1event_sigbkg_2,Form("rewinfo_template_1event_sigbkg_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_bkgsig_1",&rewinfo_template_1event_bkgsig_1,Form("rewinfo_template_1event_bkgsig_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_bkgsig_2",&rewinfo_template_1event_bkgsig_2,Form("rewinfo_template_1event_bkgsig_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_bkgbkg_1",&rewinfo_template_1event_bkgbkg_1,Form("rewinfo_template_1event_bkgbkg_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_1event_bkgbkg_2",&rewinfo_template_1event_bkgbkg_2,Form("rewinfo_template_1event_bkgbkg_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_sigsig_1",&rewinfo_template_2events_sigsig_1,Form("rewinfo_template_2events_sigsig_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_sigsig_2",&rewinfo_template_2events_sigsig_2,Form("rewinfo_template_2events_sigsig_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_sigbkg_1",&rewinfo_template_2events_sigbkg_1,Form("rewinfo_template_2events_sigbkg_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_sigbkg_2",&rewinfo_template_2events_sigbkg_2,Form("rewinfo_template_2events_sigbkg_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_bkgsig_1",&rewinfo_template_2events_bkgsig_1,Form("rewinfo_template_2events_bkgsig_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_bkgsig_2",&rewinfo_template_2events_bkgsig_2,Form("rewinfo_template_2events_bkgsig_2[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_bkgbkg_1",&rewinfo_template_2events_bkgbkg_1,Form("rewinfo_template_2events_bkgbkg_1[%d]/F",nclosest*6));
  OutputTree[i]->Branch("rewinfo_template_2events_bkgbkg_2",&rewinfo_template_2events_bkgbkg_2,Form("rewinfo_template_2events_bkgbkg_2[%d]/F",nclosest*6));

  OutputTree[i]->Branch("vetoobjects_count",&vetoobjects_count,"vetoobjects_count/I");
  OutputTree[i]->Branch("vetoobjects_pt",&vetoobjects_pt,"vetoobjects_pt[vetoobjects_count]/F");
  OutputTree[i]->Branch("vetoobjects_eta",&vetoobjects_eta,"vetoobjects_eta[vetoobjects_count]/F");
  OutputTree[i]->Branch("vetoobjects_phi",&vetoobjects_phi,"vetoobjects_phi[vetoobjects_count]/F");
  OutputTree[i]->Branch("vetoobjects_type",&vetoobjects_type,"vetoobjects_type[vetoobjects_count]/I");

  }

  inputtree_isinitialized = false;

  LightTreeGenReco = new TTree("LightTreeGenReco","LightTreeGenReco");

  LightTreeGenReco->Branch("event_run",&event_run,"event_run/I");
  LightTreeGenReco->Branch("event_lumi",&event_lumi,"event_lumi/I");
  LightTreeGenReco->Branch("event_number",&event_number,"event_number/I");

  LightTreeGenReco->Branch("event_luminormfactor",&event_luminormfactor,"event_luminormfactor/F");
  LightTreeGenReco->Branch("event_Kfactor",&event_Kfactor,"event_Kfactor/F");
  LightTreeGenReco->Branch("event_weight",&event_weight,"event_weight/F");

  LightTreeGenReco->Branch("pholead_pt",&pholead_pt,"pholead_pt/F");
  LightTreeGenReco->Branch("photrail_pt",&photrail_pt,"photrail_pt/F");
  LightTreeGenReco->Branch("pholead_eta",&pholead_eta,"pholead_eta/F");
  LightTreeGenReco->Branch("photrail_eta",&photrail_eta,"photrail_eta/F");
  LightTreeGenReco->Branch("pholead_phi",&pholead_phi,"pholead_phi/F");
  LightTreeGenReco->Branch("photrail_phi",&photrail_phi,"photrail_phi/F");
  LightTreeGenReco->Branch("pholead_SCeta",&pholead_SCeta,"pholead_SCeta/F");
  LightTreeGenReco->Branch("photrail_SCeta",&photrail_SCeta,"photrail_SCeta/F");
  LightTreeGenReco->Branch("pholead_SCphi",&pholead_SCphi,"pholead_SCphi/F");
  LightTreeGenReco->Branch("photrail_SCphi",&photrail_SCphi,"photrail_SCphi/F");
  LightTreeGenReco->Branch("pholead_r9",&pholead_r9,"pholead_r9/F");
  LightTreeGenReco->Branch("photrail_r9",&photrail_r9,"photrail_r9/F");

  LightTreeGenReco->Branch("pholead_GEN_pt",&pholead_GEN_pt,"pholead_GEN_pt/F");
  LightTreeGenReco->Branch("photrail_GEN_pt",&photrail_GEN_pt,"photrail_GEN_pt/F");
  LightTreeGenReco->Branch("pholead_GEN_eta",&pholead_GEN_eta,"pholead_GEN_eta/F");
  LightTreeGenReco->Branch("photrail_GEN_eta",&photrail_GEN_eta,"photrail_GEN_eta/F");
  LightTreeGenReco->Branch("pholead_GEN_phi",&pholead_GEN_phi,"pholead_GEN_phi/F");
  LightTreeGenReco->Branch("photrail_GEN_phi",&photrail_GEN_phi,"photrail_GEN_phi/F");

  LightTreeGenReco->Branch("reco_has_matched_gen_no_acceptance",&tree_reco_has_matched_gen_no_acceptance,"reco_has_matched_gen_no_acceptance/O");
  LightTreeGenReco->Branch("reco_has_matched_gen_within_acceptance",&tree_reco_has_matched_gen_within_acceptance,"reco_has_matched_gen_within_acceptance/O");
  LightTreeGenReco->Branch("reco_has_matched_gen_outside_acceptance",&tree_reco_has_matched_gen_outside_acceptance,"reco_has_matched_gen_outside_acceptance/O");
  LightTreeGenReco->Branch("gen_in_acc",&tree_gen_in_acc,"gen_in_acc/O");
  LightTreeGenReco->Branch("gen_in_acc_has_matched_reco",&tree_gen_in_acc_has_matched_reco,"gen_in_acc_has_matched_reco/O");
  LightTreeGenReco->Branch("gen_in_acc_has_no_matched_reco",&tree_gen_in_acc_has_no_matched_reco,"gen_in_acc_has_no_matched_reco/O");

  fHNumPU = new TH1F("NumPU_rew","NumPU_rew",50,0,50);
  fHNumPU_noweight = new TH1F("NumPU_noweight","NumPU_noweight",50,0,50);
  fHNumPUTrue = new TH1F("NumPUTrue_rew","NumPUTrue_rew",50,0,50);
  fHNumPUTrue_noweight = new TH1F("NumPUTrue_noweight","NumPUTrue_noweight",50,0,50);
  fHNumVtx = new TH1F("NumVtx_rew","NumVtx_rew",50,0,50);
	

	
  cout << "Trees and histos created" << endl;



}

void DiPhotonMiniTree::Analyze(){

  //  cout << endl << "-----------------------" << endl;

  float weight;
  if (!isdata) weight = GetPUWeight(fTR->PUnumInteractions);
  else weight=1;
  
  event_luminormfactor=AddWeight;

  event_fileuuid = uuid;

  event_run = fTR->Run;
  event_lumi = fTR->LumiSection;
  event_number = fTR->Event;

  if (!isdata) {
    fHNumPU->Fill(fTR->PUnumInteractions,weight);
    fHNumPU_noweight->Fill(fTR->PUnumInteractions);
    fHNumPUTrue->Fill(fTR->PUnumTrueInteractions,weight);
    fHNumPUTrue_noweight->Fill(fTR->PUnumTrueInteractions);
  }
  fHNumVtx->Fill(fTR->NVrtx,weight);

  //  return; // RUNNING FOR PU FILE


  // FILTERS TO FIX
  // scraping veto done at ntuplizer level
  if (!PassPrimaryVertexFilter()) return;
  if (!fTR->HBHENoiseFilterResult) return;
  if (!fTR->CSCTightHaloID) return;
  //  if (!fTR->EcalDeadTPFilterFlag) return;
  //  if (!fTR->RA2TrackingFailureFilterFlag) return;

  event_weight = weight;
  event_rho = fTR->Rho;
  event_sigma = fTR->Sigma;
  if (!isdata) {
    event_nPU = fTR->PUnumInteractions;
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

  std::vector<int> passing_selection[18];

  bool pass[18];
  int pass12_whoissiglike;

  for (int sel_cat=0; sel_cat<18; sel_cat++){

    if (sel_cat!=11) continue;

    if (isstep2 && sel_cat!=0) continue;

    if (sel_cat!=10 && !passtrigger) continue; // no trigger for Zmumu selection

    if (isdata){
      if (sel_cat>=11) continue;
      if (sel_cat>=9) continue;
      if (sel_cat==4 || sel_cat==5) continue;
    }
    else {
      if (sel_cat==9 || sel_cat==10) continue;
      if (sel_cat==4 || sel_cat==5) continue;
    }

    pass[sel_cat]=false;

    std::vector<int> passing;

    if (sel_cat==10){
      for (int i=0; i<fTR->NMus; i++){
        passing.push_back(i);
      }
      passing = MuonSelection(fTR,passing);
    }
    else {
      for (int i=0; i<fTR->NPhotons; i++){
	passing.push_back(i);
      }
      for (int i=0; i<(int)(passing.size())-1; i++){
	assert(fTR->PhoPt[passing.at(i)]>=fTR->PhoPt[passing.at(i+1)]);
      }
      passing = PhotonPreSelection(fTR,passing);
    }


    if (sel_cat==0){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==1){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==2){
      //      std::vector<int> passing_sel = PhotonSelection(fTR,passing);
      //      std::vector<int> passing_presel = passing;
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      //      pass[sel_cat] = (passing_sel.size()>=1) && SinglePhotonEventSelection(fTR,passing);
      //      pass[sel_cat] = (passing_presel.size()>=2) && SinglePhotonEventSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==3){
      passing = PhotonSelection(fTR,passing,"revert_pixel_veto"); // revert pixel veto already done
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==4){
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==5){
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==6){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==7){
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> passing_sel = PhotonSelection(fTR,passing);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sel.size()>=1){
	int fondo = passing_bkg[0];
	int forcone = passing_sel[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==8){
      passing = PhotonSelection(fTR,passing,"invert_sieie_cut");
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==9){
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==10){
      pass[sel_cat] = DiMuonFromZSelection(fTR,passing);
    }
    else if (sel_cat==11){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==12){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = SinglePhotonEventSelection(fTR,passing);
    }
    else if (sel_cat==13){
      passing = SignalSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==14){
      passing = PhotonSelection(fTR,passing);
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==15){
      passing = BackgroundSelection(fTR,passing);
      passing = PhotonSelection(fTR,passing);
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==16){
      std::vector<int> passing_bkg = BackgroundSelection(fTR,passing);
      passing_bkg = PhotonSelection(fTR,passing_bkg);
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      passing_sig = PhotonSelection(fTR,passing_sig);
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int forcone = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[forcone]) {newpassing.push_back(fondo); newpassing.push_back(forcone); pass12_whoissiglike=1;}
	else {newpassing.push_back(forcone); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);
    }
    else if (sel_cat==17){
      std::vector<int> passing_sig = SignalSelection(fTR,passing);
      passing_sig = PhotonSelection(fTR,passing_sig);
      std::vector<int> passing_bkg = PhotonSelection(fTR,passing,"invert_sieie_cut");
      std::vector<int> newpassing;
      if (passing_bkg.size()>=1 && passing_sig.size()>=1){
	int fondo = passing_bkg[0];
	int prompt = passing_sig[0];
	if (fTR->PhoPt[fondo] > fTR->PhoPt[prompt]) {newpassing.push_back(fondo); newpassing.push_back(prompt); pass12_whoissiglike=1;}
	else {newpassing.push_back(prompt); newpassing.push_back(fondo); pass12_whoissiglike=0;}
      }
      passing=newpassing;
      pass[sel_cat] = StandardEventSelection(fTR,passing);      
    }

    passing_selection[sel_cat] = passing;

  }



  //  cout << "D" << endl;

  for (int sel_cat=0; sel_cat<18; sel_cat++){

    if (sel_cat!=11) continue;

    if (isstep2 && sel_cat!=0) continue;

    if (sel_cat!=10 && !passtrigger) continue; // no trigger for Zmumu selection

    if (isdata){
      if (sel_cat>=11) continue;
      if (sel_cat>=9) continue;
      if (sel_cat==4 || sel_cat==5) continue;
    }
    else {
      if (sel_cat==9 || sel_cat==10) continue;
      if (sel_cat==4 || sel_cat==5) continue;
    }

    if (!pass[sel_cat]) continue;

    std::vector<int> passing = passing_selection[sel_cat];
    int minsize = (is2d[sel_cat]) ? 2 : 1;

    if (passing_selection[sel_cat].size()<minsize){
      std::cout << "Error!!!" << std::endl;
      continue;
    }

    if (sel_cat==10){
      for (int i=0; i<passing.size(); i++){
	ResetVars();
	FillMuonInfo(passing.at(i));
	OutputTree[sel_cat]->Fill();
      }
    }
    else if (is2d[sel_cat]){
      ResetVars();
      FillLead(passing.at(0));
      FillTrail(passing.at(1));
      bool dofill=true;

      if (sel_cat==7 || sel_cat==14 || sel_cat==16 || sel_cat==17) event_pass12whoissiglike=pass12_whoissiglike;

      if (sel_cat==6 || ((sel_cat==7 || sel_cat==16) && pass12_whoissiglike==0)) {

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
	    pholead_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(0)];
	    pholead_PhoSCRemovalPFIsoCombined = pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
	    if (pholead_PhoSCRemovalPFIsoCharged==999 || pholead_PhoSCRemovalPFIsoNeutral==999 || pholead_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	  }
      }
      if (sel_cat==6 || ((sel_cat==7 || sel_cat==16) && pass12_whoissiglike==1)) {

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
	  photrail_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(1)];
	  photrail_PhoSCRemovalPFIsoCombined = photrail_PhoSCRemovalPFIsoCharged+photrail_PhoSCRemovalPFIsoNeutral+photrail_PhoSCRemovalPFIsoPhoton;
	  if (photrail_PhoSCRemovalPFIsoCharged==999 || photrail_PhoSCRemovalPFIsoNeutral==999 || photrail_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	}

	}


      if (sel_cat==0 && isstep2){ // new templates from event mixing

	if (!inputtree_isinitialized) InitInputTree();

	Long64_t index_matchingtree = matchingtree->GetEntryNumberWithIndex(event_run*10000+event_lumi,event_number);
	if (index_matchingtree<0) {cout << "NO MATCHING FOUND!!!" << endl; cout << event_run << " " << event_lumi << " " << event_number << endl; dofill=false;}
	matchingtree->GetEntry(index_matchingtree);
	if (event_run!=matchingtree_event_run || event_lumi!=matchingtree_event_lumi || event_number!=matchingtree_event_number){cout << "WRONG MATCHING (including under/overflow)" << endl; dofill=false;}

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

      if (sel_cat==7) {
	std::vector<int> removals = GetPFCandInsideFootprint(fTR,passing.at(!pass12_whoissiglike),0,"photon");
	int index=0;
	for (int k=0; k<fTR->NPfCand; k++){
	  if (index==global_maxN_photonpfcandidates) {std::cout << "Too many pfcandidates" << std::endl; dofill=false; break;}
	  if (fTR->PfCandPdgId[k]!=22) continue;
	  float eta = fabs(fTR->PfCandEta[k]);
	  if (eta>1.4442 && eta<1.566) continue;
	  if (eta>2.5) continue;
	  if (fTR->PhoisPFPhoton[passing.at(!pass12_whoissiglike)] && fTR->PhoMatchedPFPhotonCand[passing.at(!pass12_whoissiglike)]==k) continue;	
	  bool removed = false;
	  for (int j=0; j<removals.size(); j++) if (k==removals.at(j)) removed=true;
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

	FillVetoObjects(fTR,passing.at(!pass12_whoissiglike),TString("exclude_object_itself"));

      }


      float invmass = (CorrPhoton(fTR,passing.at(0))+CorrPhoton(fTR,passing.at(1))).M();
      dipho_mgg_photon = invmass;
      if (dofill) OutputTree[sel_cat]->Fill();
    }

    else if (!is2d[sel_cat]){

      for (int i=0; i<passing.size(); i++){
      ResetVars();
      FillLead(passing.at(i));
      bool dofill = true;

      if (sel_cat==1) {


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
	pholead_PhoSCRemovalPFIsoCharged = fTR->PhoSCRemovalPFIsoChargedRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoNeutral = fTR->PhoSCRemovalPFIsoNeutralRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoPhoton = fTR->PhoSCRemovalPFIsoPhotonRCone[passing.at(i)];
	pholead_PhoSCRemovalPFIsoCombined = pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
	if (pholead_PhoSCRemovalPFIsoCharged==999 || pholead_PhoSCRemovalPFIsoNeutral==999 || pholead_PhoSCRemovalPFIsoPhoton==999) dofill=false;
	}

      }

//      if (sel_cat==11 || sel_cat==12) for (int k=0; k<50; k++) {
//	pholead_test_rotatedphotoniso[k]=PFIsolation(passing.at(i),0.025*k,"photon",NULL,NULL,NULL,NULL,NULL,NULL);
//	if (!FindCloseJetsAndPhotons(fTR,0.025*k,passing.at(i),"")) pholead_test_rotatedwithcheckphotoniso[k]=PFIsolation(passing.at(i),0.025*k,"photon",NULL,NULL,NULL,NULL,NULL,NULL);
//      }


      if (sel_cat==1 || sel_cat==2) {
	std::vector<int> removals = GetPFCandInsideFootprint(fTR,passing.at(i),0,"photon");
	int index=0;
	for (int k=0; k<fTR->NPfCand; k++){
	  //	  if (index==global_maxN_photonpfcandidates) {std::cout << "Too many pfcandidates" << std::endl; dofill=false; break;}
	  if (fTR->PfCandPdgId[k]!=22) continue;
	  //	  cout << fTR->PfCandPt[k] << " " << fTR->PfCandEta[k] << " " << fTR->PfCandPhi[k] << endl;
	  float eta = fabs(fTR->PfCandEta[k]);
	  if (eta>1.4442 && eta<1.566) continue;
	  if (eta>2.5) continue;
	  if (fTR->PhoisPFPhoton[passing.at(i)] && fTR->PhoMatchedPFPhotonCand[passing.at(i)]==k) continue;	
	  bool removed = false;
	  for (int j=0; j<removals.size(); j++) if (k==removals.at(j)) removed=true;
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
	allphotonpfcand_count = index;

	FillVetoObjects(fTR,passing.at(i),(sel_cat==1) ? TString("") : TString("exclude_object_itself"));

      }

      if (dofill) OutputTree[sel_cat]->Fill();
      }

    }

    
    }
  

  { // lightweight tree for efficiency and unfolding studies

    ResetVars(); 

    std::vector<int> passing;
    for (int i=0; i<fTR->NPhotons; i++){
      passing.push_back(i);
    }
    for (int i=0; i<(int)(passing.size())-1; i++){
      assert(fTR->PhoPt[passing.at(i)]>=fTR->PhoPt[passing.at(i+1)]);
    }

    std::vector<int> passing_gen;
    for (int i=0; i<fTR->NGenPhotons; i++){
      passing_gen.push_back(i);
    }

    bool found_reco = false;
    bool gen_in_acc = false;

    passing = PhotonPreSelection(fTR,passing);
    passing = PhotonSelection(fTR,passing);
    found_reco = passtrigger && StandardEventSelection(fTR,passing);

    bool reco_has_matched_gen_no_acceptance = false;
    bool reco_has_matched_gen_within_acceptance = false;
    bool reco_has_matched_gen_outside_acceptance = false;

    if (found_reco) {
 
      bool match0 = ( (fTR->PhoMCmatchexitcode[passing.at(0)]==1 || fTR->PhoMCmatchexitcode[passing.at(0)]==2) && fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[passing.at(0)]]<5 );
      bool match1 = ( (fTR->PhoMCmatchexitcode[passing.at(1)]==1 || fTR->PhoMCmatchexitcode[passing.at(1)]==2) && fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[passing.at(1)]]<5 ); 
      reco_has_matched_gen_no_acceptance = match0 && match1;

      if (reco_has_matched_gen_no_acceptance){
	int m0 = fTR->PhoMCmatchindex[passing.at(0)];
	int m1 = fTR->PhoMCmatchindex[passing.at(1)];
	bool match0acc = ( (fabs(fTR->GenPhotonEta[m0])<1.4442) || (fabs(fTR->GenPhotonEta[m0])>1.566 && fabs(fTR->GenPhotonEta[m0])<2.5) ) && (fTR->GenPhotonPt[m0]>25);
	bool match1acc = ( (fabs(fTR->GenPhotonEta[m1])<1.4442) || (fabs(fTR->GenPhotonEta[m1])>1.566 && fabs(fTR->GenPhotonEta[m1])<2.5) ) && (fTR->GenPhotonPt[m1]>25);
	reco_has_matched_gen_within_acceptance = match0acc && match1acc && ((fTR->GenPhotonPt[m0]>40) || (fTR->GenPhotonPt[m1]>40)) && (Util::GetDeltaR(fTR->GenPhotonEta[m0],fTR->GenPhotonEta[m1],fTR->GenPhotonPhi[m0],fTR->GenPhotonPhi[m1])>global_dR_cut_acceptance);
	reco_has_matched_gen_outside_acceptance = !reco_has_matched_gen_within_acceptance;
      }

    }

    for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); ){
      bool pass=0;
      if ( (fabs(fTR->GenPhotonEta[*it])<1.4442) || (fabs(fTR->GenPhotonEta[*it])>1.566 && fabs(fTR->GenPhotonEta[*it])<2.5)  ) pass=1;
      if (fTR->GenPhotonPt[*it]<=25) pass=0;
      if (fTR->GenPhotonIsoDR04[*it]>5) pass=0;
      if (!pass) it=passing_gen.erase(it); else it++;
    }

    for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); ){
      bool pass=0;
      int mother = fTR->GenPhotonMotherID[*it];
      if (mother>=-6 && mother<=6) pass=1;
      if (mother==21) pass=1;
      if (mother==22 && fTR->GenPhotonMotherStatus[*it]==3) pass=1;
      if (!pass) it=passing_gen.erase(it); else it++;
    }

    {
      int sizegenphotonsbefore = passing_gen.size();
      std::vector<OrderPair> passing_gen_ordered;
      for (vector<int>::iterator it = passing_gen.begin(); it != passing_gen.end(); it++){
	passing_gen_ordered.push_back(make_pair<int,float>(int(*it),float(fTR->GenPhotonPt[*it])));
      }
      std::sort(passing_gen_ordered.begin(),passing_gen_ordered.end(),indexComparator);
      passing_gen.clear();
      for (vector<OrderPair>::iterator it = passing_gen_ordered.begin(); it != passing_gen_ordered.end(); it++) passing_gen.push_back(it->first);
      assert(sizegenphotonsbefore==passing_gen.size());
    }
    for (int i=0; i<(int)(passing_gen.size())-1; i++){
      assert(fTR->GenPhotonPt[passing_gen.at(i)]>=fTR->GenPhotonPt[passing_gen.at(i+1)]);
    }

    if (passing_gen.size()>=2 && fTR->GenPhotonPt[passing_gen.at(0)]>40 && \
	Util::GetDeltaR(fTR->GenPhotonEta[passing_gen.at(0)],fTR->GenPhotonEta[passing_gen.at(1)],fTR->GenPhotonPhi[passing_gen.at(0)],fTR->GenPhotonPhi[passing_gen.at(1)])>global_dR_cut_acceptance) gen_in_acc=true;

    bool gen_in_acc_has_matched_reco = false;

    if (found_reco && gen_in_acc) {
      int found = 0;
      for (std::vector<int>::iterator it = passing.begin(); it != passing.end(); it++){
	if (fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)
	  if (fTR->PhoMCmatchindex[*it]==passing_gen.at(0) || fTR->PhoMCmatchindex[*it]==passing_gen.at(1)) found++;
      }
      if (found==2 && (fTR->PhoMCmatchindex[passing.at(0)] != fTR->PhoMCmatchindex[passing.at(1)])) gen_in_acc_has_matched_reco=true;
    }

    bool gen_in_acc_has_no_matched_reco = gen_in_acc && !gen_in_acc_has_matched_reco;

    //CHECK
    if (gen_in_acc && reco_has_matched_gen_no_acceptance){
      bool ok = false;
      if (passing_gen.at(0)==fTR->PhoMCmatchindex[passing.at(0)] && passing_gen.at(1)==fTR->PhoMCmatchindex[passing.at(1)]) ok = true;
      if (passing_gen.at(0)==fTR->PhoMCmatchindex[passing.at(1)] && passing_gen.at(1)==fTR->PhoMCmatchindex[passing.at(0)]) ok = true;
      if (!ok) {std::cout << "PROBLEM! MISMATCH IN GEN OBJECTS! Throwing away the event." << std::endl; reco_has_matched_gen_no_acceptance=false; gen_in_acc=false;}
    }

    if (reco_has_matched_gen_no_acceptance) {
      FillLead(passing.at(0));
      FillTrail(passing.at(1));
    }

    {
      int m0=-1; int m1=-1;
      if (gen_in_acc) {m0 = passing_gen.at(0); m1=passing_gen.at(1);}
      else if (reco_has_matched_gen_no_acceptance){
	m0 = fTR->PhoMCmatchindex[passing.at(0)];
	m1 = fTR->PhoMCmatchindex[passing.at(1)];
	if (fTR->GenPhotonPt[m0]<fTR->GenPhotonPt[m1]) {int temp=m1; m1=m0; m0=temp;}
      }
      if (gen_in_acc || reco_has_matched_gen_no_acceptance) {
	TLorentzVector genphotons[2];
	genphotons[0].SetPtEtaPhiM(fTR->GenPhotonPt[m0],fTR->GenPhotonEta[m0],fTR->GenPhotonPhi[m0],0);
	genphotons[1].SetPtEtaPhiM(fTR->GenPhotonPt[m1],fTR->GenPhotonEta[m1],fTR->GenPhotonPhi[m1],0);
	pholead_GEN_eta =     genphotons[0].Eta();
	photrail_GEN_eta =    genphotons[1].Eta();
	pholead_GEN_phi =     genphotons[0].Phi();
	photrail_GEN_phi =    genphotons[1].Phi();
	pholead_GEN_pt =      genphotons[0].Pt();
	photrail_GEN_pt =     genphotons[1].Pt();
      }
    }

    tree_reco_has_matched_gen_no_acceptance = reco_has_matched_gen_no_acceptance;
    tree_reco_has_matched_gen_within_acceptance = reco_has_matched_gen_within_acceptance;
    tree_reco_has_matched_gen_outside_acceptance = reco_has_matched_gen_outside_acceptance;
    tree_gen_in_acc = gen_in_acc;
    tree_gen_in_acc_has_matched_reco = gen_in_acc_has_matched_reco;
    tree_gen_in_acc_has_no_matched_reco = gen_in_acc_has_no_matched_reco;


    if (reco_has_matched_gen_no_acceptance || gen_in_acc) LightTreeGenReco->Fill();

  }




  
};

void DiPhotonMiniTree::End(){
  fOutputFile->cd();
  for (int i=0; i<18; i++) OutputTree[i]->Write();	
  LightTreeGenReco->Write();
  fHNumPU->Write();
  fHNumPU_noweight->Write();
  fHNumPUTrue->Write();
  fHNumPUTrue_noweight->Write();
  fHNumVtx->Write();
	
  fOutputFile->Close();

}


void DiPhotonMiniTree::FillPhoIso_NewTemplates(TreeReader *fTR, Int_t *n1_arr, Int_t *n2_arr, std::vector<int> passing, SigBkgMode mode, ChoiceMixingTemplates mixing){

  //  if (mixing!=k1Event || mode!=kSigBkg) return;

  //  cout << "EVENT " << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]] << " " << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]] << endl;

  int m1 = (mode==kSigSig || mode==kSigBkg) ? 0 : 1;
  int m2 = (mode==kSigSig || mode==kBkgSig) ? 0 : 1;
  
  int found = 0;

	for (int l=0; l<nclosest_inputmatching; l++){
	  if (found==nclosest) break;
	  int n1 = n1_arr[l];
	  int n2 = n2_arr[l];
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
	    if (InputTree[m1]->GetEntry(n1)<=0) continue;
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
	    
	    if (skip_EBEE) {
	      cout << "EB/EE migration, skipping1 " << m1 << " " << m2 << " " << l << endl; 
	      cout << fTR->Run << " " << fTR->LumiSection << " " << fTR->Event << endl;
	      cout << fTR->SCEta[fTR->PhotSCindex[passing.at(0)]] << " " << input_pholead_SCeta << endl;
	      continue;
	    }
	  }
	  
	  skip_EBEE = false;

	  if (n2>=0){
	    if (InputTree[m2]->GetEntry(n2)<=0) continue;
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
	    
	    if (skip_EBEE) {
	      cout << "EB/EE migration, skipping2 " << m1 << " " << m2 << " " << l << endl; 
	      cout << fTR->Run << " " << fTR->LumiSection << " " << fTR->Event << endl;
	      cout << fTR->SCEta[fTR->PhotSCindex[passing.at(1)]] << " " << input_pholead_SCeta << endl;
	      continue;
	    }
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


TLorentzVector DiPhotonMiniTree::CorrPhoton(TreeReader *fTR, int i){
  TLorentzVector corr;
  corr.SetPtEtaPhiE(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i],fTR->PhoEnergy[i]); 
  float corre = fTR->PhoRegrEnergy[i];
  corr.SetE(corre);
  corr.SetRho(corre);
  return corr;
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

std::vector<int> DiPhotonMiniTree::PhotonPreSelection(TreeReader *fTR, std::vector<int> passing){

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
    if (fabs(eta)<1.4442) energy*=phocorr->getEtaCorrectionBarrel(eta);
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
    float r9=R9Rescale(fTR->SCR9[fTR->PhotSCindex[*it]]);
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
    float sieie=SieieRescale(fTR->PhoSigmaIetaIeta[*it],(bool)(fabs(eta)<1.4442));
    bool pass=0;
    if (fabs(eta)<1.4442 && sieie<0.014 && sieie>0.001) pass=1; // to add sigmaiphiphi>0.001 in the future
    if (fabs(eta)>1.566 && sieie<0.034) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }
	
  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (trigger)
    float r9=R9Rescale(fTR->SCR9[fTR->PhotSCindex[*it]]);
    bool pass=0;
    float etcorrecaliso=fTR->PhoIso03Ecal[*it]-0.012*fTR->PhoPt[*it];
    float etcorrhcaliso=fTR->PhoIso03Hcal[*it]-0.005*fTR->PhoPt[*it];
    float etcorrtrkiso=fTR->PhoIso03TrkHollow[*it]-0.002*fTR->PhoPt[*it];
    if (r9<0.9 && etcorrecaliso<4 && etcorrhcaliso<4 && etcorrtrkiso<4) pass=1;
    if (r9>0.9 && etcorrecaliso<50 && etcorrhcaliso<50 && etcorrtrkiso<50) pass=1;
    if (!pass) it=passing.erase(it); else it++;
  }

//  XXX EM ENRICHMENT PRESELECTION NOT APPLIED
//  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){ // isolation cuts (filter)
//    bool pass=0;
//    if(fTR->pho_Cone02ChargedHadronIso_dR02_dz02_dxy01[*it]<4) pass=1;
//    if (!pass) it=passing.erase(it); else it++;
//  }

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

  // Genlevel isolation cut for MC
  // selects only photons matched to gen level, with gen-level iso < 5 GeV

  if (isdata) {
    std::cout << "Calling gen level isolation cut on data!!!" << std::endl;
    return passing;
  }
  
  for (std::vector<int>::iterator it = passing.begin(); it != passing.end(); ){
      bool pass=0;
      if (fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)
	if(fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[*it]]<5)
	  pass=1;
      if (!pass) it=passing.erase(it); else it++;
    }

  return passing;

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
    float sieie=SieieRescale(fTR->PhoSigmaIetaIeta[*it],(bool)(fabs(eta)<1.4442));
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
    if (fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2) pass=1;
    if (pass) {
      float dR = Util::GetDeltaR(fTR->GenPhotonEta[fTR->PhoMCmatchindex[*it]],fTR->PhoEta[*it],fTR->GenPhotonPhi[fTR->PhoMCmatchindex[*it]],fTR->PhoPhi[*it]);
      if (dR>0.1) cout << "PATHOLOGICAL MATCHING DR " << dR << endl;
    }
    if (!pass) it=passing.erase(it); else it++;
  }

  std::vector<int> passing2 = GenLevelIsolationCut(fTR,passing);

  return passing2;

};

std::vector<int> DiPhotonMiniTree::BackgroundSelection(TreeReader *fTR, std::vector<int> passing){

  for (vector<int>::iterator it = passing.begin(); it != passing.end(); ){
    bool pass=0;
    if (!(fTR->PhoMCmatchexitcode[*it]==1 || fTR->PhoMCmatchexitcode[*it]==2)) {
      pass=1;
    }
    else {
      if (fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[*it]]>5) pass=1;
    }
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

  float invmass0 = (CorrPhoton(fTR,passing.at(0))+CorrPhoton(fTR,passing.at(1))).M();
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

bool DiPhotonMiniTree::StandardEventSelection(TreeReader *fTR, std::vector<int> &passing){

  if (passing.size()<2) return false;

  passing.resize(2); // keep only the first two

  float invmass0 = (CorrPhoton(fTR,passing.at(0))+CorrPhoton(fTR,passing.at(1))).M();
  float deta=fTR->PhoEta[passing.at(0)]-fTR->PhoEta[passing.at(1)];
  float dphi=Util::DeltaPhi(fTR->PhoPhi[passing.at(0)],fTR->PhoPhi[passing.at(1)]);
  double dR=sqrt(dphi*dphi+deta*deta);

  if (fTR->PhoPt[passing.at(0)]<40) return false;
  if (fTR->PhoPt[passing.at(1)]<25) return false;
  //  if (invmass0<80) return false;
  if (dR<global_dR_cut_acceptance) return false;

  return true;
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
  
  const bool debug=false;
  if (debug) std::cout << "calling FindCloseJetsAndPhotons eta=" << eta << " phi=" << phi << std::endl;

  const float mindR = 0.8;
  bool found=false;

  for (int i=0; i<fTR->NJets; i++){
    if (fTR->JPt[i]<20) continue;
    //    cout << "vj " << fTR->JPt[i] << " " << fTR->JEta[i] << " " << fTR->JPhi[i] << endl;
    //    if (!(fTR->JPassPileupIDT0[fTR->JVrtxListStart[i]+0])) continue;
    float dR = Util::GetDeltaR(eta,fTR->JEta[i],phi,fTR->JPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    if (dR<mindR) found=true;
    if (debug) if (dR<mindR) std::cout << "Found jet eta=" << fTR->JEta[i] << " phi=" << fTR->JPhi[i] << std::endl;
  }

  for (int i=0; i<fTR->NPhotons; i++){
    if (fTR->PhoPt[i]<10) continue;
    //    cout << "vg " << fTR->PhoPt[i] << " " << fTR->PhoEta[i] << " " << fTR->PhoPhi[i] << endl;
    float dR = Util::GetDeltaR(eta,fTR->PhoEta[i],phi,fTR->PhoPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    if (dR<mindR) found=true;
    if (debug) if (dR<mindR) std::cout << "Found phot eta=" << fTR->PhoEta[i] << " phi=" << fTR->PhoPhi[i] << std::endl;
  }

  if (debug) std::cout << "returning " << found << std::endl;
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

  std::vector<std::pair<TVector3,int>> obj;

  if (mod!="" && mod!="exclude_object_itself") {std::cout << "error" << std::endl;}

  TVector3 photon_position = TVector3(fTR->SCX[fTR->PhotSCindex[phoqi]],fTR->SCY[fTR->PhotSCindex[phoqi]],fTR->SCZ[fTR->PhotSCindex[phoqi]]);

  double eta = photon_position.Eta();
  double phi = photon_position.Phi();
  
  for (int i=0; i<fTR->NJets; i++){
    if (fTR->JPt[i]<20) continue;
    float dR = Util::GetDeltaR(eta,fTR->JEta[i],phi,fTR->JPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    TVector3 a;
    a.SetPtEtaPhi(fTR->JPt[i],fTR->JEta[i],fTR->JPhi[i]);
    obj.push_back(std::pair<TVector3,int>(a,0));
  }

  for (int i=0; i<fTR->NPhotons; i++){
    if (fTR->PhoPt[i]<10) continue;
    float dR = Util::GetDeltaR(eta,fTR->PhoEta[i],phi,fTR->PhoPhi[i]);
    if (mod=="exclude_object_itself") if (dR<0.2) continue;
    TVector3 a;
    a.SetPtEtaPhi(fTR->PhoPt[i],fTR->PhoEta[i],fTR->PhoPhi[i]);
    obj.push_back(std::pair<TVector3,int>(a,1));
  }

  if (obj.size()>global_maxN_vetoobjects) {std::cout << "MaxN vetoobjects reached" << std::endl; obj.resize(global_maxN_vetoobjects);}
  for (int i=0; i<obj.size(); i++){
    vetoobjects_pt[i]=obj.at(i).first.Pt();
    vetoobjects_eta[i]=obj.at(i).first.Eta();
    vetoobjects_phi[i]=obj.at(i).first.Phi();
    vetoobjects_type[i]=obj.at(i).second;
  }
  vetoobjects_count = obj.size();

};

std::vector<int> DiPhotonMiniTree::GetPFCandIDedRemovals(TreeReader *fTR, int phoqi){
  std::vector<int> out;
  if (fTR->PhoisPFPhoton[phoqi]) out.push_back(fTR->PhoMatchedPFPhotonCand[phoqi]);
  return out;
};

std::vector<int> DiPhotonMiniTree::GetPFCandInsideFootprint(TreeReader *fTR, int phoqi, float rotation_phi, TString component){
  return GetPFCandWithFootprintRemoval(fTR,phoqi,rotation_phi,false,component);
};

std::vector<int> DiPhotonMiniTree::GetPFCandInsideFootprint(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, TString component){
  return GetPFCandWithFootprintRemoval(fTR,pfcands,phoqi,rotation_phi,false,component);
};

std::vector<int> DiPhotonMiniTree::GetPFCandWithFootprintRemoval(TreeReader *fTR, int phoqi, float rotation_phi, bool outoffootprint, TString component){

  if (component!="photon"){
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return std::vector<int>();
  }

  if (component!="neutral" && component!="charged" && component!="photon" && component!="combined") {
    std::cout << "Wrong choice for component" << std::endl;
    return std::vector<int>();
  }

  int scindex = fTR->PhotSCindex[phoqi];
  
  if (scindex<0) {
    std::cout << "Error in GetPFCandOverlappingSC" << std::endl;
    std::cout << scindex << " " << phoqi << std::endl;
    return std::vector<int>();
  }

  bool isbarrel = fTR->PhoisEB[phoqi];
  int nxtals = fTR->SCNXtals[scindex];

  std::vector<int> result;

  for (int i=0; i<fTR->NPfCand; i++){

    int type = FindPFCandType(fTR->PfCandPdgId[i]);
    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    TVector3 sc_position = TVector3(fTR->SCX[scindex],fTR->SCY[scindex],fTR->SCZ[scindex]);

    bool inside=false;

    if (fTR->PhoisPFPhoton[phoqi] && fTR->PhoMatchedPFPhotonCand[phoqi]==i) continue;

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
    if (inside) result.push_back(i);

  }

  return result;

};

std::vector<int> DiPhotonMiniTree::GetPFCandWithFootprintRemoval(TreeReader *fTR, pfcandidates_struct *pfcands, int phoqi, float rotation_phi, bool outoffootprint, TString component){

  assert(rotation_phi==0);

  if (component!="photon"){
    std::cout << "propagation not implemented for non photon objects!!!" << std::endl;
    return std::vector<int>();
  }

  if (phoqi<0) return std::vector<int>();

  int scindex = fTR->PhotSCindex[phoqi];
  
  if (scindex<0) {
    std::cout << "Error in GetPFCandOverlappingSC" << std::endl;
    std::cout << scindex << " " << phoqi << std::endl;
    return std::vector<int>();
  }

  bool isbarrel = fTR->PhoisEB[phoqi];
  int nxtals = fTR->SCNXtals[scindex];

  std::vector<int> result;

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
    if (inside) result.push_back(i);

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
  
  std::vector<int> footprint;
  if (doremoval1) {std::vector<int> footprint1 = GetPFCandInsideFootprint(fTR,pfcands,phoqi1,0,"photon"); footprint.insert(footprint.end(),footprint1.begin(),footprint1.end());} // phoqi<0 ritorna vuoto
  if (doremoval2) {std::vector<int> footprint2 = GetPFCandInsideFootprint(fTR,pfcands,phoqi2,0,"photon"); footprint.insert(footprint.end(),footprint2.begin(),footprint2.end());} // phoqi<0 ritorna vuoto

  for (int i=0; i<pfcands->PfCandPt.size(); i++){

    bool removed=false;
    for (int j=0; j<footprint.size(); j++) {
      if (i==footprint.at(j)) removed=true;
    }
    if (removed) continue;

    if (phoqi1>=0) if (GetPFCandDeltaRFromSC(fTR,pfcands,phoqi1,i,matched_eta1).dR<0.4) result1+=pfcands->PfCandPt[i];
    if (phoqi2>=0) if (GetPFCandDeltaRFromSC(fTR,pfcands,phoqi2,i,matched_eta2).dR<0.4) result2+=pfcands->PfCandPt[i];
    

  } // end pf cand loop


  return std::make_pair<float,float>(float(result1),float(result2));

};

float DiPhotonMiniTree::PFIsolation(int phoqi, float rotation_phi, TString component, int *counter, std::vector<float> *energies, std::vector<float> *ets, std::vector<float> *detas, std::vector<float> *dphis, float* newphi, std::vector<int> removals){

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
    std::vector<int> footprint = GetPFCandInsideFootprint(fTR,phoqi,rotation_phi,"photon");
    for (int i=0; i<footprint.size(); i++) removals.push_back(footprint.at(i));
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

    if (fTR->PhoisPFPhoton[phoqi] && fTR->PhoMatchedPFPhotonCand[phoqi]==i) continue;
 
    int type = FindPFCandType(fTR->PfCandPdgId[i]);

    if (!(type==0 || type==1 || type==2)) continue;

    if (component=="neutral" && type!=0) continue;
    if (component=="charged" && type!=1) continue;
    if (component=="photon" && type!=2) continue;

    bool removed = false;
    for (int j=0; j<removals.size(); j++) {
      if (i==removals.at(j)) removed=true;
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


void DiPhotonMiniTree::FillLead(int index){

  pholead_eta = fTR->PhoEta[index];
  pholead_phi = fTR->PhoPhi[index];
  pholead_pt = fTR->PhoPt[index];
  pholead_energy = fTR->PhoEnergy[index];
  pholead_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  pholead_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  pholead_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  pholead_r9 = R9Rescale(fTR->SCR9[fTR->PhotSCindex[index]]);
  pholead_sieie = SieieRescale(fTR->PhoSigmaIetaIeta[index],(bool)(fabs(fTR->SCEta[fTR->PhotSCindex[index]])<1.4442));
  pholead_hoe = fTR->PhoHoverE[index];

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
    pholead_PhoSCRemovalPFIsoCharged=fTR->PhoSCRemovalPFIsoCharged[index];
    pholead_PhoSCRemovalPFIsoNeutral=fTR->PhoSCRemovalPFIsoNeutral[index];
    pholead_PhoSCRemovalPFIsoPhoton=fTR->PhoSCRemovalPFIsoPhoton[index];
    pholead_PhoSCRemovalPFIsoCombined=pholead_PhoSCRemovalPFIsoCharged+pholead_PhoSCRemovalPFIsoNeutral+pholead_PhoSCRemovalPFIsoPhoton;
  }
  pholead_PhoPassConversionVeto=fTR->PhoPassConversionVeto[index];
  if (fTR->PhoMCmatchindex[index]>=0) pholead_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
  pholead_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
//  pholead_scarea = scarea[fTR->PhotSCindex[index]];
//  pholead_scareaSF = scareaSF[fTR->PhotSCindex[index]];

  {
  jetmatching_struct m = PFMatchPhotonToJet(index);
  cout << "Matching: " << m.m_jet << endl;
  cout << "Footprint/RECO: " << m.phopt_footprint_total/pholead_pt << endl;
  cout << "Fraction of footprint: " << m.phopt_footprint_m_frac << endl;
  cout << "Fraction of jet pt: " << m.jetpt_m_frac << endl;
  cout << "Pho: " << pholead_pt << " " << pholead_eta << " " << pholead_phi << endl;
  if (m.m_jet>=0) cout << "Jet: " << m.jetpt_pf << " " <<  fTR->JEta[m.m_jet] << " " << fTR->JPhi[m.m_jet] << endl;
  }

};

float DiPhotonMiniTree::SieieRescale(float sieie, bool isbarrel){
  if (isdata) return sieie; // rescale sieie only in MC
  return isbarrel ? 0.87*sieie+0.0011 : 0.99*sieie;
};

float DiPhotonMiniTree::R9Rescale(float r9){
  if (isdata) return r9; // rescale r9 only in MC
  return 1.005*r9;
};

void DiPhotonMiniTree::FillTrail(int index){

  photrail_eta = fTR->PhoEta[index];
  photrail_phi = fTR->PhoPhi[index];
  photrail_pt = fTR->PhoPt[index];
  photrail_energy = fTR->PhoEnergy[index];
  photrail_SCeta = fTR->SCEta[fTR->PhotSCindex[index]];
  photrail_SCphi = fTR->SCPhi[fTR->PhotSCindex[index]];
  photrail_PhoHasPixSeed=fTR->PhoHasPixSeed[index];
  photrail_r9 = R9Rescale(fTR->SCR9[fTR->PhotSCindex[index]]);
  photrail_sieie = SieieRescale(fTR->PhoSigmaIetaIeta[index],(bool)(fabs(fTR->SCEta[fTR->PhotSCindex[index]])<1.4442));
  photrail_hoe = fTR->PhoHoverE[index];

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
    photrail_PhoSCRemovalPFIsoCharged=fTR->PhoSCRemovalPFIsoCharged[index];
    photrail_PhoSCRemovalPFIsoNeutral=fTR->PhoSCRemovalPFIsoNeutral[index];
    photrail_PhoSCRemovalPFIsoPhoton=fTR->PhoSCRemovalPFIsoPhoton[index];
    photrail_PhoSCRemovalPFIsoCombined=photrail_PhoSCRemovalPFIsoCharged+photrail_PhoSCRemovalPFIsoNeutral+photrail_PhoSCRemovalPFIsoPhoton;
  }

  photrail_PhoPassConversionVeto=fTR->PhoPassConversionVeto[index];
  if (fTR->PhoMCmatchindex[index]>=0) photrail_GenPhotonIsoDR04=fTR->GenPhotonIsoDR04[fTR->PhoMCmatchindex[index]];
  photrail_PhoMCmatchexitcode=fTR->PhoMCmatchexitcode[index];
//  photrail_scarea = scarea[fTR->PhotSCindex[index]];
//  photrail_scareaSF = scareaSF[fTR->PhotSCindex[index]];

};

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
  tree_reco_has_matched_gen_no_acceptance = false;
  tree_reco_has_matched_gen_within_acceptance = false;
  tree_reco_has_matched_gen_outside_acceptance = false;
  tree_gen_in_acc = false;
  tree_gen_in_acc_has_matched_reco = false;
  tree_gen_in_acc_has_no_matched_reco = false;

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


bool DiPhotonMiniTree::FindImpingingTrack(TreeReader *fTR, int phoqi, int &reference_index_found, bool dofootprintremoval, std::vector<int> removals){

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
//    if (fTR->PhoisPFPhoton[phoqi] && fTR->PhoMatchedPFPhotonCand[phoqi]==i) continue;
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

int DiPhotonMiniTree::CountChargedHadronsInCone(TreeReader *fTR, int phoqi, std::vector<int> removals, bool skipvetocones){

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

std::vector<int> DiPhotonMiniTree::NChargedHadronsInConeSelection(TreeReader *fTR, std::vector<int> passing, int minimum, int maximum){

  return std::vector<int>();

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

void DiPhotonMiniTree::BugfixWrongGenPhotonsMatching(TreeReader *fTR){   // redo the matching to fix bug dphi (no fabs in matching)

  for (int i=0; i<fTR->NPhotons; i++){

    if (fTR->PhoMCmatchexitcode[i]==-1 || fTR->PhoMCmatchexitcode[i]==0) continue;

    int matched=-1;
    for (int k=0; k<fTR->NGenPhotons; k++){
      if (fabs(fTR->GenPhotonPt[k]-fTR->PhoPt[i])/fTR->GenPhotonPt[k]>2) continue;
      if (Util::GetDeltaR(fTR->GenPhotonEta[k],fTR->PhoEta[i],fTR->GenPhotonPhi[k],fTR->PhoPhi[i])>0.1) continue;
      matched=k; break;
    }

    if (matched==-1) {fTR->PhoMCmatchexitcode[i]=-2; fTR->PhoMCmatchindex[i]=-999;}
    else {

      fTR->PhoMCmatchindex[i]=matched;

      int mother = fTR->GenPhotonMotherID[matched];      
      if ((mother>=-6 && mother<=6) || (mother==21)) fTR->PhoMCmatchexitcode[i]=1;
      else if (mother==22 && fTR->GenPhotonMotherStatus[matched]==3) fTR->PhoMCmatchexitcode[i]=2;
      else fTR->PhoMCmatchexitcode[i]=3;


    }
  }

};

void DiPhotonMiniTree::InitInputTree(){
  if (isstep2){
    TString title = Form("%s/matchingtree_%u.root",input_filename.Data(),uuid);
    cout << "opening " << title.Data() << endl;
    f_input = TFile::Open(title.Data(),"read");
    f_input->GetObject("Tree_1Drandomcone_template",InputTree[0]);
    f_input->GetObject("Tree_2Drandomconesideband_template",InputTree[1]);
    for (int m=0; m<2; m++){
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

    matchingtree->BuildIndex("matchingtree_event_run*10000+matchingtree_event_lumi","matchingtree_event_number");

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

  // prepare list of pfcands to represent the photon deposit
  std::vector<int> pfcands = GetPFCandInsideFootprint(fTR,phoqi,0,"photon");
  if (fTR->PhoisPFPhoton[phoqi] && fTR->PhoMatchedPFPhotonCand[phoqi]>=0) {
    int m = fTR->PhoMatchedPFPhotonCand[phoqi];
    for (int i=0; i<pfcands.size(); i++) assert(pfcands.at(i)!=m); // this should never happen
    pfcands.push_back(m);
  }

  // init ranking
  std::vector<std::pair<int,float> > ranking;
  for (int i=0; i<fTR->NJets; i++) ranking.push_back(std::pair<int,float>(i,0));
  if (ranking.size()==0) {
    cout << "PFMatchPhotonToJet: no jets in the event! Returning error state" << endl;
    return out;
  }

  // loop on candidates
  out.phopt_footprint_total=0;
  for (int i=0; i<pfcands.size(); i++){
    float pt = fTR->PfCandPt.at(pfcands.at(i));
    out.phopt_footprint_total+=pt;
    int j = fTR->PfCandBelongsToJet.at(pfcands.at(i));
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
