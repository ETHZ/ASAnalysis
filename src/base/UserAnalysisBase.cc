#include <cstdlib>

#include <TH1I.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TTree.h>

#include "helper/pdgparticle.hh"
#include "helper/Monitor.hh"
#include "base/TreeReader.hh"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "base/UserAnalysisBase.hh"


using namespace std;

UserAnalysisBase::UserAnalysisBase(TreeReader *tr, bool isData, string globaltag){
  fTR = tr;
  fTlat = new TLatex();
  fVerbose = 0;
  fDoPileUpReweight = false;


  if(globaltag == "") globaltag = "START53_V7A"; // need a default GT, otherwise the next line crashes
  fMetCorrector = new OnTheFlyCorrections(globaltag, isData);
}

UserAnalysisBase::~UserAnalysisBase(){
  if(fDoPileUpReweight) delete fPUWeight;
  delete fMetCorrector;

}
std::pair<float ,float> UserAnalysisBase::GetOnTheFlyCorrections(){
  float corrMetx = fTR->PFMETpx;
  float corrMety = fTR->PFMETpy;
  for (int ind = 0; ind<fTR->JMetCorrRawEta.size(); ++ind) {

    std::pair <float, float> corr = fMetCorrector->getCorrections(fTR->JMetCorrRawPt[ind],fTR->JMetCorrRawEta[ind],fTR->JMetCorrNoMuPt[ind],
                                                                  fTR->JMetCorrPhi[ind],fTR->JMetCorrEMF[ind],fTR->Rho, fTR->JMetCorrArea[ind]);
    corrMetx += corr.first ;
    corrMety += corr.second;
  }
  float newmet    = sqrt(corrMetx*corrMetx + corrMety*corrMety);
  float newmetphi = atan2(corrMety, corrMetx);
  return make_pair(newmet, newmetphi);
  
}

void UserAnalysisBase::BeginRun(Int_t& run) {
  // Called when run changes
  // Need to re-create HLT map in case it changed
    
  GetHLTNames();
}

void UserAnalysisBase::ReadPDGTable(const char* filename){
  // Fills the fPDGMap map from a textfile to associate pdgids with names
  int pdgid(0), type(0);
  string Name, Texname, Typename;
  ifstream IN(filename);
  char buff1[200], buff2[200], buff3[200];
  char readbuff[200];
  // Loop over lines of datafile
  while( IN.getline(readbuff, 200, '\n') ){
    if (readbuff[0] == '#') {continue;} // Skip lines commented with '#'
    sscanf(readbuff, "%d %s %d %s %s", &type, buff1, &pdgid, buff2, buff3);
    // Convert chararrays to strings
    Typename = string(buff1); Name = string(buff2); Texname = string(buff3);
    pdgparticle *p = new pdgparticle(pdgid, Name, Texname, type, Typename);
    // Fill map
    fPDGMap[pdgid] = *p;
  }
}

int UserAnalysisBase::GetPDGParticle(pdgparticle &part, int id){
  if( fPDGMap.empty() ){
    if(fVerbose > 0) cout << "UserAnalysisBase::GetPDGParticle ==> PDGMap not filled!" << endl;
    return -1;
  }
  else{
    map<int, pdgparticle>::iterator it = fPDGMap.find(id);
    if(it == fPDGMap.end()){
      if(fVerbose > 1) cout << "UserAnalysisBase::GetPDGParticle ==> PDGParticle with ID " << id << " not found!" << endl;
      return -1;
    }
    else{
      part = it->second;
      return 1;
    }
  }
  return 0;
}

void UserAnalysisBase::GetHLTNames(){
	
  fHLTLabelMap.clear();
  fHLTLabels.clear();
  fHLTLabels.reserve(fTR->HLTNames.size());
  for( size_t i=0; i < fTR->HLTNames.size(); i++ ){
    fHLTLabelMap[fTR->HLTNames[i]] = i; 
    fHLTLabels.push_back(fTR->HLTNames[i]);
    if (fVerbose>3) cout << " " << i << " " << fTR->HLTNames[i] << endl;
  }
}

int UserAnalysisBase::GetHLTBit(string theHltName){
  if( fHLTLabelMap.empty() ) return -1;
  else{
    map<string,int>::iterator it = fHLTLabelMap.find(theHltName);
    if(it == fHLTLabelMap.end()){
      if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTBit ==> Bit with name " << theHltName << " not found!" << endl;
      return -1;
    }
    else{
      return it->second;
    }
  }
}

bool UserAnalysisBase::GetHLTResult(string theHltName){
  if( fHLTLabelMap.empty() ) return false;
  else{
    int bit = GetHLTBit(theHltName);
    if (bit < 0 ) {
      if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTResult ==> Bit with name " << theHltName << " not found!" << endl;
      return false;
    }
    else return (bool)fTR->HLTResults[bit];
  }
}

int UserAnalysisBase::GetHLTPrescale(string theHltName) {
  if( fHLTLabelMap.empty() ) return 0;
  else{
    int bit = GetHLTBit(theHltName);
    if (bit < 0 ) {
      if(fVerbose > 1) cout << "UserAnalysisBase::GetHLTPrescale ==> Bit with name " << theHltName << " not found!" << endl;
      return 0;
    }
    else return fTR->HLTPrescale[bit];
  }
}

int UserAnalysisBase::getSusyMass(int pdgid, int round){
  float mpart(0.);
  int npart(0);
  for (int i = 0; i < fTR->nGenParticles; ++i){
    if (fabs(fTR->genInfoId[i]) != pdgid || fTR->genInfoStatus[i] != 3) continue;
    else {
      mpart += fTR->genInfoM[i];
      npart++;
    }
  }
  mpart = mpart/npart; // averaging over both particles in the event
  int roundmass = round * (int) (mpart/round + 0.5);
  return roundmass;
}

float UserAnalysisBase::getISRWeight(float susyPt, int flag){ // flag==0: central flag==1: up flag==2: down
  if ( susyPt <= 120.){ 
    return 1.0;
  }
  else if ( susyPt > 120. && susyPt <= 150.){ 
    if (flag == 0) return 0.95;
    if (flag == 1) return 1.00;
    if (flag == 2) return 0.90;
  }
  else if ( susyPt > 150. && susyPt <= 250.){ 
    if (flag == 0) return 0.90;
    if (flag == 1) return 1.00;
    if (flag == 2) return 0.80;
  }
  else if ( susyPt > 250.                  ){ 
    if (flag == 0) return 0.80;
    if (flag == 1) return 1.00;
    if (flag == 2) return 0.60;
  }
  else {
    cout << "SOMETHING WENT WRONG IN THE ISR SYSTEMATIC!! APPARENTLY YOU HAVE A NEGATIVE SUSY-PT VALUE..?" << endl;
    exit(-1);
  }
  return -9999999.9;
}

int UserAnalysisBase::getNParticle(int pdgid, int status){
	int npart(0);
	for (int i = 0; i < fTR->nGenParticles; ++i){
		if (fabs(fTR->genInfoId[i]) != pdgid || fTR->genInfoStatus[i] != status) continue;
		else {
			npart++;
		}
	}
	return npart;
}

float UserAnalysisBase::getSusySystemPt(int pdgid1, int pdgid2){
	int npart(0);
	TLorentzVector tmp, final;
	if (pdgid2 == -1) { // quick fix for ewino. if pdgid2 isn't set, just search for pair-produced particles with pdgid1 as before
		for (int i = 0; i < fTR->nGenParticles; ++i){
			if (fabs(fTR->genInfoId[i]) != pdgid1 || fTR->genInfoStatus[i] != 3) continue;
			else {
				npart++;
				tmp.SetPtEtaPhiM(fTR->genInfoPt[i], fTR->genInfoEta[i], fTR->genInfoPhi[i], fTR->genInfoM[i]);
				final += tmp;
			}
		}
	}
	else {
		for (int i = 0; i < fTR->nGenParticles; ++i){
			if ((fabs(fTR->genInfoId[i]) == pdgid1 && fTR->genInfoStatus[i] == 3 ) ||
				(fabs(fTR->genInfoId[i]) == pdgid2 && fTR->genInfoStatus[i] == 3 )  ) {
				npart++;
				tmp.SetPtEtaPhiM(fTR->genInfoPt[i], fTR->genInfoEta[i], fTR->genInfoPhi[i], fTR->genInfoM[i]);
				final += tmp;
			}
		}
	}

	if (npart != 2) {
		cout << " MORE OR LESS THAN TWO OF YOUR DESIRED SUSY INITIAL PARTICLES FOUND!! CHECK UP ON THAT!!" << endl;
		exit(-10);
	}
	return final.Pt();
}

///////////////////////////////////////////////////////////////
// Object selections:
// JETS

bool UserAnalysisBase::IsGoodBasicPFJet(int index){
  // Basic PF jet cleaning and ID cuts
  // Loose PF jet ID (WARNING: HF not included in our ntuple)
  // See PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
  if ( !(fTR->JNConstituents[index] > 1) )    return false;
  if ( !(fTR->JNeutralEmFrac[index]     < 0.99) ) return false;
  if ( !(fTR->JNeutralHadFrac[index]    < 0.99) ) return false;
  if (fabs(fTR->JEta[index]) < 2.4 ) { // Cuts for |eta|<2.4
    if ( !(fTR->JChargedEmFrac[index]  < 0.99) )  return false;
    if ( !(fTR->JChargedHadFrac[index] > 0.00) )  return false;
    if ( !(fTR->JNAssoTracks[index]    > 0   ) )  return false;
  }
  return true;
}

bool UserAnalysisBase::IsGoodPFJetMedium(int index) {
  // Medium PF JID
  if ( ! IsGoodBasicPFJet(index) ) return false;
  if ( !(fTR->JNeutralHadFrac[index] < 0.95)         ) return false;
  if ( !(fTR->JNeutralEmFrac[index]  < 0.95)         ) return false;
  return true;
}

bool UserAnalysisBase::IsGoodPFJetTight(int index) {
  // Tight PF JID
  if ( ! IsGoodBasicPFJet(index) ) return false;
  if ( !(fTR->JNeutralHadFrac[index] < 0.90)         ) return false;
  if ( !(fTR->JNeutralEmFrac[index]  < 0.90)         ) return false;
  return true;
}

// MUONS
float UserAnalysisBase::MuPFIso(int index){
  double neutral = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index] - 0.5*fTR->MuPfIsoR03SumPUPt[index] );
  float iso = ( fTR->MuPfIsoR03ChHad[index] + TMath::Max(0., neutral) ) / fTR->MuPt[index];
  return iso;
}

float UserAnalysisBase::MuPFIso04(int index){
  double neutral = (fTR->MuPfIsoR04NeHad[index] + fTR->MuPfIsoR04Photon[index] - 0.5*fTR->MuPfIsoR04SumPUPt[index] );
  float iso = ( fTR->MuPfIsoR04ChHad[index] + TMath::Max(0., neutral) ) / fTR->MuPt[index];
  return iso;
}

float UserAnalysisBase::MuRadIso(int index){
  return (fTR->MumuonRadPFIsoChHad03[index] + fTR->MumuonRadPFIsoNHad03[index] + fTR->MumuonRadPFIsoPhoton03[index]) / fTR->MuPt[index] ;
}

// ELECTRONS

bool UserAnalysisBase::ElPassesConvRej(int index){
  if (fTR->ElNumberOfMissingInnerHits[index] > 0 ) return false;
  if (!fTR->ElPassConversionVeto[index] ) return false;
  return true;
}

float UserAnalysisBase::Aeff(float eta) {
  float abseta = fabs(eta); // making sure we're looking at |eta|
  // now for HCP dataset for cone of dR < 0.3
  if(abseta < 1.0)   return 0.13;
  if(abseta < 1.479) return 0.14;
  if(abseta < 2.0)   return 0.07;
  if(abseta < 2.2)   return 0.09;
  if(abseta < 2.3)   return 0.11;
  if(abseta < 2.4)   return 0.11;
  return 0.14;
}

float UserAnalysisBase::ElPFIso(int index){
  double neutral = fTR->ElEventelPFIsoValueNeutral03PFIdStandard[index] + fTR->ElEventelPFIsoValueGamma03PFIdStandard[index];
  double rhocorr = fTR->RhoForIso * Aeff(fTR->ElSCEta[index]);
  double iso = ( fTR->ElEventelPFIsoValueCharged03PFIdStandard[index] + TMath::Max(0., neutral - rhocorr) )/ fTR->ElPt[index];
  return iso;
}

float UserAnalysisBase::ElRadIso(int index){
  double iso = ( fTR->ElelectronRadPFIsoChHad03[index] + fTR->ElelectronRadPFIsoNHad03[index] + fTR->ElelectronRadPFIsoPhoton03[index] ) / fTR->ElPt[index];
  return iso;
}

float UserAnalysisBase::ElRelDetIso(int index){
  // Apply 1 GeV subtraction in ECAL Barrel isolation
  if( fabs(fTR->ElEta[index]) < 1.479 ){ // Barrel
    double iso = ( fTR->ElDR03TkSumPt[index] + TMath::Max(0., fTR->ElDR03EcalRecHitSumEt[index] - 1.) + fTR->ElDR03HcalTowerSumEt[index] ) / fTR->ElPt[index];
    return iso;
  }
  return fTR->ElRelIso03[index]; // EndCap
}

///////////////////////////////////////////////////////////////
// Event selections:
bool UserAnalysisBase::IsGoodEvent(){
  // Some cuts on the primary vertex
  double pvx = fTR->PrimVtxx;
  double pvy = fTR->PrimVtxy;
  double pvz = fTR->PrimVtxz;
  if (fTR->PrimVtxIsFake) return false;
  if (fabs(pvz) > 24.) return false;
  if (sqrt(pvx*pvx + pvy*pvy) > 2.0) return false; // Wrt 0,0
  if (fTR->PrimVtxNdof < 4.0) return false; // this is a float, cut value at 4.0
  return true;
}

///////////////////////////////////////////////////////////////
// JEC
float UserAnalysisBase::getNewJetInfo(int ind, string which){
  std::vector<float> newjet = fMetCorrector->getCorrPtECorr(fTR->JPt[ind]   , fTR->JEta[ind]  , fTR->JE[ind], 
                                                            fTR->JEcorr[ind], fTR->JArea[ind] , fTR->Rho);
  if (which == "pt"  ) return newjet[0]; // new jet pt
  if (which == "e"   ) return newjet[1]; // new jet energy
  if (which == "corr") return newjet[2]; // new jet correction
}
float UserAnalysisBase::GetJetPtNoResidual(int ind){
  float pt = fMetCorrector->getJetPtNoResidual(fTR->JPt[ind], fTR->JEta[ind], fTR->JEcorr[ind], 
                                               fTR->JArea[ind], fTR->Rho);
  return pt;
}
float UserAnalysisBase::GetJECUncert(float pt, float eta){
  return fMetCorrector->getJECUncertainty(pt, eta);
}


// ---------------------------------------------
// Pile Up Reweighting
void UserAnalysisBase::SetPileUpSrc(string data_PileUp, string mc_PileUp){
  if(data_PileUp.size() == 0 || mc_PileUp.size() == 0 ) return;
  if(fDoPileUpReweight  == 1 ) {cout << "ERROR in SetPileUpSrc: fPUWeight already initialized" << endl; return; }
    
  TFile *fpu = new TFile(data_PileUp.c_str());
  string pileupUp="pileupUp";
  string pileupDown="pileupDown";
  if(!((TH1D*)fpu->Get("pileupUp"))||!((TH1D*)fpu->Get("pileupDown"))) {
    //the user supplied a file with only the standard histogram, which does not contain the up and down shapes
    cout << "Did not find up and down shapes in provided file " << endl;
    pileupUp="pileup";
    pileupDown="pileup";
  }
  fpu->Close();
    
  fPUWeight     = new reweight::LumiReWeighting( mc_PileUp, data_PileUp,"pileup","pileup");
  fPUWeightUp   = new reweight::LumiReWeighting( mc_PileUp, data_PileUp,"pileup",pileupUp);
  fPUWeightDown = new reweight::LumiReWeighting( mc_PileUp, data_PileUp,"pileup",pileupDown);
  fDoPileUpReweight = true;
}


//2012 PU reweighting: USE NUMBER OF TRUE INTERACTIONS!
float UserAnalysisBase::GetPUWeight(float nPUTrueinteractions){
  if(! fDoPileUpReweight) return -999.99;
  else return fPUWeight->weight( nPUTrueinteractions); 
}

float UserAnalysisBase::GetPUWeightUp(float nPUTrueinteractions){
  if(! fDoPileUpReweight) return -999.99;
  else return fPUWeightUp->weight( nPUTrueinteractions); 
}

float UserAnalysisBase::GetPUWeightDown(float nPUTrueinteractions){
  if(! fDoPileUpReweight) return -999.99;
  else return fPUWeightDown->weight( nPUTrueinteractions); 
}

