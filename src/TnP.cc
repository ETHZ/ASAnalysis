#include "TnP.hh"

#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TRandom3.h"

#include "TWbox.h"
#include "TMultiGraph.h"
#include "TGaxis.h"


#include "TF1.h"
#include "TH1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TStyle.h"

// stuff for fitting
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"


#include <iostream>
#include <iomanip>
#include <time.h> // access to date/time

#include <typeinfo>


using namespace std;
using namespace RooFit;

int TnP::Wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0; 
}


//____________________________________________________________________________
TnP::TnP(TString inputfile, bool createHistos){
	fCreateHistos = createHistos;
	fInputFile = inputfile;

	fOutputFileName = "tnp_output.root";

	checkFlavor();
	if (!fIsData) 
		fPUWeight = new reweight::LumiReWeighting("/shome/mdunser/puhistos/MC2012PU.root", "/shome/mdunser/puhistos/may21/DataPUTrue_may21_upDown.root", "pileup", "pileup");

	if(fIsMu){
		fIsoCut = 0.10;
		fD0Cut  = 0.005;
	}
	else{
		fIsoCut = 0.09;
		fD0Cut  = 0.01;
	}
	
	effs = new eff[fnBins];
	// dynamic creation of the arrays
	fPassIDoverAll   = new TH1F*[fnBins];
	fPassIPoverID    = new TH1F*[fnBins];
	fPassISOoverIDIP = new TH1F*[fnBins];

	fFailIDoverAll   = new TH1F*[fnBins];
	fFailIPoverID    = new TH1F*[fnBins];
	fFailISOoverIDIP = new TH1F*[fnBins];

	if(fCreateHistos){
		for(int i=0; i<fnBins; ++i){
			fPassIDoverAll   [i] = new TH1F(Form("PassIDoverAll_%d"  ,i), Form("PassIDoverAll_%d"  ,i), 150, 0., 150.);
			fPassIPoverID    [i] = new TH1F(Form("PassIPoverID_%d"   ,i), Form("PassIPoverID_%d"   ,i), 150, 0., 150.);
			fPassISOoverIDIP [i] = new TH1F(Form("PassISOoverIDIP_%d",i), Form("PassISOoverIDIP_%d",i), 150, 0., 150.);
                                                                                                  
			fFailIDoverAll   [i] = new TH1F(Form("FailIDoverAll_%d"  ,i), Form("FailIDoverAll_%d"  ,i), 150, 0., 150.);
			fFailIPoverID    [i] = new TH1F(Form("FailIPoverID_%d"   ,i), Form("FailIPoverID_%d"   ,i), 150, 0., 150.);
			fFailISOoverIDIP [i] = new TH1F(Form("FailISOoverIDIP_%d",i), Form("FailISOoverIDIP_%d",i), 150, 0., 150.);

			fPassIDoverAll   [i] ->Sumw2();
			fPassIPoverID    [i] ->Sumw2();
			fPassISOoverIDIP [i] ->Sumw2();

			fFailIDoverAll   [i] ->Sumw2();
			fFailIPoverID    [i] ->Sumw2();
			fFailISOoverIDIP [i] ->Sumw2();
		}
	} // end if createhistos
	else{
		TFile *inFile = TFile::Open(fInputFile, "READ");
		if(inFile == NULL){
			cout << "File " << fInputFile << " does not exist!" << endl;
			exit(1);
		}
		else{
			std::cout << "reading file " + fInputFile << std::endl;
		}
		inFile->cd();

		for(int i=0; i<fnBins; ++i){
			if(fVerbose > 2) inFile->ls();
			getObjectSafe(inFile, Form("PassIDoverAll_%d"  ,i), fPassIDoverAll   [i]);
			getObjectSafe(inFile, Form("PassIPoverID_%d"   ,i), fPassIPoverID    [i]);
			getObjectSafe(inFile, Form("PassISOoverIDIP_%d",i), fPassISOoverIDIP [i]);

			getObjectSafe(inFile, Form("FailIDoverAll_%d"  ,i), fFailIDoverAll   [i]);
			getObjectSafe(inFile, Form("FailIPoverID_%d"   ,i), fFailIPoverID    [i]);
			getObjectSafe(inFile, Form("FailISOoverIDIP_%d",i), fFailISOoverIDIP [i]);
		}
	} // end else createhistos
}

void TnP::checkFlavor(){
	// check in the first event of the tree whether it's muon or electron tree
	TFile *pFile = TFile::Open(fInputFile);
	TTree *tree; getObjectSafe(pFile, "probeTree", tree);
	int isMu;
	tree->SetBranchAddress("isMuEvent", &isMu);
	tree->GetEntry(1);
	if(isMu == 1) {
		fIsMu  = true;
		fnBins = 13; // set the proper nbins
	}
	else {
		fIsMu  = false;
		fnBins = 25; // set the proper nbins
	}
	TString flav = fIsMu ? "muons" : "electrons";
	cout << "We're dealing with "+flav+" here!!" << endl;


	//checking also if it is data
	int run;
	tree->SetBranchAddress("Run", &run);
	tree->GetEntry(1);
	if(run==1){
	  fIsData = false;
	}else{
	  fIsData = true;
	}

	pFile->Close();
	delete pFile, tree;

}

//____________________________________________________________________________
TnP::~TnP(){
}

void TnP::doFitting(){

	if(fCreateHistos) {
		cout << "i'm going to do create all the histograms first, that will take a while... " << endl;
		fillHistos();
		writeHistos();
	}

	for(int i=0; i<fnBins; ++i){
	  //if(i>0) continue;
	  //if(i!=3) continue;
	  simFitPassFail(fPassIDoverAll  [i], fFailIDoverAll  [i], 0, i); // ID
	  simFitPassFail(fPassIPoverID   [i], fFailIPoverID   [i], 1, i); // IP
	  simFitPassFail(fPassISOoverIDIP[i], fFailISOoverIDIP[i], 2, i); // ISO
	}

	if (fIsMu) printMuTable();
	else       printElTable();

}
void TnP::simFitPassFail(TH1F* passHisto, TH1F* failHisto, int flag, int bin){

  if(passHisto == NULL) cout << "histogram doesn't exist!!!!!" << endl;

  double rangeFitMin(50.), rangeFitMax(140.);
  if(bin%6 == 0 ) rangeFitMin = 40.;
  if(bin%6 == 1 ) rangeFitMin = 45.;
  if(!fIsMu && bin == 19 ) rangeFitMin = 50.;
  if(!fIsMu && bin == 14 ) rangeFitMin = 40.;
  if(!fIsMu && bin == 18 ) rangeFitMin = 35.;
  if(!fIsMu && fIsData && bin == 7)    rangeFitMin = 40.;
  if(fIsMu && fIsData && bin == 0)    rangeFitMin = 35.;
  if(fIsMu && fIsData && bin == 0 && flag ==0 )    rangeFitMin = 40.;
  if(fIsMu && fIsData && bin == 1 )    rangeFitMin = 55.;
  if(fIsMu && fIsData && bin == 6 )    rangeFitMin = 30.;
  if(fIsMu && fIsData && bin == 3 )    {rangeFitMin = 55.;rangeFitMax = 120.;}
  RooRealVar px("px", "px", rangeFitMin,rangeFitMax) ;

    


  //Defining Double CB for SIGNAL
  RooRealVar pmean ("pmean" , "pass mean of gaussian" , 90, 85 , 95);
  RooRealVar psigma("psigma", "pass width of gaussian",  2, 0. ,  5); 
  RooRealVar a("a","a",2.,1.0,10);
  RooRealVar n("n","n",2.,0.5,10.);   
  RooRealVar aDx("aDx","aDx",3.,3,10.);
  RooRealVar nDx("nDx","nDx",5.,0.,10.);   
  RooDoubleCB func1("cb","cb PDF", px, pmean, psigma, a, n, aDx, nDx) ;


  //Defining Double CB for BACKGROUND
  RooRealVar mean_b ("mean_b" , "pass mean of gaussian" , 60, rangeFitMin*0.90 , rangeFitMax*1.5);
  RooRealVar sigma_b("sigma_b", "pass width of gaussian",  10, 2. ,  50); 
  RooRealVar a_b("a_b","a_b",4.,4.0,10);
  RooRealVar n_b("n_b","n_b",2.,0.5,10.);   
  RooRealVar aDx_b("aDx_b","aDx_b",2.,0.3,10.);
  RooRealVar nDx_b("nDx_b","nDx_b",5.,0.,20.);   
  RooDoubleCB cbBkg("cbBkg","cbBkg", px, mean_b, sigma_b, a_b, n_b, aDx_b, nDx_b) ;


  // --- this depend on the pT binning/cut on the probe
  if(bin%6 == 0) mean_b.setVal(45.);
  if(bin%6 == 1) mean_b.setVal(55.);
  if(bin%6 == 2) mean_b.setVal(70.);
  if(bin%6 == 3) mean_b.setVal(85.);
  if(bin%6 == 4) mean_b.setVal(90.);
  if(bin%6 == 5) mean_b.setVal(95.);
  // ---


  //Defining polynomial
  //RooRealVar a0("a0", "", 0., -1., 1.);
  //RooRealVar a1("a1", "", 0., -1., 1.);
  //RooRealVar a2("a2", "",  0., -1., 1.);
  //RooChebychev pol("pol","pol",px,RooArgList(a1,a2));
  RooRealVar a1("a1", "", 60.,-100,100);
  RooRealVar a2("a2", "", -1., -10., 0.);
  RooPolynomial pol("pol","pol",px,RooArgList(a1,a2));


  //Defining exponential function
  RooRealVar lambda("lambda", "slope", -0.1, -5., 0.);
  RooExponential func3("expo", "exponential PDF", px, lambda);


  //Defining gauss function  
  //RooRealVar gmean("gmean","gmean",60.,50.,80.);
  RooRealVar gmean("gmean","gmean",60.,50.,180.);
  RooRealVar gsigma("gsigma","gsigma",10.,5.,50.);   

  RooGaussian pgauss("pgauss", "pgaussian PDF", px, gmean, gsigma) ;

  //Defining sum of signal+bkg pdf functions
  //double initialSigValue = passHisto->Integral(rangeFitMin,rangeFitMax);
  double initialSigValue = passHisto->Integral(passHisto->FindBin(rangeFitMin),passHisto->FindBin(rangeFitMax));
  RooRealVar s("s", "signal yield", initialSigValue, initialSigValue/10.,initialSigValue*10.);
  RooRealVar b("b", "background yield", initialSigValue/10.,initialSigValue/10000.,initialSigValue*1.);  
  RooRealVar b1("b1", "background1 yield", initialSigValue*0.5,initialSigValue/10000.,initialSigValue*1.);
  RooRealVar b2("b2", "background2 yield", initialSigValue*0.25,initialSigValue/10000.,initialSigValue*1.);


  //RooAddPdf sum("sum", "DoubleCB plus Bkg",
  //		RooArgList(func1, func3), RooArgList(s, b));

  //RooAddPdf sum("sum", "as for failing probes",
  //		RooArgList(func1, pgauss, func3), RooArgList(s, b1,b2));

  RooAddPdf sum("sum", "DoubleCB plus Bkg",
  		RooArgList(func1, cbBkg), RooArgList(s, b));


  RooDataHist passData("passData", "passData", px, Import(*passHisto) );
  RooDataHist failData("failData", "failData", px, Import(*failHisto) );
  RooPlot* passFrame = px.frame(Title("pass TH1 with Poisson error bars"));
  RooPlot* failFrame = px.frame(Title("fail TH1 with Poisson error bars"));
  passData.plotOn(passFrame);
  failData.plotOn(failFrame);


  //func3.fitTo (passData );
  //func3.plotOn(passFrame);

  sum.fitTo (passData );
  sum.plotOn(passFrame);
  sum.plotOn(passFrame,RooFit::Components(cbBkg),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
  //sum.plotOn(passFrame,RooFit::Components(func3),RooFit::LineStyle(kDashed));
  //sum.plotOn(passFrame,RooFit::Components(pgauss),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
  sum.plotOn(passFrame,RooFit::Components(func1),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));


  float npass    = s.getVal();
  float npassErr = s.getError();
  


  // drawing things
  TLatex *lat = new TLatex();
  lat->SetNDC(kTRUE);
  lat->SetTextColor(kBlack);
  lat->SetTextSize(0.05);

  TString indText;
  if     (flag==0) indText = "ID over all";
  else if(flag==1) indText = "IP over ID";
  else if(flag==2) indText = "ISO over ID+IP";
  
  //TCanvas* c = new TCanvas("foo", "bar", 800, 675);
  //c->Divide(1,2);
  //c->cd(1);
  TCanvas* c = new TCanvas("foo", "bar", 800, 800);
  c->Divide(1,2);
  c->cd(1);
  
  //------------ Fitting failing-probe dilepton pairs
  bool fixSignalParameters=true;
  if(fixSignalParameters){
    if(fIsMu && fIsData && (bin == 8 || bin == 3) && flag == 0 )
      a.setRange(a.getVal()*0.60,a.getVal()*1.5);
    else
      a.setRange(a.getVal()*0.75,a.getVal()*1.5);
    n.setRange(n.getVal()*0.75,n.getVal()*1.5);
    aDx.setRange(aDx.getVal()*0.97,aDx.getVal()*1.03);
    nDx.setRange(nDx.getVal()*0.97,nDx.getVal()*1.03);
    pmean.setRange(pmean.getVal()*0.97,pmean.getVal()*1.03);
    psigma.setRange(psigma.getVal()*0.97,psigma.getVal()*1.3);
  }

  //px.setRange(50.,80.);
  initialSigValue = failHisto->Integral(failHisto->FindBin(rangeFitMin),failHisto->FindBin(rangeFitMax));
  //RooRealVar b1("b1", "background1 yield", initialSigValue*0.25,initialSigValue/10000.,initialSigValue*1.);
  //RooRealVar b2("b2", "background2 yield", initialSigValue*0.7,initialSigValue/10000.,initialSigValue*1.);
  b1.setVal(initialSigValue*0.25); b1.setRange(initialSigValue/10000.,initialSigValue*1.);
  b2.setVal(initialSigValue*0.7); b2.setRange(initialSigValue/10000.,initialSigValue*1.);
  b.setVal(initialSigValue*0.10); b.setRange(initialSigValue/10000.,initialSigValue*1.);
  s.setVal(initialSigValue/2.); s.setRange(initialSigValue/100.,initialSigValue*1.);
  //b1.setRange(initialSigValue/10000.,initialSigValue/2.);
  //b2.setRange(initialSigValue/10000.,initialSigValue/2.);

  //RooAddPdf failPdf("sum", "DoubleCB plus Poly plus Exp",
  //		    RooArgList(pol, func3), RooArgList(b1, b2));


  if(fIsMu && (bin==5 || bin==11)) {
    //gmean.setVal(120.); gmean.setRange(100.,150.);
    //gsigma.setVal(30.); gsigma.setRange(20.,50.);
    b1.setVal(initialSigValue*0.5);
    b.setVal(initialSigValue*0.5);
    s.setRange(initialSigValue/200.,initialSigValue*1.);
    pmean.setRange(pmean.getVal()*1.01,pmean.getVal()*0.99);

  }

  //
  if((!fIsMu) && (bin==0 || bin==1 || bin==6 || bin==7 ||  bin==12 || bin==13 || bin==18 || bin==19) ){
    b.setVal(initialSigValue*0.8);    
    s.setVal(initialSigValue*0.05);    
    s.setRange(initialSigValue*0.01,initialSigValue*0.70);    
    //cout << "initialSigValue, s, b: " << initialSigValue << " , " << s.getVal() << " , " << b.getVal() << endl;
    //Wait();
  }
  
  if(!fIsMu && bin == 24)  mean_b.setVal(55.);
  if(fIsMu && fIsData && bin == 11 && flag == 0 )    {
    mean_b.setRange(90.,110.);
    sigma_b.setRange(0.,15.);
  }
  if(fIsMu && fIsData && bin == 3 && flag == 0 )    {
    sigma_b.setRange(0.,8.);
  }


  //RooAddPdf failPdf("sum", "DoubleCB plus Bkg1 plus Bkg2",
  // 		    RooArgList(func1, pgauss, func3), RooArgList(s, b1, b2));

  RooAddPdf failPdf("sum", "DoubleCB plus Bkg",
  		    RooArgList(func1, cbBkg), RooArgList(s, b));



  passFrame->Draw();
  lat->DrawLatex(0.10,0.92, indText+" passing");
  lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));

  //cout << "Here I am " << endl;
  //gPad->Update(); Wait();
  
  c->cd(2);
  failPdf.fitTo (failData );
  failPdf.plotOn(failFrame);
  failPdf.plotOn(failFrame,RooFit::Components(cbBkg),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
  //failPdf.plotOn(failFrame,RooFit::Components(func3),RooFit::LineStyle(kDashed));
  //failPdf.plotOn(failFrame,RooFit::Components(pgauss),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
  failPdf.plotOn(failFrame,RooFit::Components(func1),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));

  float nfail    = s.getVal();
  float nfailErr = s.getError();

  float eff = (npass+nfail) > 0 ? npass/(npass+nfail) : 0.;
  //float effErr = (npass+nfail) > 0 ? npass/(npass+nfail)*TMath::Sqrt( (npassErr/npass)*(npassErr/npass) + (nfailErr/nfail)*(nfailErr/nfail) ) : 0.;

  //binomial uncertainty:
  float effErr = sqrt( npass/(npass+nfail)*(1- npass/(npass+nfail)) / (npass+nfail) );
  // float effErr = (npass+nfail) > 0 ? npass/(npass+nfail)*TMath::Sqrt( (npassErr/npass)*(npassErr/npass) + (nfailErr/nfail)*(nfailErr/nfail) ) : 0.;

  /*
  cout << "cb mean,sigma, a,n aD, nD: "
       << pmean.getVal() << " , " 
       << psigma.getVal() << " , " 
       << a.getVal() << " , " 
       << n.getVal() << " , "
       << aDx.getVal() << " , " 
       << nDx.getVal() << endl;


  cout << "sig, bkg: "
       << s.getVal() << " , "
       << b.getVal() << endl;

  Wait();
  */

  failFrame->Draw();
  lat->DrawLatex(0.10,0.92, indText+" failing");
  lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));


  // save as root and pdf
  c->Print(Form("tnp/passAndFail_bin%d_flag%d.pdf" ,bin,flag));
  c->Print(Form("tnp/passAndFail_bin%d_flag%d.png" ,bin,flag));
  //c->Print(Form("passAndFail_bin%d_flag%d.root",bin,flag));
  
  delete passFrame, c, lat;
  delete failFrame;


  // store the efficiencies:
  effs[bin].bin = bin;
  if(flag == 0) {
    effs[bin].idEff = eff;
    effs[bin].idEffErr = effErr;
  }
  if(flag == 1) {
    effs[bin].ipEff = eff;
    effs[bin].ipEffErr = effErr;
  }
  if(flag == 2) {
    effs[bin].isoEff = eff;
    effs[bin].isoEffErr = effErr;
  };
  

}


void TnP::fillHistos(){

	// load the tree

	TFile *pFile = TFile::Open(fInputFile);
	TTree *tree; getObjectSafe(pFile, "probeTree", tree);
	float tagPt, tagEta, tagPhi, tagIsoRel, tagD0, tagDz;
	float tagPassID;
	float probePt, probeEta, probePhi, probeIsoRel, probeD0, probeDz;
	float probePassID;
	float mll, deltaR, nvtx, pfmet;// , rho, ht, nTrueInt;
	int ntrue;
	// tag variables
	tree->SetBranchAddress("tagPt"       , &tagPt       );
	tree->SetBranchAddress("tagEta"      , &tagEta      );
	tree->SetBranchAddress("tagPhi"      , &tagPhi      );
	tree->SetBranchAddress("tagIsoRel"   , &tagIsoRel   );
	tree->SetBranchAddress("tagD0"       , &tagD0       );
	tree->SetBranchAddress("tagDz"       , &tagDz       );
	tree->SetBranchAddress("tagPassID"   , &tagPassID   );
	// probe variables
	tree->SetBranchAddress("probePt"     , &probePt     );
	tree->SetBranchAddress("probeEta"    , &probeEta    );
	tree->SetBranchAddress("probePhi"    , &probePhi    );
	tree->SetBranchAddress("probeIsoRel" , &probeIsoRel );
	tree->SetBranchAddress("probeD0"     , &probeD0     );
	tree->SetBranchAddress("probeDz"     , &probeDz     );
	tree->SetBranchAddress("probePassID" , &probePassID );
	// other variables
	tree->SetBranchAddress("mass2L"      , &mll         );
	tree->SetBranchAddress("deltaR"      , &deltaR      );
	tree->SetBranchAddress("nvtx"        , &nvtx        );
	tree->SetBranchAddress("pfmet"       , &pfmet       );
	tree->SetBranchAddress("nTrueInt"    , &ntrue       );

	long nentries = tree->GetEntries();
	cout << "looping over " << nentries << " events" << endl;

	int index = -1;	
	
	float puWeight(1.);

	for( int i = 0; i < nentries; i++ ){
		tree->GetEntry(i);
		if(!fIsData) puWeight = fPUWeight->weight(ntrue); // PU weight only for MC
		if(!passesAll(tagIsoRel, tagD0, tagDz, tagPassID)) continue; // make sure the tag passes everything
		index = fIsMu ? getPtEtaIndexMu(probePt, probeEta) : getPtEtaIndexEl(probePt, probeEta);
		if(!fIsMu and !checkElEta(probeEta)) continue;
		if(probePassID == 0){                               // probe fails ID
			fFailIDoverAll[index]   ->Fill(mll, puWeight);            // fill the histo
			fFailIDoverAll[fnBins-1]->Fill(mll, puWeight);            // fill the inclusive histo
		}
		else {                                              // probe passes ID
			fPassIDoverAll[index]   ->Fill(mll, puWeight);            // fill the histo
			fPassIDoverAll[fnBins-1]->Fill(mll, puWeight);            // fill the inclusive histo
			if(!passesIP(probeD0, probeDz)){                // probe fails IP but passes ID
				fFailIPoverID[index]   ->Fill(mll, puWeight);         // fill the histo
				fFailIPoverID[fnBins-1]->Fill(mll, puWeight);         // fill the inclusive histo
			}
			else {                                          // probe passes ID+IP
				fPassIPoverID[index]   ->Fill(mll, puWeight);         // fill the histo
				fPassIPoverID[fnBins-1]->Fill(mll, puWeight);         // fill the inclusive histo
				if(!passesIso(probeIsoRel)){                // probe fails ISO but passes ID+IP
					fFailISOoverIDIP[index]   ->Fill(mll, puWeight);  // fill the histo
					fFailISOoverIDIP[fnBins-1]->Fill(mll, puWeight);  // fill the inclusive histo
				}
				else{                                       // probe passes ID+IP+ISO
					fPassISOoverIDIP[index]   ->Fill(mll, puWeight);  // fill the histo
					fPassISOoverIDIP[fnBins-1]->Fill(mll, puWeight);  // fill the inclusive histo
				}
			}
		}
	}
	pFile->Close();

}
bool TnP::passesAll(float iso, float d0, float dz, int passid){
	if(iso       > fIsoCut) return false;
	if(fabs(d0)  > fD0Cut ) return false;
	if(fabs(dz)  > 0.10  ) return false;
	if(passid   == 0     ) return false;
	return true;
}
bool TnP::passesID(int passid){
	if(passid   == 0) return false;
	return true;
}
bool TnP::passesIP(float d0, float dz){
	if(fabs(d0)  > fD0Cut) return false;
	if(fabs(dz)  > 0.10 ) return false;
	return true;
}
bool TnP::passesIso(float iso){
	if(iso       > fIsoCut) return false;
	return true;
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________

void TnP::writeHistos(){
	
	TFile *pFile = TFile::Open(fInputFile, "UPDATE");
	pFile->cd();
	for(int i=0; i<fnBins; ++i){
		fPassIDoverAll   [i]->Write(fPassIDoverAll   [i]->GetName(), TObject::kOverwrite);
		fFailIDoverAll   [i]->Write(fFailIDoverAll   [i]->GetName(), TObject::kOverwrite);

		fPassIPoverID    [i]->Write(fPassIPoverID    [i]->GetName(), TObject::kOverwrite);
		fFailIPoverID    [i]->Write(fFailIPoverID    [i]->GetName(), TObject::kOverwrite);

		fPassISOoverIDIP [i]->Write(fPassISOoverIDIP [i]->GetName(), TObject::kOverwrite);
		fFailISOoverIDIP [i]->Write(fFailISOoverIDIP [i]->GetName(), TObject::kOverwrite);
	}
	pFile->Close();
}
void TnP::bookHistos(){
}
void TnP::deleteHistos(){
}
int  TnP::readHistos(TString filename){
	std::cout << "this is the filename of which to read the histos: " << filename << std::endl;
	return 0;
}

// =====================================================
// MUON FUNCTIONS
// =====================================================

int TnP::getPtEtaIndexMu(float pt, float eta){
	float aeta = fabs(eta);
	if(pt > 10. and pt <= 15.){
		if(aeta <= 1.2) return 0;
		if(aeta >  1.2) return 6;
	}
	else if(pt > 15. and pt <= 20.){
		if(aeta <= 1.2) return 1;
		if(aeta >  1.2) return 7;
	}
	else if(pt > 20. and pt <= 30.){
		if(aeta <= 1.2) return 2;
		if(aeta >  1.2) return 8;
	}
	else if(pt > 30. and pt <= 40.){
		if(aeta <= 1.2) return 3;
		if(aeta >  1.2) return 9;
	}
	else if(pt > 40. and pt <= 50.){
		if(aeta <= 1.2) return 4;
		if(aeta >  1.2) return 10;
	}
	else if(pt > 50.){
		if(aeta <= 1.2) return 5;
		if(aeta >  1.2) return 11;
	}
	else{
		cout << " ======== ERROR ================= " << endl;
		cout << " something went wrong with the index return for your muon... " << endl;
		cout << " exiting... " << endl;
		exit(-1);
	}
	return -1;
}
TString TnP::getPtEtaFromIndexMu(int i){
	TString str;
	switch (i) {
		case  0: str = "10 < pT < 15 GeV  ; |eta| <= 1.2"; break;
		case  6: str = "10 < pT < 15 GeV  ; |eta| >  1.2"; break;
		case  1: str = "15 < pT < 20 GeV  ; |eta| <= 1.2"; break;
		case  7: str = "15 < pT < 20 GeV  ; |eta| >  1.2"; break;
		case  2: str = "20 < pT < 30 GeV  ; |eta| <= 1.2"; break;
		case  8: str = "20 < pT < 30 GeV  ; |eta| >  1.2"; break;
		case  3: str = "30 < pT < 40 GeV  ; |eta| <= 1.2"; break;
		case  9: str = "30 < pT < 40 GeV  ; |eta| >  1.2"; break;
		case  4: str = "40 < pT < 50 GeV  ; |eta| <= 1.2"; break;
		case 10: str = "40 < pT < 50 GeV  ; |eta| >  1.2"; break;
		case  5: str = "     pT > 50 GeV  ; |eta| <= 1.2"; break;
		case 11: str = "     pT > 50 GeV  ; |eta| >  1.2"; break;
		case 12: str = "    inclusive                   "; break;
	}
	return str;
}
void TnP::printMuTable(){
  TFile *oFile = TFile::Open(fOutputFileName, "UPDATE");
  oFile->cd();
  
  float xbins[7]={10,15,20,30,40,50,80};
  
  string dataString;
  if(fIsData) dataString="Data"; else dataString="Sim";
  
  string flvString;
  if(fIsMu) flvString="Mu"; else flvString="El";
  

  string nameID  = flvString+dataString+"EffID";
  string nameIP  = flvString+dataString+"EffIP";
  string nameISO = flvString+dataString+"EffISO";
  string nameALL = flvString+dataString+"EffAll";
  
  TH1F* hEta1ID  = new TH1F(("eta1"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta1IP  = new TH1F(("eta1"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta1ISO = new TH1F(("eta1"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta1ALL = new TH1F(("eta1"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

  TH1F* hEta2ID  = new TH1F(("eta2"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta2IP  = new TH1F(("eta2"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta2ISO = new TH1F(("eta2"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta2ALL = new TH1F(("eta2"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

    
  for(int i=0; i<fnBins; ++i){
    if( i  <= 5){ 
      hEta1ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta1ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta1IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta1IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta1ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta1ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta1ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta1ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    


    if( i > 5 && i <= 11){ 
      hEta2ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta2ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta2IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta2IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta2ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta2ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta2ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta2ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    


  };
  hEta1ID  -> SetMaximum(1.05) ; hEta1ID  -> SetMinimum(0.20) ; hEta1ID  -> Write(hEta1ID  -> GetName(),TObject::kOverwrite)   ;
  hEta1IP  -> SetMaximum(1.05) ; hEta1IP  -> SetMinimum(0.20) ; hEta1IP  -> Write(hEta1IP  -> GetName(),TObject::kOverwrite)   ;
  hEta1ISO -> SetMaximum(1.05) ; hEta1ISO -> SetMinimum(0.20) ; hEta1ISO -> Write(hEta1ISO -> GetName(),TObject::kOverwrite) ;
  hEta1ALL -> SetMaximum(1.05) ; hEta1ALL -> SetMinimum(0.20) ; hEta1ALL -> Write(hEta1ALL -> GetName(),TObject::kOverwrite) ;

  hEta2ID  -> SetMaximum(1.05) ; hEta2ID  -> SetMinimum(0.20) ; hEta2ID  -> Write(hEta2ID  -> GetName(),TObject::kOverwrite)   ;
  hEta2IP  -> SetMaximum(1.05) ; hEta2IP  -> SetMinimum(0.20) ; hEta2IP  -> Write(hEta2IP  -> GetName(),TObject::kOverwrite)   ;
  hEta2ISO -> SetMaximum(1.05) ; hEta2ISO -> SetMinimum(0.20) ; hEta2ISO -> Write(hEta2ISO -> GetName(),TObject::kOverwrite) ;
  hEta2ALL -> SetMaximum(1.05) ; hEta2ALL -> SetMinimum(0.20) ; hEta2ALL -> Write(hEta2ALL -> GetName(),TObject::kOverwrite) ;

  oFile->Close();

  printText();

}
void TnP::printText(){
  ofstream OUT(fIsMu ? "muEfficiencies.txt" : "elEfficiencies.txt", ios::trunc);
  OUT << Form("========================================================") << endl;
  OUT << Form("========================================================") << endl;
  for(int i=0; i<fnBins; ++i){
    OUT << Form("At bin %2d: "+getPtEtaFromIndexMu(i), i+1) << endl;
    OUT << Form("========================") << endl;
    OUT << Form("ID-eff : %.3f +- %.3f", effs[i].idEff , effs[i].idEffErr ) << endl;
    OUT << Form("IP-eff : %.3f +- %.3f", effs[i].ipEff , effs[i].ipEffErr ) << endl;
    OUT << Form("ISO-eff: %.3f +- %.3f", effs[i].isoEff, effs[i].isoEffErr) << endl;
    OUT << Form("=============================================") << endl;
  }
  OUT.close();

  ofstream OUTTEX(fIsMu ? "muEfficiencies.tex" : "elEfficiencies.tex", ios::trunc);
  OUTTEX << Form("========================================================") << endl;
  OUTTEX << Form("                     10 < pT < 15   | 15 < pT < 20   & 20 < pT < 30     &    30 < pT 40  &   40 < pT < 50     & 50 < pT ") << endl;
  OUTTEX << Form("-------------------------------------------------------------------------------------------------------------------------") << endl;
  //OUT << Form("    |eta| < 1.2 : %.2f +- %.2f    |   %.2f +- %.2f   |%.2f +- %.2f    |%.2f +- %.2f    |%.2f +- %.2f    |%.2f +- %.2f    |")
  //OUT << Form("    |eta| > 1.2 : %.2f +- %.2f    |   %.2f +- %.2f   |%.2f +- %.2f    |%.2f +- %.2f    |%.2f +- %.2f    |%.2f +- %.2f    |")

  OUTTEX << Form("========================================================") << endl;

 
}
// =====================================================
// ELECTRON FUNCTIONS
// =====================================================
bool TnP::checkElEta(float eta){
	float aeta = fabs(eta);
	if(aeta > 2.4) return false;
	if(aeta > 1.4442 and aeta < 1.566) return false;
	return true;
}

int TnP::getPtEtaIndexEl(float pt, float eta){
	float aeta = fabs(eta);
	// cout << Form("pt: %.2f   eta: %.2f", pt, aeta) << endl;
	if     (pt >= 10. and pt < 15.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  0;
		if( 0.8   <= aeta and aeta <  1.4442) return  6;
		if( 1.566 <= aeta and aeta <  2.0   ) return 12;
		if( 2.0   <= aeta and aeta <  2.4   ) return 18;
	}
	else if(pt >= 15. and pt < 20.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  1;
		if( 0.8   <= aeta and aeta <  1.4442) return  7;
		if( 1.566 <= aeta and aeta <  2.0   ) return 13;
		if( 2.0   <= aeta and aeta <  2.4   ) return 19;
	}
	else if(pt >= 20. and pt < 30.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  2;
		if( 0.8   <= aeta and aeta <  1.4442) return  8;
		if( 1.566 <= aeta and aeta <  2.0   ) return 14;
		if( 2.0   <= aeta and aeta <  2.4   ) return 20;
	}
	else if(pt >= 30. and pt < 40.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  3;
		if( 0.8   <= aeta and aeta <  1.4442) return  9;
		if( 1.566 <= aeta and aeta <  2.0   ) return 15;
		if( 2.0   <= aeta and aeta <  2.4   ) return 21;
	}
	else if(pt >= 40. and pt < 50.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  4;
		if( 0.8   <= aeta and aeta <  1.4442) return 10;
		if( 1.566 <= aeta and aeta <  2.0   ) return 16;
		if( 2.0   <= aeta and aeta <  2.4   ) return 22;
	}
	else if(pt >= 50.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  5;
		if( 0.8   <= aeta and aeta <  1.4442) return 11;
		if( 1.566 <= aeta and aeta <  2.0   ) return 17;
		if( 2.0   <= aeta and aeta <  2.4   ) return 23;
	}
	else{
		cout << " ======== ERROR ================= " << endl;
		cout << " something went wrong with the index return for your muon... " << endl;
		cout << " exiting... " << endl;
		exit(-1);
	}
	return -1;
}
TString TnP::getPtEtaFromIndexEl(int i){
	TString str;
	switch (i) {
		case  0: str = "10 < pT < 15 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  6: str = "10 < pT < 15 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 12: str = "10 < pT < 15 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 18: str = "10 < pT < 15 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  1: str = "15 < pT < 20 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  7: str = "15 < pT < 20 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 13: str = "15 < pT < 20 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 19: str = "15 < pT < 20 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  2: str = "20 < pT < 30 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  8: str = "20 < pT < 30 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 14: str = "20 < pT < 30 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 20: str = "20 < pT < 30 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  3: str = "30 < pT < 40 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  9: str = "30 < pT < 40 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 15: str = "30 < pT < 40 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 21: str = "30 < pT < 40 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  4: str = "40 < pT < 50 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case 10: str = "40 < pT < 50 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 16: str = "40 < pT < 50 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 22: str = "40 < pT < 50 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  5: str = "     pT > 50 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case 11: str = "     pT > 50 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 17: str = "     pT > 50 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 23: str = "     pT > 50 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case 24: str = "               inclusive                      "; break;
	}
	return str;
}
void TnP::printElTable(){
  TFile *oFile = TFile::Open(fOutputFileName, "UPDATE");
  oFile->cd();
  
  float xbins[7]={10,15,20,30,40,50,80};
  
  string dataString;
  if(fIsData) dataString="Data"; else dataString="Sim";
  
  string flvString;
  if(fIsMu) flvString="Mu"; else flvString="El";
  

  string nameID  = flvString+dataString+"EffID";
  string nameIP  = flvString+dataString+"EffIP";
  string nameISO = flvString+dataString+"EffISO";
  string nameALL = flvString+dataString+"EffAll";
  
  TH1F* hEta1ID  = new TH1F(("eta1"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta1IP  = new TH1F(("eta1"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta1ISO = new TH1F(("eta1"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta1ALL = new TH1F(("eta1"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

  TH1F* hEta2ID  = new TH1F(("eta2"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta2IP  = new TH1F(("eta2"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta2ISO = new TH1F(("eta2"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta2ALL = new TH1F(("eta2"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

  TH1F* hEta3ID  = new TH1F(("eta3"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta3IP  = new TH1F(("eta3"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta3ISO = new TH1F(("eta3"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta3ALL = new TH1F(("eta3"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

  TH1F* hEta4ID  = new TH1F(("eta4"+nameID).c_str(), nameID.c_str() ,6,xbins);
  TH1F* hEta4IP  = new TH1F(("eta4"+nameIP).c_str(), nameIP.c_str() ,6,xbins);
  TH1F* hEta4ISO = new TH1F(("eta4"+nameISO).c_str(), nameISO.c_str() ,6,xbins);
  TH1F* hEta4ALL = new TH1F(("eta4"+nameALL).c_str(), nameALL.c_str() ,6,xbins);

    
  for(int i=0; i<fnBins; ++i){
    if( i  <= 5){ 
      hEta1ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta1ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta1IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta1IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta1ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta1ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta1ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta1ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    


    if( i > 5 && i <= 11){ 
      hEta2ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta2ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta2IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta2IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta2ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta2ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta2ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta2ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    

    if( i > 11 && i <= 17){ 
      hEta3ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta3ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta3IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta3IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta3ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta3ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta3ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta3ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    

    if( i > 17 && i <= 23){ 
      hEta4ID->SetBinContent(i%6 +1,effs[i].idEff);
      hEta4ID->SetBinError(i%6 +1,effs[i].idEffErr);
      
      hEta4IP->SetBinContent(i%6 +1,effs[i].ipEff);
      hEta4IP->SetBinError(i%6 +1,effs[i].ipEffErr);
    
      hEta4ISO->SetBinContent(i%6 +1,effs[i].isoEff);
      hEta4ISO->SetBinError(i%6 +1,effs[i].isoEffErr);
    
      hEta4ALL->SetBinContent(i%6 +1,effs[i].idEff*effs[i].ipEff*effs[i].isoEff);
      hEta4ALL->SetBinError(i%6 +1,effs[i].idEffErr+effs[i].ipEffErr+effs[i].isoEffErr);
    }    


  };
  hEta1ID->SetMaximum(1.05); hEta1ID->SetMinimum(0.20); hEta1ID->Write(hEta1ID->GetName(),TObject::kOverwrite);
  hEta1IP->SetMaximum(1.05); hEta1IP->SetMinimum(0.20); hEta1IP->Write(hEta1IP->GetName(),TObject::kOverwrite);
  hEta1ISO->SetMaximum(1.05); hEta1ISO->SetMinimum(0.20); hEta1ISO->Write(hEta1ISO->GetName(),TObject::kOverwrite);
  hEta1ALL->SetMaximum(1.05); hEta1ALL->SetMinimum(0.20); hEta1ALL->Write(hEta1ALL->GetName(),TObject::kOverwrite);

  hEta2ID->SetMaximum(1.05); hEta2ID->SetMinimum(0.20); hEta2ID->Write(hEta2ID->GetName(),TObject::kOverwrite);
  hEta2IP->SetMaximum(1.05); hEta2IP->SetMinimum(0.20); hEta2IP->Write(hEta2IP->GetName(),TObject::kOverwrite);
  hEta2ISO->SetMaximum(1.05); hEta2ISO->SetMinimum(0.20); hEta2ISO->Write(hEta2ISO->GetName(),TObject::kOverwrite);
  hEta2ALL->SetMaximum(1.05); hEta2ALL->SetMinimum(0.20); hEta2ALL->Write(hEta2ALL->GetName(),TObject::kOverwrite);

  hEta3ID->SetMaximum(1.05); hEta3ID->SetMinimum(0.20); hEta3ID->Write(hEta3ID->GetName(),TObject::kOverwrite);
  hEta3IP->SetMaximum(1.05); hEta3IP->SetMinimum(0.20); hEta3IP->Write(hEta3IP->GetName(),TObject::kOverwrite);
  hEta3ISO->SetMaximum(1.05); hEta3ISO->SetMinimum(0.20); hEta3ISO->Write(hEta3ISO->GetName(),TObject::kOverwrite);
  hEta3ALL->SetMaximum(1.05); hEta3ALL->SetMinimum(0.20); hEta3ALL->Write(hEta3ALL->GetName(),TObject::kOverwrite);

  hEta4ID->SetMaximum(1.05); hEta4ID->SetMinimum(0.20); hEta4ID->Write(hEta4ID->GetName(),TObject::kOverwrite);
  hEta4IP->SetMaximum(1.05); hEta4IP->SetMinimum(0.20); hEta4IP->Write(hEta4IP->GetName(),TObject::kOverwrite);
  hEta4ISO->SetMaximum(1.05); hEta4ISO->SetMinimum(0.20); hEta4ISO->Write(hEta4ISO->GetName(),TObject::kOverwrite);
  hEta4ALL->SetMaximum(1.05); hEta4ALL->SetMinimum(0.20); hEta4ALL->Write(hEta4ALL->GetName(),TObject::kOverwrite);

  oFile->Close();
  
  printText();
}

// ------------------------------------------
// --  STYLE STUFF  -------------------------
// ------------------------------------------

void TnP::setHistoStyle(TH1F* histo){
	histo->SetMarkerStyle(20);
	histo->SetMarkerSize(0.9);
	histo->SetLineWidth(0.8);
	gStyle->SetOptStat(1111111);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
}

// ------------------------------------------
// --  PU WEIGHT STUFF  ---------------------
// ------------------------------------------


