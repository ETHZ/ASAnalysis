#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>


float getSLumi(string* sname){
	if      (*sname == "HTauTau" ) return 5.48751e+06;
	else if (*sname == "HWW"     ) return 769036.;
	else if (*sname == "HZZ"     ) return 1.56315e+07;
	else if (*sname == "TTbarG"  ) return 778052.;
	else if (*sname == "TTbarWW" ) return 1.06784e+08;
	else if (*sname == "TTbarW"  ) return 821246.;
	else if (*sname == "TTbarZ"  ) return 996427.;
	else if (*sname == "TTbarZ"  ) return 996427.;
	else if (*sname == "TbZ"     ) return 1.30003e+07;
	else if (*sname == "W+W+"    ) return 400290.;
	else if (*sname == "W-W-"    ) return 1.08112e+06;
	else if (*sname == "WWG"     ) return 406856.;
	else if (*sname == "WWW"     ) return 2.70596e+06;
	else if (*sname == "WWZ"     ) return 3.7532e+06;
	else if (*sname == "WZTo3LNu") return 1.90712e+06;
	else if (*sname == "WZZ"     ) return 1.11552e+07;
	else if (*sname == "ZZTo4L"  ) return 2.71718e+07;
	else if (*sname == "ZZZ"     ) return 4.06372e+07;
	else if (*sname == "TTJets"  ) return 87666.1;
	else if (*sname == "TTJets_madgraph_v1"  ) return 29465.8;
	else if (*sname == "TTJets_madgraph_v2"  ) return 5831.55;
	else if (*sname == "TTJets_v1"     ) return 26338.3;
	else if (*sname == "SingleT_tW"    ) return 44416.4;
	else if (*sname == "SingleT_s"     ) return 68273.6;
	else if (*sname == "SingleT_t"     ) return 1776.88;
	else if (*sname == "SingleTbar_tW" ) return 44040.8;
	else if (*sname == "SingleTbar_s"  ) return 78945.4;
	else if (*sname == "SingleTbar_t"  ) return 64453.3;
	else if (*sname == "WJets"         ) return 480.526;
	else if (*sname == "DYJets"        ) return 8470.7 ;
	else if (*sname == "GJets200"      ) return 57068.5;
	else if (*sname == "GJets400"      ) return 371646.;
	else if (*sname == "WWTo2L2Nu"     ) return 332404.;
	else if (*sname == "SingleTbar_t"  ) return 64453.3;
	else if (*sname == "DPSWW"         ) return 1.41816e+06;
	else if (*sname == "DPSWW"         ) return 1.41816e+06;
	else if (*sname == "WGstarMu"      ) return 156712.; // ARE THOSE REALLY 8 TEV SAMPLES?
	else if (*sname == "WGstarTau"     ) return 148723.; // ARE THOSE REALLY 8 TEV SAMPLES?
	return 0.;  // this will protect against having forgotten a sample by nan

}

template <class T> inline void getObjectSafe(TFile* pFile, TString name, T*& object){
    pFile->GetObject(name, object);
    if(!object){
        std::cout << name + " not found!" << std::endl;
        exit(-1);
    }
    return;
};


int loopTree(){

	gStyle->SetStatX(0.8);
	gStyle->SetStatY(0.5);

	float intlumi = 19500.;

	TString inputFile = "SSDLYields.root";
	TFile *pFile = TFile::Open(inputFile);
	TTree *sigtree; getObjectSafe(pFile, "SigEvents", sigtree);
	
	string *sname = 0;
	int flag;
	int   SType, Flavor, TLCat, NJ, NbJ, NbJmed;
	float puweight, pT1, pT2, HT, MET, SLumi, HLTSF;
	float mll, jet0ptb, jet1ptb;
	int   passZVeto;
	
	sigtree->SetBranchAddress("SystFlag", &flag);
	// sigtree->SetBranchAddress("Run",      &run);
	sigtree->SetBranchAddress("SName",    &sname);
	sigtree->SetBranchAddress("SType",    &SType);
	sigtree->SetBranchAddress("TLCat",    &TLCat);
	sigtree->SetBranchAddress("HLTSF",    &HLTSF);
	sigtree->SetBranchAddress("PUWeight", &puweight);
	sigtree->SetBranchAddress("SLumi",    &SLumi);
	sigtree->SetBranchAddress("Flavor",   &Flavor);
	sigtree->SetBranchAddress("pT1",      &pT1);
	sigtree->SetBranchAddress("pT2",      &pT2);
	sigtree->SetBranchAddress("HT",       &HT);
	sigtree->SetBranchAddress("MET",      &MET);
	sigtree->SetBranchAddress("NJ",       &NJ);
	sigtree->SetBranchAddress("NbJ",      &NbJ);
	sigtree->SetBranchAddress("NbJmed",   &NbJmed);
	sigtree->SetBranchAddress("PassZVeto",&passZVeto);
	sigtree->SetBranchAddress("Mll",&mll);
	sigtree->SetBranchAddress("Jet0PtB",&jet0ptb);
	sigtree->SetBranchAddress("Jet1PtB",&jet1ptb);

	float minPt1( 20.), minPt2(20.);
	float minHT ( 80.), maxHT (8000.);
	float minMET( 40.), maxMET(8000.);
	float minNJ (  2 );
	float minNB (  2 );
	bool  ss( false );


	// invariant mass of the leptons
	TH1D * MLL_tt_dat = new TH1D("mll_tt_dat", "mll_tt_dat", 25, 0., 250.); MLL_tt_dat->Sumw2(); MLL_tt_dat->SetMarkerStyle(20); MLL_tt_dat->SetMarkerSize(1.2);
	TH1D * MLL_tt_ttj = new TH1D("mll_tt_ttj", "mll_tt_ttj", 25, 0., 250.); MLL_tt_ttj->Sumw2(); MLL_tt_ttj->SetFillColor(46); 
	TH1D * MLL_tt_rar = new TH1D("mll_tt_rar", "mll_tt_rar", 25, 0., 250.); MLL_tt_rar->Sumw2(); MLL_tt_rar->SetFillColor(38);
	TH1D * MLL_tt_ttw = new TH1D("mll_tt_ttw", "mll_tt_ttw", 25, 0., 250.); MLL_tt_ttw->Sumw2(); MLL_tt_ttw->SetFillColor(44);

	// pt of all b jets
	TH1D * BPT_tt_dat = new TH1D("bpt_tt_dat", "bpt_tt_dat", 25, 0., 250.); BPT_tt_dat->Sumw2(); BPT_tt_dat->SetMarkerStyle(20); BPT_tt_dat->SetMarkerSize(1.2);
	TH1D * BPT_tt_ttj = new TH1D("bpt_tt_ttj", "bpt_tt_ttj", 25, 0., 250.); BPT_tt_ttj->Sumw2(); BPT_tt_ttj->SetFillColor(46); 
	TH1D * BPT_tt_rar = new TH1D("bpt_tt_rar", "bpt_tt_rar", 25, 0., 250.); BPT_tt_rar->Sumw2(); BPT_tt_rar->SetFillColor(38);
	TH1D * BPT_tt_ttw = new TH1D("bpt_tt_ttw", "bpt_tt_ttw", 25, 0., 250.); BPT_tt_ttw->Sumw2(); BPT_tt_ttw->SetFillColor(44);

	// pt of the first bjet
	TH1D * B1PT_tt_dat = new TH1D("b1pt_tt_dat", "b1pt_tt_dat", 25, 0., 250.); B1PT_tt_dat->Sumw2(); B1PT_tt_dat->SetMarkerStyle(20); B1PT_tt_dat->SetMarkerSize(1.2);
	TH1D * B1PT_tt_ttj = new TH1D("b1pt_tt_ttj", "b1pt_tt_ttj", 25, 0., 250.); B1PT_tt_ttj->Sumw2(); B1PT_tt_ttj->SetFillColor(46); 
	TH1D * B1PT_tt_rar = new TH1D("b1pt_tt_rar", "b1pt_tt_rar", 25, 0., 250.); B1PT_tt_rar->Sumw2(); B1PT_tt_rar->SetFillColor(38);
	TH1D * B1PT_tt_ttw = new TH1D("b1pt_tt_ttw", "b1pt_tt_ttw", 25, 0., 250.); B1PT_tt_ttw->Sumw2(); B1PT_tt_ttw->SetFillColor(44);

	// pt of the second bjet
	TH1D * B2PT_tt_dat = new TH1D("b2pt_tt_dat", "b2pt_tt_dat", 25, 0., 250.); B2PT_tt_dat->Sumw2(); B2PT_tt_dat->SetMarkerStyle(20); B2PT_tt_dat->SetMarkerSize(1.2);
	TH1D * B2PT_tt_ttj = new TH1D("b2pt_tt_ttj", "b2pt_tt_ttj", 25, 0., 250.); B2PT_tt_ttj->Sumw2(); B2PT_tt_ttj->SetFillColor(46); 
	TH1D * B2PT_tt_rar = new TH1D("b2pt_tt_rar", "b2pt_tt_rar", 25, 0., 250.); B2PT_tt_rar->Sumw2(); B2PT_tt_rar->SetFillColor(38);
	TH1D * B2PT_tt_ttw = new TH1D("b2pt_tt_ttw", "b2pt_tt_ttw", 25, 0., 250.); B2PT_tt_ttw->Sumw2(); B2PT_tt_ttw->SetFillColor(44);


	int nttw(0), nttj(0), nrare(0);

	// looping on the tree
	for( int i = 0; i < sigtree->GetEntries(); i++ ){
		sigtree->GetEntry(i);

		if (*sname == "TTbarWNLO" ) continue; // avoid double counting here...

		if ( !passZVeto)  continue; // z-veto
		if ( flag != 0 )  continue; // no systematics
		if ( mll    < 8.) continue; // trigger cut
		// selection cuts
		if ( HT     < minHT  || HT  > maxHT ) continue;
		if ( MET    < minMET || MET > maxMET) continue;
		if ( NJ     < minNJ)    continue;
		if ( NbJmed < minNB)    continue;

		// cut on the pts. taking emu into account
		if(Flavor == 1 || Flavor == 4){
			if(pT1 > pT2){
				if(pT1 < minPt1) continue;
				if(pT2 < minPt2) continue;
			}
			if(pT1 < pT2){
				if(pT1 < minPt2) continue;
				if(pT2 < minPt1) continue;
			}
		}
		else{
			if(pT1 < minPt1) continue;
			if(pT2 < minPt2) continue;
		}
		//------------------

		if (TLCat != 0) continue; // only tight tight for now


		if ( ss && Flavor  > 2) continue; // only same-sign
		if (!ss && Flavor != 4) continue; // only e-mu for OS

		float scale(1.); if (SType  > 3) scale = puweight*HLTSF*intlumi / getSLumi(sname); // event weight for mc events

	// cout << *sname << Form(":   puweight: %.2f   HLTSF: %.2f   slumi: %.2f  ", puweight, HLTSF, SLumi) << endl;

		if (SType < 3){								// data
			MLL_tt_dat  -> Fill(mll);
			BPT_tt_dat  -> Fill(jet0ptb); BPT_tt_dat  -> Fill(jet1ptb);
			B1PT_tt_dat -> Fill(jet0ptb);
			B2PT_tt_dat -> Fill(jet1ptb);
		}
		else if (SType == 15 && *sname != "TTbarW"){ // rares without ttw
			MLL_tt_rar  -> Fill(mll, scale);
			BPT_tt_rar  -> Fill(jet0ptb, scale); BPT_tt_rar  -> Fill(jet1ptb, scale);
			B1PT_tt_rar -> Fill(jet0ptb, scale);
			B2PT_tt_rar -> Fill(jet1ptb, scale);
			nrare++;
		}
		else if (*sname == "TTJets"             ||
		         *sname == "TTJets_v1"          ||
		         *sname == "TTJets_madgraph_v1" ||
		         *sname == "TTJets_madgraph_v2" ||
		         *sname == "SingleT_t"          ||
		         *sname == "SingleT_tW"         ||
		         *sname == "SingleT_s"          ||
		         *sname == "SingleTbar_t"       ||
		         *sname == "SingleTbar_tW"      ||
		         *sname == "SingleTbar_s"       ){  // top samples
			MLL_tt_ttj  -> Fill(mll, scale);
			BPT_tt_ttj  -> Fill(jet0ptb, scale); BPT_tt_ttj  -> Fill(jet1ptb, scale);
			B1PT_tt_ttj -> Fill(jet0ptb, scale);
			B2PT_tt_ttj -> Fill(jet1ptb, scale);
			nttj++;
		}
		else if (*sname == "TTbarW"){ 				// ttw
			MLL_tt_ttw  -> Fill(mll, scale);
			BPT_tt_ttw  -> Fill(jet0ptb, scale); BPT_tt_ttw  -> Fill(jet1ptb, scale);
			B1PT_tt_ttw -> Fill(jet0ptb, scale);
			B2PT_tt_ttw -> Fill(jet1ptb, scale);
			nttw++;
		}
			
	}
	
	cout << Form("nttj: %d   nttw: %d   nrare: %d",nttj, nttw, nrare) << endl;

	TLegend *leg = new TLegend(0.70,0.63,0.90,0.90);
	leg->SetFillColor(0);
	leg->AddEntry(MLL_tt_dat, "Data","p");
	leg->AddEntry(MLL_tt_ttj, "Top","f");
	leg->AddEntry(MLL_tt_rar, "Rares","f");
	leg->AddEntry(MLL_tt_ttw, "t#bar{t} + W","f");

	cout << Form("integrals: \n\t ttj: %.2f \n\t ttw: %.2f \n\t rar: %.2f \n\t data: %d", MLL_tt_ttj->Integral(), MLL_tt_ttw->Integral(), MLL_tt_rar->Integral(), (int) MLL_tt_dat->Integral()) << endl;

	TCanvas * c = new TCanvas("foo", "bar");

	// stacking and plotting stuff
	//----------------------------
	// mll
	THStack * MLL_tt_stack = new THStack("mll_stack", "m_{ll}");
	MLL_tt_stack->Add(MLL_tt_ttj);
	MLL_tt_stack->Add(MLL_tt_rar);
	MLL_tt_stack->Add(MLL_tt_ttw);
	MLL_tt_stack->Draw("hist");
	MLL_tt_dat  ->Draw("P same");
	MLL_tt_stack->GetXaxis()->SetTitle("m_{ll}");
	MLL_tt_stack->SetMaximum(1.2*TMath::Max(MLL_tt_stack->GetMaximum(), MLL_tt_dat->GetMaximum() ) );
	leg->Draw("same");
	c->SaveAs("mll.pdf");
	
	// pt of all b jets
	THStack * BPT_tt_stack = new THStack("bpt_stack", "p_{T}^{b}");
	BPT_tt_ttj->Scale(BPT_tt_dat->Integral()/BPT_tt_ttj->Integral());
	BPT_tt_stack->Add(BPT_tt_ttj);
	BPT_tt_stack->Add(BPT_tt_rar);
	BPT_tt_stack->Add(BPT_tt_ttw);
	BPT_tt_stack->Draw("hist");
	BPT_tt_dat  ->Draw("P sames");
	BPT_tt_stack->GetXaxis()->SetTitle("p_{T}^{b}");
	BPT_tt_stack->SetMaximum(1.2*TMath::Max(BPT_tt_stack->GetMaximum(), BPT_tt_dat->GetMaximum() ) );
	leg->Draw("same");
	c->SaveAs("bpt.pdf");
	
	// pt of the first b jet
	THStack * B1PT_tt_stack = new THStack("b1pt_stack", "p_{T}^{b1}");
	B1PT_tt_ttj->Scale(B1PT_tt_dat->Integral()/B1PT_tt_ttj->Integral());
	B1PT_tt_stack->Add(B1PT_tt_ttj);
	B1PT_tt_stack->Add(B1PT_tt_rar);
	B1PT_tt_stack->Add(B1PT_tt_ttw);
	B1PT_tt_stack->Draw("hist");
	B1PT_tt_dat  ->Draw("P sames");
	B1PT_tt_stack->GetXaxis()->SetTitle("p_{T}^{b1}");
	B1PT_tt_stack->SetMaximum(1.2*TMath::Max(B1PT_tt_stack->GetMaximum(), B1PT_tt_dat->GetMaximum() ) );
	leg->Draw("same");
	c->SaveAs("b1pt.pdf");

	// pt of the second b jet
	THStack * B2PT_tt_stack = new THStack("b2pt_stack", "p_{T}^{b2}");
	B2PT_tt_ttj->Scale(B2PT_tt_dat->Integral()/B2PT_tt_ttj->Integral());
	B2PT_tt_stack->Add(B2PT_tt_ttj);
	B2PT_tt_stack->Add(B2PT_tt_rar);
	B2PT_tt_stack->Add(B2PT_tt_ttw);
	B2PT_tt_stack->Draw("hist");
	B2PT_tt_dat  ->Draw("P sames");
	B2PT_tt_stack->GetXaxis()->SetTitle("p_{T}^{b2}");
	B2PT_tt_stack->SetMaximum(1.2*TMath::Max(B2PT_tt_stack->GetMaximum(), B2PT_tt_dat->GetMaximum() ) );
	leg->Draw("same");
	c->SaveAs("b2pt.pdf");


	return 1;
}
