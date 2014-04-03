#! /usr/bin/python
import ROOT
import helper


class ttwplot :

	def __init__(self, path, lumi = 19500.) :
		self.path = path
		helper.mkdir(self.path)
		self.lumi = lumi

		# colors
		self.colors = {}
		self.colors['obs'  ] = 1
		self.colors['fake' ] = 46
		self.colors['chmid'] = 49
		self.colors['rare' ] = 38
		self.colors['wz'   ] = 39
		self.colors['ttz'  ] = 42
		self.colors['ttw'  ] = 44

		# process names
		self.process_names = {}
		self.process_names['obs'  ] = 'Observed'
		self.process_names['fake' ] = 'Non-prompt lepton'
		self.process_names['chmid'] = 'Charge MisID'
		self.process_names['rare' ] = 'Rare SM'
		self.process_names['wz'   ] = 'WZ'
		self.process_names['ttz'  ] = 't#bar{t} + Z'
		self.process_names['ttw'  ] = 't#bar{t} + W'
		self.process_names['bgtot'] = 'Backgrounds'

		# variable names
		self.var_names = {}
		self.var_names['HT']  = 'H_{T} [GeV]'
		self.var_names['Mll'] = 'm_{ll} [GeV]'


	def get_fillColor(self, process) :
		if process in self.colors.keys() :
			return self.colors[process]
		return 0


	def get_varName(self, var) :
		if var in self.var_names.keys() :
			return self.var_names[var]
		return ''


	def save_plot(self, histos, path, var) :
		'''save plot with observation and predictions'''

		gr_obs  = helper.getGraphPoissonErrors(histos['obs'])
		hs_pred = ROOT.THStack('hs_pred', 'hs_pred')

		# adding predictions to stack
		hs_pred.Add(histos['fake' ])
		hs_pred.Add(histos['chmid'])
		hs_pred.Add(histos['rare' ])
		hs_pred.Add(histos['wz'   ])
		hs_pred.Add(histos['ttz'  ])
		hs_pred.Add(histos['ttw'  ])

		# set minimum and maximum
		maximum = 1.8 * max(histos['obs'].GetMaximum(), histos['pred'].GetMaximum())
		hs_pred.SetMaximum(maximum)

		for process, histo in histos.iteritems() :
			histo.SetMaximum(maximum)
			histo.SetLineColor(1)
			histo.SetLineWidth(1)
			histo.SetFillColor(self.get_fillColor(process))

		# data histogram settings
		histos['obs'  ].SetMarkerStyle(20)
		histos['obs'  ].SetMarkerSize(1.1)
		histos['obs'  ].SetLineWidth(2)
		gr_obs.SetMarkerStyle(20)
		gr_obs.SetMarkerSize(1.1)
		gr_obs.SetLineWidth(2)
		
		# special settings for total BG line and uncertainty
		histos['bgtot'].SetLineWidth(3)
		histos['bgtot'].SetLineColor(1)
#		histos['bgtot'].SetFillColor(12)
		histos['bgtot'].SetFillStyle(0)

		histos['pred' ].SetLineWidth(0)
		histos['pred' ].SetFillColor(12)
		histos['pred' ].SetFillStyle(3005)

		# legend
		leg = ROOT.TLegend()
		leg.AddEntry(histos['obs'  ], self.process_names['obs'  ], 'lp')
		leg.AddEntry(histos['fake' ], self.process_names['fake' ], 'f')
		leg.AddEntry(histos['chmid'], self.process_names['chmid'], 'f')
		leg.AddEntry(histos['rare' ], self.process_names['rare' ], 'f')
		leg.AddEntry(histos['wz'   ], self.process_names['wz'   ], 'f')
		leg.AddEntry(histos['ttz'  ], self.process_names['ttz'  ], 'f')
		leg.AddEntry(histos['ttw'  ], self.process_names['ttw'  ], 'f')
		leg.AddEntry(histos['bgtot'], self.process_names['bgtot'], 'l')
		leg.AddEntry(histos['pred' ], 'BG uncertainty', 'fl')
		# set position
		width = 0.17
		x = 0.65
		y = 0.93
		leg.SetX1NDC(x)
		leg.SetX2NDC(x+width)
		leg.SetY1NDC(y-leg.GetNRows()*0.25*width)
		leg.SetY2NDC(y)
		leg.SetFillStyle(0)
		leg.SetTextFont(42)
		leg.SetTextSize(0.03)
		leg.SetBorderSize(0)
		leg.SetTextAlign(12)

		canvas = ROOT.TCanvas('C_ObsPred', 'Observed vs Predicted', 0, 0, 600, 600)
#		canvas.SetLeftMargin(0.12)
#		canvas.SetRightMargin(0.04)
		print 'left'  , canvas.GetLeftMargin()
		print 'right' , canvas.GetRightMargin()
		print 'top'   , canvas.GetTopMargin()
		print 'bottom', canvas.GetBottomMargin()
		canvas.SetLeftMargin(0.12)
		canvas.SetRightMargin(0.04)
		canvas.SetTopMargin(0.04)
		canvas.SetBottomMargin(0.12)
		canvas.cd()

#		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptTitle(0)
		ROOT.gPad.SetTicks(1,1)

		# draw and set axis titles
		hs_pred.Draw()
		hs_pred.GetXaxis().SetTitle(self.get_varName(var))
		hs_pred.GetYaxis().SetTitle('Events')
		hs_pred.GetXaxis().SetTitleOffset(1.25)
		hs_pred.GetYaxis().SetTitleOffset(1.25)
		hs_pred.GetXaxis().SetTitleSize(0.04)
		hs_pred.GetYaxis().SetTitleSize(0.04)
		hs_pred.GetXaxis().SetLabelSize(0.04)
		hs_pred.GetYaxis().SetLabelSize(0.04)
		print hs_pred.GetXaxis().GetTitleSize()
		##		hs_pred[var]->Draw("goff");
		##		hs_pred[var]->GetXaxis()->SetTitle(xAxisTitle[var].Data());
		##		hs_pred[var]->GetXaxis()->SetTitleOffset(1.07);
		##		hs_pred[var]->GetYaxis()->SetTitle(yAxisTitle[var].Data());
		##		hs_pred[var]->GetYaxis()->SetTitleSize(0.045);
		##		hs_pred[var]->GetXaxis()->SetTitleSize(0.045);
		##		hs_pred[var]->GetYaxis()->SetLabelSize(0.045);
		##		hs_pred[var]->GetXaxis()->SetLabelSize(0.045);
		##		hs_pred[var]->GetYaxis()->SetTitleOffset(1.25);
		##		hs_pred[var]->GetXaxis()->SetTitleOffset(1.065);
		##		if (diffVarName[var] == "NJ" || diffVarName[var] == "NbJmed"){
		##			for(size_t i = 1; i <= nbins[var]; ++i) hs_pred[var]->GetXaxis()->SetBinLabel(i, Form("%d", (int)bins[var][i-1]));
		##			hs_pred[var]->GetXaxis()->SetLabelSize(0.07);
		##			hs_pred[var]->GetXaxis()->SetTitleSize(0.045);
		##			hs_pred[var]->GetXaxis()->SetTitleOffset(1.07);
		##		}
		##		if (diffVarName[var] == "Int") {
		##			for(size_t i = 1; i <= nbins[var]; ++i) {
		##				TString binlabel = "?";
		##				if (i == 1) binlabel = "e"  + chargeSignString + "e"  + chargeSignString;
		##				if (i == 2) binlabel = "#mu"+ chargeSignString + "#mu"+ chargeSignString;
		##				if (i == 3) binlabel = "e"  + chargeSignString + "#mu"+ chargeSignString;
		##				hs_pred[var]->GetXaxis()->SetBinLabel(i, binlabel);
		##			}
		##		}
		leg.Draw()
		histos['pred'].Draw('0 E2 same')
		histos['bgtot'].Draw('hist same')
#		gr_obs.Draw('P same')
		histos['obs'  ].Draw('PE X0 same')
#		self.drawTopLine(0.56, 0.8)
		self.draw_cmsLine()
		raw_input('ok? ')

		canvas.Print(path + 'HT' + '.pdf')


	def drawTopLine(self, rightedge = 0.60, scale = 1., leftedge = 0.13) :
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(scale * 0.05)
		latex.DrawLatex(leftedge, 0.92, 'CMS Preliminary')
		latex.SetTextFont(42)
		latex.SetTextSize(scale * 0.04)
		latex.DrawLatex(rightedge, 0.92, '#sqrt{s} = 8 TeV, L_{int} = %4.1f fb^{-1}' % (self.lumi/1000.))


	def draw_cmsLine(self) :
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(0.04)
		latex.SetTextAlign(13)
		latex.DrawLatex(0.15, 0.93, 'CMS')
		latex.SetTextFont(42)
		latex.SetTextSize(0.03)
		latex.DrawLatex(0.15, 0.88, '#sqrt{s} = 8 TeV, L_{int} = %4.1f fb^{-1}' % (self.lumi/1000.))
