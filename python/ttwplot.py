#! /usr/bin/python
import ROOT
import helper


class ttwplot :

	def __init__(self, path) :
		self.path = path
		helper.mkdir(self.path)

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


	def get_fillColor(self, process) :
		if process in self.colors.keys() :
			return self.colors[process]
		else :
			return 0


	def save_plot(self, histos, path, var_name) :
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
		maximum = 1.2 * max(histos['obs'].GetMaximum(), histos['pred'].GetMaximum())
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

		histos['pred' ].SetLineWidth(3)
		histos['pred' ].SetFillColor(12)
		histos['pred' ].SetFillStyle(3005)

		# legend
		leg = ROOT.TLegend()
		leg.AddEntry(histos['obs'  ], self.process_names['obs'  ], 'p')
		leg.AddEntry(histos['fake' ], self.process_names['fake' ], 'f')
		leg.AddEntry(histos['chmid'], self.process_names['chmid'], 'f')
		leg.AddEntry(histos['rare' ], self.process_names['rare' ], 'f')
		leg.AddEntry(histos['wz'   ], self.process_names['wz'   ], 'f')
		leg.AddEntry(histos['ttz'  ], self.process_names['ttz'  ], 'f')
		leg.AddEntry(histos['ttw'  ], self.process_names['ttw'  ], 'f')
		# set position
		width = 0.17
		x = 0.68
		y = 0.88
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
		canvas.SetLeftMargin(0.12)
		canvas.SetRightMargin(0.04)
		canvas.cd()

		hs_pred.Draw()
		leg.Draw()
		histos['pred'].Draw('0 E2 same')
		histos['bgtot'].Draw('same')
#		gr_obs.Draw('P same')
		histos['obs'  ].Draw('PE X0 same')
#		histos['fake'].Draw('same')
		raw_input('ok? ')

		canvas.Print(path + 'HT' + '.pdf')
