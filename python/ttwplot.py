#! /usr/bin/python
import ROOT
import helper


class ttwplot :

	def __init__(self) :
		foo = 0


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
		hs_pred        .SetMaximum(maximum)

		for process, histo in histos.iteritems() :
			histo.SetMaximum(maximum)
			histo.SetLineColor(1)
			histo.SetLineWidth(1)

		# set colors
		histos['fake' ].SetFillColor(46)
		histos['chmid'].SetFillColor(49)
		histos['rare' ].SetFillColor(38)
		histos['wz'   ].SetFillColor(39)
		histos['ttz'  ].SetFillColor(42)
		histos['ttw'  ].SetFillColor(44)

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
		leg.AddEntry(histos['obs'  ], 'Observed'         , 'p')
		leg.AddEntry(histos['fake' ], 'Non-prompt lepton', 'f')
		leg.AddEntry(histos['chmid'], 'Charge MisID'     , 'f')
		leg.AddEntry(histos['rare' ], 'Rare SM'          , 'f')
		leg.AddEntry(histos['wz'   ], 'WZ'               , 'f')
		leg.AddEntry(histos['ttz'  ], 't#bar{t} + Z'     , 'f')
		leg.AddEntry(histos['ttw'  ], 't#bar{t} + W'     , 'f')
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
