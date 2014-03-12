#!/usr/bin/python
import ROOT

pdfsets = []
pdfsets.append(('mstw2008', 0.00008988, 0.00000184, 0.00000205))
pdfsets.append(('nnpdf20 ', 0.00009225, 0.00000227, 0.00000222))
pdfsets.append(('ct10    ', 0.00009001, 0.00000319, 0.00000255))

central = 0.0000857803883223637

# print latex table
print '\\begin{tabular}{l|rrrrr}'
print '\t\\hline\\hline'
print 'PDF set & A_0 & \\Delta A^+ & \\Delta A^- & A_0 + \\Delta A^+ & A_0 - \\Delta A^- \\\\'
print '\\hline'
print 'cteq6ll & %.8f & & & & \\\\' %(central)
print '\\hline'
for pdfset in pdfsets :
	print '%s & %.8f & %.8f & %.8f' %(pdfset), '& %.8f & %.8f \\\\' %(pdfset[1]+pdfset[2], pdfset[1]-pdfset[3])
print '\\hline\\hline'

histo = ROOT.TH1D('acceptances', 'acceptances', 3, 0., 3.)
print range(histo.GetNbinsX())
for i in range(histo.GetNbinsX()) :
	histo.SetBinContent(i+1, pdfsets[i][1])

asym = ROOT.TGraphAsymmErrors(histo)
for i in range(histo.GetNbinsX()) :
	asym.SetPointEXhigh(i+1, pdfsets[i][2])
	asym.SetPointEXlow (i+1, pdfsets[i][3])

histo.Draw()
asym.Draw('e2same')
raw_input('ok? ')
