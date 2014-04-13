#! /usr/bin/python
import os, sys
import plotter
import ratios
import selection
import helper


class optcuts :

	def __init__(self, path, cardfile, cutspath) :
		self.path = path + '/'
		self.cardfile = cardfile
		self.cutspath = cutspath + '/'


	def optimize(self, channel) :
		for eff in range(35, 75, 5) :
			self.do_analysis(eff)


	def do_analysis(self, eff) :
		cuts = opt.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (eff))

		pl = plotter.plotter(self.path, self.cardfile)

		# selections
		sels = {}
		sels['Normal'] = self.set_selection('eff%d' % eff, cuts, pl.systematics['Normal'])

		# charge mis-ID scale factor
		pl.chmid_sf = 1.62

		# get fake and prompt ratios
		EWK_SF = {}
		EWK_SF['el']   = pl.get_EWK_SF('el')
		EWK_SF['mu17'] = pl.get_EWK_SF('mu17')
		EWK_SF['mu24'] = pl.get_EWK_SF('mu24')
		pl.fpr.fill_ratios(pl.get_samples('SingleDoubleMu'), pl.get_samples('DoubleEle'), 0, True, EWK_SF)

		results = {}
		resultspath = self.path + 'optcuts/results_eff%d.pkl' % eff

		if os.path.exists(resultspath) :
			print '[status] loading results of predictions from %s..' % (resultspath)
			results = helper.load_object(resultspath)

		else :
			results['Normal'] = pl.make_IntPredictions(sels['Normal'], self.path + 'SSDLYields_skim_Normal.root')
			helper.save_object(results, resultspath)

		# make datacards for each charge-flavor channel
		for ch_str in results['Normal'] :
			for chan in results['Normal'][ch_str] :
				if ch_str != '++' : continue
				datacard = pl.make_datacard(results, chan, ch_str, 'optcuts', 'eff%d' % eff)

#				with open(datacard, 'r+') as file :
#					if not any(['LepUp' in line for line in file.readlines()]) :
#						file.write('adding LepUp')  # TODO: add default line here

#		pl.make_datacard(results, 'al', '++', 'optcuts', 'eff%d' % (eff))
#		pl.make_datacard(results, 'al', 'al', 'optcuts', 'eff%d' % (eff))


	def read_cuts(self, cutfile) :
#		'''reads file with cuts'''
#		if not os.path.exists(cutfile) :
#			print '[ERROR] %s does not exist!' % (cutfile)
#			sys.exit(1)
#		print '[status] reading cuts from %s' % (cutfile)
#HT 378.54 1e+30
#pT1 61.7743 1e+30
#pT2 61.7743 1e+30
#NJ 3 100000.
#NbJ 1 100000.
#NbJmed 1 100000.
#Charge 1 10
		cuts = {}
		with open(cutfile, 'r') as file :
			for i, line in enumerate(file.readlines()) :
				if line[0] is '#' : continue
				splitline = line.split()
				if (len(splitline) != 3) :
					print '[ERROR] check format of %s (line %d)!' % (cutfile, i+1)
					sys.exit(1)
				[var, min, max] = splitline
				min = float(min)
				max = float(max)
				cuts[var] = [min, max]
				print splitline
#				[name, inputfile, datamc, channel] = splitline[:4]
#				datamc = int(datamc)
#				channel = int(channel)
#				xsec = -1.
#				if len(splitline) == 5 :
#					xsec = float(splitline[4])
#				samples[name] = sample.sample(name = name, datamc = datamc, channel = channel, xsec = xsec, ngen = -1)
#				if verbose > 0 : print samples[name]
#		return samples
		return cuts


	def set_selection(self, name, cuts, systflag = 0) :
		'''sets selcetion'''

		sel = selection.selection(name)
		sel.systflag = systflag

		for var, cut in cuts.iteritems() :
			if var == 'HT' :
				sel.minHT = float(cut[0])
				sel.maxHT = float(cut[1])
			if var == 'pT1' :
				sel.minPt1 = float(cut[0])
				sel.maxPt1 = float(cut[1])
			if var == 'pT2' :
				sel.minPt2 = float(cut[0])
				sel.maxPt2 = float(cut[1])
			if var == 'NJ' :
				sel.minNjets = int(cut[0])
				sel.maxNjets = int(cut[1])
			if var == 'NbJmed' :
				sel.minNbjetsM = int(cut[0])
				sel.maxNbjetsM = int(cut[1])
#			if var == 'Charge' :
#				sel.charge = int(cut[0])

		return sel

#NbJmed 1 100000.
#Charge 1 10

if __name__ == '__main__' :
	args = sys.argv

	if ('--help' in args) or ('-h' in args) or ('-d' not in args) or ('-c' not in args) or ('--cuts' not in args) or ('--channel' not in args) :
		print 'usage: ./optcuts.py -d <INPUTDIR> -c <DATACARD> -s <CUTSDIR>'
		print ''
		print '       -d Directory in which the SSDLYYields.root output from the SDLDumper is.'
		print ''
		print '       -c Datacard with list of samples and corresponding cross sections.'
		print ''
		print '       -cutsDir Directory in which the output files of TMVA wich selection cuts are.'
		sys.exit(1)

	if ('-d' in args) and (args[args.index('-d')+1] != '') :
		path = str(args[args.index('-d')+1])
		print path

	if ('-c' in args) and (args[args.index('-c')+1] != '') :
		cardfile = str(args[args.index('-c')+1])
		print cardfile

	if ('--cuts' in args) and (args[args.index('--cuts')+1] != '') :
		cutspath = str(args[args.index('--cuts')+1])

	if ('--channel' in args) and (args[args.index('--channel')+1] != '') :
		channel = str(args[args.index('--channel')+1])

	opt = optcuts(path, cardfile, cutspath)
	opt.optimize(channel)
#	plotter = plotter(path)
#	plotter.do_analysis(cardfile)
