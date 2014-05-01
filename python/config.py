#! /usr/bin/python
import selection


def get_histoBins(var, sel = '', total_bin = False) :
		histo_settings = {}
		if   var is 'HT'     : histo_settings['nbins'] =  5; histo_settings['min'] = 100.; histo_settings['max'] = 600.;
		elif var is 'MET'    : histo_settings['nbins'] =  6; histo_settings['min'] =   0.; histo_settings['max'] = 120.;
		elif var is 'NJ'     : histo_settings['nbins'] =  4; histo_settings['min'] =   2.; histo_settings['max'] =   6.;
		elif var is 'NbJmed' : histo_settings['nbins'] =  4; histo_settings['min'] =   0.; histo_settings['max'] =   4.;
		elif var is 'pT1'    : histo_settings['nbins'] =  9; histo_settings['min'] =  20.; histo_settings['max'] = 200.;
		elif var is 'pT2'    : histo_settings['nbins'] =  8; histo_settings['min'] =  20.; histo_settings['max'] = 100.;
		elif var is 'Mll'    : histo_settings['nbins'] = 14; histo_settings['min'] =  20.; histo_settings['max'] = 300.;
		elif var is 'NVrtx'  : histo_settings['nbins'] = 40; histo_settings['min'] =   0.; histo_settings['max'] =  40.;
		elif var is 'minMT'  : histo_settings['nbins'] = 20; histo_settings['min'] =   0.; histo_settings['max'] = 400.;
		elif var is 'M3'     : histo_settings['nbins'] = 12; histo_settings['min'] =   0.; histo_settings['max'] = 600.;
		elif var is 'Int'    : histo_settings['nbins'] =  3; histo_settings['min'] =   0.; histo_settings['max'] =   3.;
		elif var is 'CFChan' and not total_bin : histo_settings['nbins'] =  6; histo_settings['min'] =   0.; histo_settings['max'] =   6.;  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		elif var is 'CFChan' and     total_bin : histo_settings['nbins'] =  7; histo_settings['min'] =  -1.; histo_settings['max'] =   6.;  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		elif var is 'Charge' : histo_settings['nbins'] =  2; histo_settings['min'] =  -2.; histo_settings['max'] =   2.;  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.

		if sel != '' :
			if sel.name.startswith('1J') :
				if   var is 'NJ'     : histo_settings['nbins'] =  5; histo_settings['min'] =   1.; histo_settings['max'] =   6.;
			if sel.name.startswith('3J1bJ') or sel.name.startswith('final') :
				if   var is 'NJ'     : histo_settings['nbins'] =  3; histo_settings['min'] =   3.; histo_settings['max'] =   6.;
				elif var is 'NbJmed' : histo_settings['nbins'] =  3; histo_settings['min'] =   1.; histo_settings['max'] =   4.;
				if sel.name.startswith('final') :
					if   var is 'pT1'    : histo_settings['nbins'] =  8; histo_settings['min'] =  40.; histo_settings['max'] = 200.;
					elif var is 'pT2'    : histo_settings['nbins'] =  6; histo_settings['min'] =  40.; histo_settings['max'] = 100.;

		return histo_settings
