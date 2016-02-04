#! /usr/bin/python
import selection


def get_histoBins(var, sel = '') :
		histo_settings = {}
		if   var == 'HT'     : histo_settings['nbins'] =  5; histo_settings['min'] = 100.; histo_settings['max'] = 600.;
		elif var == 'MET'    : histo_settings['nbins'] =  6; histo_settings['min'] =   0.; histo_settings['max'] = 120.;
		elif var == 'NJ'     : histo_settings['nbins'] =  4; histo_settings['min'] =   2.; histo_settings['max'] =   6.;
		elif var == 'NbJmed' : histo_settings['nbins'] =  4; histo_settings['min'] =   0.; histo_settings['max'] =   4.;
		elif var == 'pT1'    : histo_settings['nbins'] =  9; histo_settings['min'] =  20.; histo_settings['max'] = 200.;
		elif var == 'pT2'    : histo_settings['nbins'] =  8; histo_settings['min'] =  20.; histo_settings['max'] = 100.;
		elif var == 'Mll'    : histo_settings['nbins'] = 14; histo_settings['min'] =  20.; histo_settings['max'] = 300.;
		elif var == 'NVrtx'  : histo_settings['nbins'] = 40; histo_settings['min'] =   0.; histo_settings['max'] =  40.;
		elif var == 'minMT'  : histo_settings['nbins'] = 20; histo_settings['min'] =   0.; histo_settings['max'] = 400.;
		elif var == 'M3'     : histo_settings['nbins'] = 12; histo_settings['min'] =   0.; histo_settings['max'] = 600.;
		elif var == 'Int'    : histo_settings['nbins'] =  3; histo_settings['min'] =   0.; histo_settings['max'] =   3.;
		elif var == 'CFChan' : histo_settings['nbins'] =  6; histo_settings['min'] =   0.; histo_settings['max'] =   6.;  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		elif var == 'Charge' : histo_settings['nbins'] =  2; histo_settings['min'] =  -2.; histo_settings['max'] =   2.;  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		elif 'PFIso' in var : histo_settings['nbins'] = 20; histo_settings['min'] = 0.; histo_settings['max'] = 1.;
		elif 'D0'    in var : histo_settings['nbins'] = 20; histo_settings['min'] = -0.01; histo_settings['max'] = 0.01;
		elif 'Dz'    in var : histo_settings['nbins'] = 20; histo_settings['min'] = -0.20; histo_settings['max'] = 0.20;

		if sel != '' :
			if sel.name.startswith('1J') :
				if   var == 'NJ'     : histo_settings['nbins'] =  5; histo_settings['min'] =   1.; histo_settings['max'] =   6.;
			if sel.name.startswith('3J1bJ') or sel.name.startswith('final') :
				if   var == 'NJ'     : histo_settings['nbins'] =  3; histo_settings['min'] =   3.; histo_settings['max'] =   6.;
				elif var == 'NbJmed' : histo_settings['nbins'] =  3; histo_settings['min'] =   1.; histo_settings['max'] =   4.;
				if sel.name.startswith('final') :
					if   var == 'pT1'    : histo_settings['nbins'] =  8; histo_settings['min'] =  40.; histo_settings['max'] = 200.;
					elif var == 'pT2'    : histo_settings['nbins'] =  6; histo_settings['min'] =  40.; histo_settings['max'] = 100.;
			if sel.name.startswith('ZElElChMisId') :
				if var == 'Mll'    : histo_settings['nbins'] = 60; histo_settings['min'] =  0.; histo_settings['max'] = 300.;
			if '0J' in sel.name :
				if   var == 'NJ'     : histo_settings['nbins'] =  6; histo_settings['min'] =   0.; histo_settings['max'] =   6.;

		return histo_settings
