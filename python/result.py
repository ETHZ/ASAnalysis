#! /usr/bin/python
import math


class result :

	def __init__(self, chan, charge, chan_str) :

		self.chan     = chan
		self.charge   = charge
		self.chan_str = chan_str

		# same-sign tight-tight, tight-loose, loose-tight and loose-loose data yields
		self.nt2  = 0.
		self.nt10 = 0.
		self.nt01 = 0.
		self.nt0  = 0.

		# FR Predictions from event-by-event weights (pre stored)
		self.npp = 0.; self.npp_staterr = 0.; self.npp_systerr = 0.;
		self.npf = 0.; self.npf_staterr = 0.; self.npf_systerr = 0.;
		self.nfp = 0.; self.nfp_staterr = 0.; self.nfp_systerr = 0.;
		self.nff = 0.; self.nff_staterr = 0.; self.nff_systerr = 0.;

		# all rares, ttW and ttZ yields (tight-tight)   #TODO: REMOVE THESE
		self.rares = {}

#		self.rares_npass = {}

		self.rares_staterr = {}

		# observations
		self.obs    = 0

		# predictions
		self.ttw         = 0.
		self.ttw_err     = 0.
		self.ttw_staterr = 0.
		self.ttw_Nmc     = 0
#		self.ttw_gen     = 0.

#		self.ttw_aMCatNLO     = 0.
#		self.ttw_aMCatNLO_gen = 0.

		self.ttz             = 0.
		self.ttz_err         = 0.
		self.ttz_staterr     = 0.
		self.ttz_Nmc         = 0

		self.tth             = 0.
		self.tth_err         = 0.

		self.ttwz            = 0.
		self.ttwz_err        = 0.
		self.ttwz_staterr    = 0.
		self.ttwz_Nmc        = 0

		self.fake            = 0.
		self.fake_err        = 0.

		self.cmid            = 0.
		self.cmid_err        = 0.

		self.wz              = 0.
		self.wz_err          = 0.

		self.rare            = 0.
		self.rare_err        = 0.

		self.tot             = 0.
		self.tot_err         = 0.


#	def init_dict(self, dict, value) :
#		dict['all'] = value
#		dict['mm']  = value
#		dict['em']  = value
#		dict['ee']  = value


	def set_observations(self) :
		self.obs    = self.nt2


	def set_fakePredictions(self) :
		self.fake  = self.npf + self.nfp + self.nff


	def set_ttwzPredictions(self) :
		self.ttwz         = self.ttw + self.ttwz
		self.ttwz_err     = math.sqrt(self.ttw_err*self.ttw_err + self.ttz_err*self.ttz_err)
		self.ttwz_staterr = math.sqrt(self.ttw_staterr*self.ttw_staterr + self.ttz_staterr*self.ttz_staterr)
		self.ttwz_Nmc     = self.ttw_Nmc + self.ttz_Nmc
#
#
#	def __str__(self) :
#		observation = 'mm: %3d | em: %3d | ee: %3d' % (self.obs_mm, self.obs_em, self.obs_ee)
