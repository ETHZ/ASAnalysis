#! /usr/bin/python
import math


class result(object) :

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

		self.fake_err        = 0.
		self.fake_staterr    = 0.

		self.cmid            = 0.
		self.cmid_err        = 0.
		self.cmid_staterr    = 0.

		self.wz              = 0.
		self.wz_err          = 0.
		self.wz_staterr      = 0.

		self.rare            = 0.
		self.rare_err        = 0.
		self.rare_staterr    = 0.

		self.tot_staterr     = 0.

		self.mc              = {}
		self.mc_staterr      = {}


	@property
	def obs(self) :
		return self.nt2


	@property
	def fake(self) :
		fake = self.npf + self.nfp + self.nff
		return fake


	@property
	def ttwz(self) :
		ttwz = self.ttw + self.ttz
		return ttwz


	@property
	def ttwz_err(self) :
		ttwz_err = math.sqrt(self.ttw_err*self.ttw_err + self.ttz_err*self.ttz_err)
		return ttwz_err


	@property
	def ttwz_staterr(self) :
		ttwz_staterr = math.sqrt(self.ttw_staterr*self.ttw_staterr + self.ttz_staterr*self.ttz_staterr)
		return ttwz_staterr


	@property
	def ttwz_Nmc(self) :
		ttwz_Nmc = self.ttw_Nmc + self.ttz_Nmc
		return ttwz_Nmc


	@property
	def tot(self) :
		tot = self.fake + self.cmid + self.rare + self.wz + self.ttz
		return tot


	@property
	def tot_err(self) :
		tot_err = math.sqrt(self.fake_err*self.fake_err + self.cmid_err*self.cmid_err + self.rare_err*self.rare_err + self.wz_err*self.wz_err + self.ttz_err*self.ttz_err)
		return tot_err


	@property
	def tot_exp(self) :
		tot_exp = self.tot + self.ttw
		return tot_exp


	@property
	def tot_exp_err(self) :
		tot_exp_err = math.sqrt(self.tot_err*self.tot_err + self.ttw_statsysterr*self.ttw_statsysterr)
		return tot_exp_err


	@property
	def ttw_statsysterr(self) :
		'''statistical and systematic uncertainty without the cross section uncertainty'''
		ttw_statsysterr = math.sqrt(self.ttw_staterr*self.ttw_staterr + 0.08*0.08*self.ttw*self.ttw)
		return ttw_statsysterr
