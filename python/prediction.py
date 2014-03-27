#! /usr/bin/python


class prediction :

	def __init__(self) :

		# same-sign tight-tight, tight-loose, loose-tight and loose-loose data yields
		self.nt2_mm  = 0.
		self.nt2_em  = 0.
		self.nt2_ee  = 0.

		self.nt10_mm = 0.
		self.nt10_em = 0.
		self.nt10_ee = 0.

		self.nt01_mm = 0.
		self.nt01_em = 0.
		self.nt01_ee = 0.

		self.nt0_mm  = 0.
		self.nt0_em  = 0.
		self.nt0_ee  = 0.

		# opposite-sign data yields (tight-tight)
		self.nt2_ee_BB_os = 0.
		self.nt2_em_BB_os = 0.
		self.nt2_ee_EE_os = 0.
		self.nt2_em_EE_os = 0.
		self.nt2_ee_EB_os = 0.

		# FR Predictions from event-by-event weights (pre stored)
		self.npp_mm = 0.
		self.npp_em = 0.
		self.npp_ee = 0.

		self.npf_mm = 0.
		self.npf_em = 0.
		self.npf_ee = 0.

		self.nfp_mm = 0.
		self.nfp_em = 0.
		self.nfp_ee = 0.

		self.nff_mm = 0.
		self.nff_em = 0.
		self.nff_ee = 0.

		# all rares, ttW and ttZ yields (tight-tight)
		rares_mm = {}
		rares_em = {}
		rares_ee = {}

		rares_mm_npass = {}
		rares_em_npass = {}
		rares_ee_npass = {}

		# observations
		self.obs    = 0
		self.obs_mm = 0
		self.obs_ee = 0
		self.obs_em = 0

		# predictions
		self.ttw             = 0.
		self.ttw_mm          = 0.
		self.ttw_ee          = 0.
		self.ttw_em          = 0.
		self.ttw_err         = 0.
		self.ttw_err_mm      = 0.
		self.ttw_err_ee      = 0.
		self.ttw_err_em      = 0.
		self.ttw_staterr     = 0.
		self.ttw_staterr_mm  = 0.
		self.ttw_staterr_ee  = 0.
		self.ttw_staterr_em  = 0.
		self.ttw_Nmc         = 0
		self.ttw_Nmc_mm      = 0
		self.ttw_Nmc_ee      = 0
		self.ttw_Nmc_em      = 0
		self.ttw_gen         = 0.

		self.ttw_aMCatNLO     = 0.
		self.ttw_aMCatNLO_gen = 0.

		self.ttz             = 0.
		self.ttz_mm          = 0.
		self.ttz_ee          = 0.
		self.ttz_em          = 0.
		self.ttz_err         = 0.
		self.ttz_err_mm      = 0.
		self.ttz_err_ee      = 0.
		self.ttz_err_em      = 0.
		self.ttz_staterr     = 0.
		self.ttz_staterr_mm  = 0.
		self.ttz_staterr_ee  = 0.
		self.ttz_staterr_em  = 0.
		self.ttz_Nmc         = 0
		self.ttz_Nmc_mm      = 0
		self.ttz_Nmc_ee      = 0
		self.ttz_Nmc_em      = 0

		self.tth             = 0.
		self.tth_mm          = 0.
		self.tth_ee          = 0.
		self.tth_em          = 0.
		self.tth_err         = 0.
		self.tth_err_mm      = 0.
		self.tth_err_ee      = 0.
		self.tth_err_em      = 0.

		self.ttwz            = 0.
		self.ttwz_mm         = 0.
		self.ttwz_ee         = 0.
		self.ttwz_em         = 0.
		self.ttwz_err        = 0.
		self.ttwz_err_mm     = 0.
		self.ttwz_err_ee     = 0.
		self.ttwz_err_em     = 0.
		self.ttwz_staterr    = 0.
		self.ttwz_staterr_mm = 0.
		self.ttwz_staterr_ee = 0.
		self.ttwz_staterr_em = 0.
		self.ttwz_Nmc        = 0
		self.ttwz_Nmc_mm     = 0
		self.ttwz_Nmc_ee     = 0
		self.ttwz_Nmc_em     = 0

		self.fake            = 0.
		self.fake_mm         = 0.
		self.fake_ee         = 0.
		self.fake_em         = 0.
		self.fake_err        = 0.
		self.fake_err_mm     = 0.
		self.fake_err_ee     = 0.
		self.fake_err_em     = 0.

		self.cmid            = 0.
		self.cmid_ee         = 0.
		self.cmid_em         = 0.
		self.cmid_err        = 0.
		self.cmid_err_ee     = 0.
		self.cmid_err_em     = 0.

		self.wz              = 0.
		self.wz_mm           = 0.
		self.wz_ee           = 0.
		self.wz_em           = 0.
		self.wz_err          = 0.
		self.wz_err_mm       = 0.
		self.wz_err_ee       = 0.
		self.wz_err_em       = 0.

		self.rare            = 0.
		self.rare_mm         = 0.
		self.rare_ee         = 0.
		self.rare_em         = 0.
		self.rare_err        = 0.
		self.rare_err_mm     = 0.
		self.rare_err_ee     = 0.
		self.rare_err_em     = 0.

		self.tot             = 0.
		self.tot_mm          = 0.
		self.tot_em          = 0.
		self.tot_ee          = 0.
		self.tot_err         = 0.
		self.tot_err_mm      = 0.
		self.tot_err_em      = 0.
		self.tot_err_ee      = 0.


	def __str__(self) :
		observation = 'mm: %3d | em: %3d | ee: %3d' % (self.obs_mm, self.obs_em, self.obs_ee)
