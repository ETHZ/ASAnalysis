#! /usr/bin/python


class prediction :

	def __init__(self) :

		self.obs    = 0
		self.obs_mm = 0
		self.obs_ee = 0
		self.obs_em = 0

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
