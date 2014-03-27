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

		# FR Predictions from event-by-event weights (pre stored)
		self.npp_mm = 0.; self.npp_staterr_mm = 0.; self.npp_systerr_mm = 0.;
		self.npp_em = 0.; self.npp_staterr_em = 0.; self.npp_systerr_em = 0.;
		self.npp_ee = 0.; self.npp_staterr_ee = 0.; self.npp_systerr_ee = 0.;

		self.npf_mm = 0.; self.npf_staterr_mm = 0.; self.npf_systerr_mm = 0.;
		self.npf_em = 0.; self.npf_staterr_em = 0.; self.npf_systerr_em = 0.;
		self.npf_ee = 0.; self.npf_staterr_ee = 0.; self.npf_systerr_ee = 0.;

		self.nfp_mm = 0.; self.nfp_staterr_mm = 0.; self.nfp_systerr_mm = 0.;
		self.nfp_em = 0.; self.nfp_staterr_em = 0.; self.nfp_systerr_em = 0.;
		self.nfp_ee = 0.; self.nfp_staterr_ee = 0.; self.nfp_systerr_ee = 0.;

		self.nff_mm = 0.; self.nff_staterr_mm = 0.; self.nff_systerr_mm = 0.;
		self.nff_em = 0.; self.nff_staterr_em = 0.; self.nff_systerr_em = 0.;
		self.nff_ee = 0.; self.nff_staterr_ee = 0.; self.nff_systerr_ee = 0.;

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


	def set_observations(self) :
		self.obs    = self.nt2_mm + self.nt2_em + self.nt2_ee
		self.obs_mm = self.nt2_mm
		self.obs_ee = self.nt2_em
		self.obs_em = self.nt2_ee


	def set_fakePredictions(self) :
		self.fake_mm  = self.npf_mm + self.nfp_mm + self.nff_mm
		self.fake_em  = self.npf_em + self.nfp_em + self.nff_em
		self.fake_ee  = self.npf_ee + self.nfp_ee + self.nff_ee
		self.fake     = self.fake_mm + self.fake_em + self.fake_ee


	def __str__(self) :
		observation = 'mm: %3d | em: %3d | ee: %3d' % (self.obs_mm, self.obs_em, self.obs_ee)


	def printout(self) :
		# PRINTOUT
		print "-------------------------------------------------------------------------------------------------------------------------------"
		print "                 |               Mu/Mu               |                E/Mu               |                E/E                ||"
		print "         YIELDS  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  ||"
		print "-------------------------------------------------------------------------------------------------------------------------------"
		print "%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f" % ("Data", self.nt2_mm, self.nt10_mm, self.nt01_mm, self.nt0_mm, self.nt2_em, self.nt10_em, self.nt01_em, self.nt0_em, self.nt2_ee, self.nt10_ee, self.nt01_ee, self.nt0_ee)



		print "  Fake Predictions:"
		print "------------------------------------------------------------------------------------------------------"
		print "                 |            Mu/Mu          |           El/Mu           |            El/El          |"
		print "------------------------------------------------------------------------------------------------------"
		print " Npp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			self.npp_mm, self.npp_staterr_mm, self.npp_systerr_mm,
			self.npp_em, self.npp_staterr_em, self.npp_systerr_em,
			self.npp_ee, self.npp_staterr_ee, self.npp_systerr_ee) 
		print " Npf             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			self.npf_mm, self.npf_staterr_mm, self.npf_systerr_mm,
			self.npf_em, self.npf_staterr_em, self.npf_systerr_em,
			self.npf_ee, self.npf_staterr_ee, self.npf_systerr_ee)
		print " Nfp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			self.nfp_mm, self.nfp_staterr_mm, self.nfp_systerr_mm,
			self.nfp_em, self.nfp_staterr_em, self.nfp_systerr_em,
			self.nfp_ee, self.nfp_staterr_ee, self.nfp_systerr_ee) 
		print " Nff             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			self.nff_mm, self.nff_staterr_mm, self.nff_systerr_mm,
			self.nff_em, self.nff_staterr_em, self.nff_systerr_em,
			self.nff_ee, self.nff_staterr_ee, self.nff_systerr_ee) 
		print "------------------------------------------------------------------------------------------------------"
		print " Total Fakes     |", "%5.1f +/- %5.1f           | %5.1f +/- %5.1f           | %5.1f +/- %5.1f           |" % (
			self.fake_mm, self.fake_err_mm,
			self.fake_em, self.fake_err_em,
			self.fake_ee, self.fake_err_ee) 
		print "------------------------------------------------------------------------------------------------------"
		print " (Value +/- E_stat +/- E_syst) "
		print "//////////////////////////////////////////////////////////////////////////////////////////"


#		print "----------------------------------------------------------------------------------------------"
#		print "       SUMMARY   ||         Mu/Mu         ||         E/Mu          ||          E/E          ||"
#		print "=============================================================================================="
#		print "%16s || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||\n" % ("pred. fakes",
#			self.nF_mm, FR.getMMTotEStat(), self.FakeESyst*nF_mm,
#			self.nF_em, FR.getEMTotEStat(), self.FakeESyst*nF_em,
#			self.nF_ee, FR.getEETotEStat(), self.FakeESyst*nF_ee)
#		print "%16s ||                       || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||\n" % ("pred. chmisid",
#			nt2_em_chmid, nt2_em_chmid_e1, nt2_em_chmid_e2, nt2_ee_chmid, nt2_ee_chmid_e1, nt2_ee_chmid_e2)
#
		print "----------------------------------------------------------------------------------------------"
		print "----------------------------------------------------------------------------------------------"
