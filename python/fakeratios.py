import ROOT
import numpy as np
import math


class FakeRatios(object) :


	def __init__(self) :
		self.fVerbose  = 0
		self.fNToyMCs  = 100
		self.fAddESyst = 0.
		self.fIsMC     = False
		self.fNGen     = 0

		self.fMMNtl = -1. * np.ones(4)
		self.fEENtl = -1. * np.ones(4)
		self.fEMNtl = -1. * np.ones(4)

		self.fMFRatio = np.zeros(2)
		self.fMPRatio = np.array([1., 0.])
		self.fEFRatio = np.zeros(2)
		self.fEPRatio = np.array([1., 0.])


	def setNToyMCs (self, n) : self.fNToyMCs = n
	def setVerbose (self, n) : self.fVerbose = n
	def setAddESyst(self, e) : self.fAddESyst = e
	def setIsMC    (self, b) : self.fIsMC = b
	def setNGen    (self, n) : self.fNGen = n


	#########
	# Input #
	#########

	def setMMNtl(self, Ntt, Ntl, Nlt, Nll) :
		self.fMMNtl[0] = Ntt
		self.fMMNtl[1] = Ntl
		self.fMMNtl[2] = Nlt
		self.fMMNtl[3] = Nll


	def setEENtl(self, Ntt, Ntl, Nlt, Nll) :
		self.fEENtl[0] = Ntt
		self.fEENtl[1] = Ntl
		self.fEENtl[2] = Nlt
		self.fEENtl[3] = Nll


	def setEMNtl(self, Ntt, Ntl, Nlt, Nll) :
		self.fEMNtl[0] = Ntt
		self.fEMNtl[1] = Ntl
		self.fEMNtl[2] = Nlt
		self.fEMNtl[3] = Nll


	def setMFRatio(self, ratio, error) :
		self.fMFRatio[0] = ratio
		self.fMFRatio[1] = error


	def setMPRatio(self, ratio, error) :
		self.fMPRatio[0] = ratio
		self.fMPRatio[1] = error


	def setEFRatio(self, ratio, error) :
		self.fEFRatio[0] = ratio
		self.fEFRatio[1] = error


	def setEPRatio(self, ratio, error) :
		self.fEPRatio[0] = ratio
		self.fEPRatio[1] = error


	###########################
	# Output and combinations #
	###########################

	## MM

	def getMMNpp(self) :
		return self.getNpp(           self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNppEStat(self) :
		return self.getNppEStat(      self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNppESyst(self) :
		fromtoys = self.getESystFromToys2(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0], self.fMFRatio[1], self.fMFRatio[1], self.fMPRatio[1], self.fMPRatio[1], self.getNpp)
		addsyst = self.fAddESyst * self.getMMNpp()
		return sqrt(fromtoys + addsyst*addsyst)

	def getMMNppETot(self) :
		return math.sqrt(self.getMMNppEStat()*self.getMMNppEStat() + self.getMMNppESyst()*self.getMMNppESyst())


	def getMMNpf(self) :
		return self.getNpf(           self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNpfEStat(self) :
		return self.getNpfEStat(      self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNpfESyst(self) :
		fromtoys = self.getESystFromToys2(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0], self.fMFRatio[1], self.fMFRatio[1], self.fMPRatio[1], self.fMPRatio[1], self.getNpf)
		addsyst = self.fAddESyst * self.getMMNpf()
		return sqrt(fromtoys + addsyst*addsyst)

	def getMMNpfETot(self) :
		return math.sqrt(self.getMMNpfEStat()*self.getMMNpfEStat() + self.getMMNpfESyst()*self.getMMNpfESyst())


	def getMMNfp(self) :
		return self.getNfp(           self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNfpEStat(self) :
		return self.getNfpEStat(      self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNfpESyst(self) :
		fromtoys = self.getESystFromToys2(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0], self.fMFRatio[1], self.fMFRatio[1], self.fMPRatio[1], self.fMPRatio[1], self.getNfp)
		addsyst = self.fAddESyst * self.getMMNfp()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getMMNfpETot(self) :
		return math.sqrt(self.getMMNfpEStat()*self.getMMNfpEStat() + self.getMMNfpESyst()*self.getMMNfpESyst())


	def getMMNff(self) :
		return self.getNff(           self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNffEStat(self) :
		return self.getNffEStat(      self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMNffESyst(self) :
		fromtoys = self.getESystFromToys2(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0], self.fMFRatio[1], self.fMFRatio[1], self.fMPRatio[1], self.fMPRatio[1], self.getNff)
		addsyst = self.fAddESyst * self.getMMNff()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getMMNffETot(self) :
		return math.sqrt(self.getMMNffEStat()*self.getMMNffEStat() + self.getMMNffESyst()*self.getMMNffESyst())


	# EE

	def getEENpp(self) :
		return self.getNpp(   self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENppEStat(self) :
		return self.getNppEStat(      self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENppESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0], self.fEFRatio[1], self.fEFRatio[1], self.fEPRatio[1], self.fEPRatio[1], self.getNpp)
		addsyst = self.fAddESyst * self.getEENpp()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEENppETot(self) :
		return math.sqrt(self.getEENppEStat()*self.getEENppEStat() + self.getEENppESyst()*self.getEENppESyst())


	def getEENpf(self) :
		return self.getNfp(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENpfEStat(self) :
		return self.getNfpEStat(      self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENpfESyst(self) :
		fromtoys =  self.getESystFromToys2(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0], self.fEFRatio[1], self.fEFRatio[1], self.fEPRatio[1], self.fEPRatio[1], self.getNfp)
		addsyst = self.fAddESyst * self.getEENpf()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEENpfETot(self) :
		return math.sqrt(self.getEENpfEStat()*self.getEENpfEStat() + self.getEENpfESyst()*self.getEENpfESyst())


	def getEENfp(self) :
		return self.getNfp(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENfpEStat(self) :
		return self.getNfpEStat(      self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENfpESyst(self) :
		fromtoys =  self.getESystFromToys2(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0], self.fEFRatio[1], self.fEFRatio[1], self.fEPRatio[1], self.fEPRatio[1], self.getNfp)
		addsyst = self.fAddESyst * self.getEENfp()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEENfpETot(self) :
		return math.sqrt(self.getEENfpEStat()*self.getEENfpEStat() + self.getEENfpESyst()*self.getEENfpESyst())


	def getEENff(self) :
		return self.getNff(   self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENffEStat(self) :
		return self.getNffEStat(      self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEENffESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0], self.fEFRatio[1], self.fEFRatio[1], self.fEPRatio[1], self.fEPRatio[1], self.getNff)
		addsyst = self.fAddESyst * self.getEENff()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEENffETot(self) :
		return math.sqrt(self.getEENffEStat()*self.getEENffEStat() + self.getEENffESyst()*self.getEENffESyst())


	# EM

	def getEMNpp(self) :
		return self.getNpp(           self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNppEStat(self) :
		return self.getNppEStat(      self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNppESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNpp)
		addsyst = self.fAddESyst * self.getEMNpp()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMNppETot(self) :
		return math.sqrt(self.getEMNppEStat()*self.getEMNppEStat() + self.getEMNppESyst()*self.getEMNppESyst())


	def getEMNpf(self) :
		return self.getNpf(           self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNpfEStat(self) :
		return self.getNpfEStat(      self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNpfESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNpf)
		addsyst = self.fAddESyst * self.getEMNpf()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMNpfETot(self) :
		return math.sqrt(self.getEMNpfEStat()*self.getEMNpfEStat() + self.getEMNpfESyst()*self.getEMNpfESyst())


	def getEMNfp(self) :
		return self.getNfp(           self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNfpEStat(self) :
		return self.getNfpEStat(      self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNfpESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNfp)
		addsyst = self.fAddESyst * self.getEMNfp()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMNfpETot(self) :
		return math.sqrt(self.getEMNfpEStat()*self.getEMNfpEStat() + self.getEMNfpESyst()*self.getEMNfpESyst())


	def getEMSingleFakes(self) :
		return self.getEMNpf() + self.getEMNfp()

	def getEMSingleEStat(self) :
		return self.getNfpNpfSumEStat(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMSingleESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNfpNpfSum)
		addsyst = self.fAddESyst * (self.getEMNfp() + self.getEMNpf())
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMSingleETot(self) :
		return math.sqrt(self.getEMSingleEStat()*self.getEMSingleEStat() + self.getEMSingleESyst()*self.getEMSingleESyst())


	def getEMNff(self) :
		return self.getNff(           self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNffEStat(self) :
		return self.getNffEStat(      self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMNffESyst(self) :
		fromtoys = self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNff)
		addsyst = self.fAddESyst * self.getEMNff()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMNffETot(self) :
		return math.sqrt(self.getEMNffEStat()*self.getEMNffEStat() + self.getEMNffESyst()*self.getEMNffESyst())


	###############################
	# Single channel combinations #
	###############################

	def getMMTotFakes(self) :
		return self.getNfpNpfNffSum(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMTotEStat(self) :
		return self.getNfpNpfNffSumEStat(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0])

	def getMMTotESyst(self) :
		fromtoys = math.sqrt(self.getESystFromToys2(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], self.fMFRatio[0], self.fMFRatio[0], self.fMPRatio[0], self.fMPRatio[0], self.fMFRatio[1], self.fMFRatio[1], self.fMPRatio[1], self.fMPRatio[1], self.getNfpNpfNffSum))
		addsyst = self.fAddESyst * self.getMMTotFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEETotFakes(self) :
		return self.getNfpNpfNffSum(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEETotEStat(self) :
		return self.getNfpNpfNffSumEStat(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0])

	def getEETotESyst(self) :
		fromtoys = math.sqrt(self.getESystFromToys2(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], self.fEFRatio[0], self.fEFRatio[0], self.fEPRatio[0], self.fEPRatio[0], self.fEFRatio[1], self.fEFRatio[1], self.fEPRatio[1], self.fEPRatio[1], self.getNfpNpfNffSum))
		addsyst = self.fAddESyst * self.getEETotFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getEMTotFakes(self) :
		return self.getNfpNpfNffSum(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMTotEStat(self) :
		return self.getNfpNpfNffSumEStat(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0])

	def getEMTotESyst(self) :
		fromtoys = math.sqrt(self.getESystFromToys2(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], self.fMFRatio[0], self.fEFRatio[0], self.fMPRatio[0], self.fEPRatio[0], self.fMFRatio[1], self.fEFRatio[1], self.fMPRatio[1], self.fEPRatio[1], self.getNfpNpfNffSum))
		addsyst = self.fAddESyst * self.getEMTotFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)


	#############################
	# Multi channel combination #
	#############################

	def getTotFakes(self) :
		# Simple sum
		return self.getMMTotFakes() + self.getEETotFakes() + self.getEMTotFakes()

	def getTotSingleFakes(self) :
		return self.getMMNpf() + self.getEENpf() + self.getEMNpf() + self.getEMNfp()

	def getTotDoubleFakes(self) :
		return self.getMMNff() + self.getEENff() + self.getEMNff()

	def getTotEStat(self) :
		'''Control yields in different channels are independent!'''
		mm = self.getMMTotEStat()
		ee = self.getEETotEStat()
		em = self.getEMTotEStat()
		return math.sqrt( mm*mm + ee*ee + em*em )

	def getTotSingleEStat(self) :
		'''Control yields in different channels are independent!'''
		mm = self.getMMNpfEStat()
		ee = self.getEENpfEStat()
		em = self.getEMSingleEStat()
		return math.sqrt( mm*mm + ee*ee + em*em )

	def getTotDoubleEStat(self) :
		'''Control yields in different channels are independent!'''
		mm = self.getMMNffEStat()
		ee = self.getEENffEStat()
		em = self.getEMNffEStat()
		return math.sqrt( mm*mm + ee*ee + em*em )

	def getTotESyst(self) :
		'''
		Assume errors of p and f are uncorrelated
		Throw toys in a gaussian around f and p with df and dp as their sigmas
		Distributions for f and p are cut off at 0 and 1
		'''
		f1  = self.fMFRatio[0]
		df1 = self.fMFRatio[1]
		f2  = self.fEFRatio[0]
		df2 = self.fEFRatio[1]
		p1  = self.fMPRatio[0]
		dp1 = self.fMPRatio[1]
		p2  = self.fEPRatio[0]
		dp2 = self.fEPRatio[1]

		if self.fVerbose > 2 : print 'getTotESyst ...'
		rand = ROOT.TRandom3()
		rand.SetSeed()
		f_results = []
		i = 0
		while i < self.fNToyMCs : # vary f1
			f1_v = rand.Gaus(f1, df1)
			if f1_v > 1. or f1_v < 0. : continue # throw again if f<0 or f>1
			j = 0
			while j < self.fNToyMCs : # vary f2
				f2_v = rand.Gaus(f2, df2)
				if f2_v > 1. or f2_v < 0. : continue # throw again if f<0 or f>1
				result = self.getTotFakes(f1_v, f2_v, p1, p2)
				f_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_f = np.std(f_results)
		if self.fVerbose > 2 : print ' RMS = %f' % rms_f
		if self.fVerbose > 2 : print ''
		p_results = []
		i = 0
		while i < self.fNToyMCs : # vary p1
			p1_v = rand.Gaus(p1, dp1)
			if p1_v > 1. or p1_v < 0. : continue # throw again if p<0 or p>1
			j = 0
			while j < self.fNToyMCs : # vary p2
				p2_v = rand.Gaus(p2, dp2)
				if p2_v > 1. or p2_v < 0. : continue # throw again if p<0 or p>1
				result = self.getTotFakes(f1, f2, p1_v, p2_v)
				p_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_p = np.std(p_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_p

		fromtoys = rms_f*rms_f + rms_p*rms_p
		addsyst = self.fAddESyst * self.getTotFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getTotSingleESyst(self) :
		'''
		Assume errors of p and f are uncorrelated
		Throw toys in a gaussian around f and p with df and dp as their sigmas
		Distributions for f and p are cut off at 0 and 1
		'''
		f1  = self.fMFRatio[0]
		df1 = self.fMFRatio[1]
		f2  = self.fEFRatio[0]
		df2 = self.fEFRatio[1]
		p1  = self.fMPRatio[0]
		dp1 = self.fMPRatio[1]
		p2  = self.fEPRatio[0]
		dp2 = self.fEPRatio[1]

		if self.fVerbose > 2 : print "getTotESyst ..."
		rand = ROOT.TRandom3()
		rand.SetSeed()
		f_results = []
		i = 0
		while i < self.fNToyMCs : # vary f1
			f1_v = rand.Gaus(f1, df1)
			if f1_v > 1. or f1_v < 0. : continue # throw again if f<0 or f>1
			j = 0
			while j < self.fNToyMCs : # vary f2
				f2_v = rand.Gaus(f2, df2)
				if f2_v > 1. or f2_v < 0. : continue # throw again if f<0 or f>1
				result = self.getTotSingleFakes(f1_v, f2_v, p1, p2)
				f_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_f = np.std(f_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_f
		if self.fVerbose > 2 : print ''
		p_results = []
		i = 0
		while i < self.fNToyMCs : # vary p1
			p1_v = rand.Gaus(p1, dp1)
			if p1_v > 1. or p1_v < 0. : continue # throw again if p<0 or p>1
			j = 0
			while j < self.fNToyMCs : # vary p2
				p2_v = rand.Gaus(p2, dp2)
				if p2_v > 1. or p2_v < 0. : continue # throw again if p<0 or p>1
				result = self.getTotSingleFakes(f1, f2, p1_v, p2_v)
				p_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_p = np.std(p_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_p

		fromtoys = rms_f*rms_f + rms_p*rms_p
		addsyst = self.fAddESyst * self.getTotSingleFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getTotDoubleESyst(self) :
		'''
		Assume errors of p and f are uncorrelated
		Throw toys in a gaussian around f and p with df and dp as their sigmas
		Distributions for f and p are cut off at 0 and 1
		'''
		f1  = self.fMFRatio[0]
		df1 = self.fMFRatio[1]
		f2  = self.fEFRatio[0]
		df2 = self.fEFRatio[1]
		p1  = self.fMPRatio[0]
		dp1 = self.fMPRatio[1]
		p2  = self.fEPRatio[0]
		dp2 = self.fEPRatio[1]

		if self.fVerbose > 2 : print "getTotESyst ..."
		rand = ROOT.TRandom3()
		rand.SetSeed()
		f_results = []
		i = 0
		while i < self.fNToyMCs : # vary f1
			f1_v = rand.Gaus(f1, df1)
			if f1_v > 1. or f1_v < 0. : continue # throw again if f<0 or f>1
			j = 0
			while j < self.fNToyMCs : # vary f2
				f2_v = rand.Gaus(f2, df2)
				if f2_v > 1. or f2_v < 0. : continue # throw again if f<0 or f>1
				result = self.getTotDoubleFakes(f1_v, f2_v, p1, p2)
				f_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_f = np.std(f_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_f
		if self.fVerbose > 2 : print ''
		p_results = []
		i = 0
		while i < self.fNToyMCs : # vary p1
			p1_v = rand.Gaus(p1, dp1)
			if p1_v > 1. or p1_v < 0. : continue # throw again if p<0 or p>1
			j = 0
			while j < self.fNToyMCs : # vary p2
				p2_v = rand.Gaus(p2, dp2)
				if p2_v > 1. or p2_v < 0. : continue # throw again if p<0 or p>1
				result = self.getTotDoubleFakes(f1, f2, p1_v, p2_v)
				p_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_p = np.std(p_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_p

		fromtoys = rms_f*rms_f + rms_p*rms_p
		addsyst = self.fAddESyst * self.getTotDoubleFakes()
		return math.sqrt(fromtoys + addsyst*addsyst)

	def getTotETot(self) :
		return math.sqrt(self.getTotEStat()*self.getTotEStat() + self.getTotESyst()*self.getTotESyst())

	def getTotSingleETot(self) :
		return math.sqrt(self.getTotSingleEStat()*self.getTotSingleEStat() + self.getTotSingleESyst()*self.getTotSingleESyst())

	def getTotDoubleETot(self) :
		return math.sqrt(self.getTotDoubleEStat()*self.getTotDoubleEStat() + self.getTotDoubleESyst()*self.getTotDoubleESyst())


	###########################################
	# Helper methods for inside the toy loops #
	###########################################

	def getTotFakes(self, f1, f2, p1, p2) :
		'''Need this with custom ratios inside the ESyst toy MC loops'''
		mm = self.getNfpNpfNffSum(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], f1, f2, p1, p2)
		ee = self.getNfpNpfNffSum(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], f1, f2, p1, p2)
		em = self.getNfpNpfNffSum(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], f1, f2, p1, p2)
		return mm + ee + em

	def getTotSingleFakes(self, f1, f2, p1, p2) :
		'''Need this with custom ratios inside the ESyst toy MC loops'''
		mm = self.getNfpNpfSum(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], f1, f2, p1, p2)
		ee = self.getNfpNpfSum(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], f1, f2, p1, p2)
		em = self.getNfpNpfSum(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], f1, f2, p1, p2)
		return mm + ee + em

	def getTotDoubleFakes(self, f1, f2, p1, p2) :
		'''Need this with custom ratios inside the ESyst toy MC loops'''
		mm = self.getNff(self.fMMNtl[0], self.fMMNtl[1], self.fMMNtl[2], self.fMMNtl[3], f1, f2, p1, p2)
		ee = self.getNff(self.fEENtl[0], self.fEENtl[1], self.fEENtl[2], self.fEENtl[3], f1, f2, p1, p2)
		em = self.getNff(self.fEMNtl[0], self.fEMNtl[1], self.fEMNtl[2], self.fEMNtl[3], f1, f2, p1, p2)
		return mm + ee + em


	##########
	# ENGINE #
	##########

	# Dilepton formulas from AN-10-261
	#  -- Note: these methods return yields in tight-tight region

	def getNpp(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return p1*p2/((f1-p1)*(f2-p2)) * ( (f1-1.)*(f2-1.)*Ntt
										 + (f1-1.)* f2    *Ntl
										 +  f1    *(f2-1.)*Nlt
										 +  f1    * f2    *Nll )

	def getNpf(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return p1*f2/((f1-p1)*(f2-p2)) * ( (f1-1.)*(1.-p2)*Ntt
										 - (f1-1.)*    p2 *Ntl
										 +  f1    *(1.-p2)*Nlt
										 -  f1    *    p2 *Nll )

	def getNfp(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return f1*p2/((f1-p1)*(f2-p2)) * ( (1.-p1)*(f2-1.)*Ntt
										 + (1.-p1)* f2    *Ntl
										 -     p1 *(f2-1.)*Nlt
										 -     p1 * f2    *Nll )

	def getNff(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return f1*f2/((f1-p1)*(f2-p2)) * ( (1.-p1)*(1.-p2)*Ntt
										 - (1.-p1)*    p2 *Ntl
										 -     p1 *(1.-p2)*Nlt
										 +     p1 *    p2 *Nll )

	def getNfpNpfNffSum(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		'''This is equivalent to p1*f2*Npf + f1*p2*Nfp + f1*f2*Nff'''
		return 1./((f1-p1)*(f2-p2)) * ( (    p1*f2*(f1-1.)*(1.-p2) + f1*p2*(1.-p1)*(f2-1.) + f1*f2*(1.-p1)*(1.-p2))*Ntt
									 +  (-1.*p1*f2*(f1-1.)*   p2   + f1*p2*(1.-p1)* f2     - f1*f2*(1.-p1)*    p2 )*Ntl
									 +  (    p1*f2* f1    *(1.-p2) - f1*p2*    p1 *(f2-1.) - f1*f2*    p1 *(1.-p2))*Nlt
									 +  (-1.*p1*f2* f1    *   p2   - f1*p2*    p1 * f2     + f1*f2*    p1 *    p2 )*Nll )

	def getNfpNpfSum(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		'''This is equivalent to p1*f2*Npf + f1*p2*Nfp'''
		return 1./((f1-p1)*(f2-p2)) * ( (    p1*f2*(f1-1.)*(1.-p2) + f1*p2*(1.-p1)*(f2-1.))*Ntt
									 +  (-1.*p1*f2*(f1-1.)*   p2   + f1*p2*(1.-p1)* f2    )*Ntl
									 +  (    p1*f2* f1    *(1.-p2) - f1*p2*    p1 *(f2-1.))*Nlt
									 +  (-1.*p1*f2* f1    *   p2   - f1*p2*    p1 * f2    )*Nll )


	##############################################################
	# Simple error propagation on above formulas for stat errors #
	##############################################################

	def getNppEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return p1*p2/((f1-p1)*(f2-p2)) * math.sqrt( (f1-1.)*(f2-1.)*(f1-1.)*(f2-1.)*self.getEStat2(Ntt)
											 + (f1-1.)* f2    *(f1-1.)* f2    *self.getEStat2(Ntl)
											 +  f1    *(f2-1.)* f1    *(f2-1.)*self.getEStat2(Nlt)
											 +  f1    * f2    * f1    * f2    *self.getEStat2(Nll) )

	def getNpfEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return p1*f2/((f1-p1)*(f2-p2)) * math.sqrt( (f1-1.)*(1.-p2)*(f1-1.)*(1.-p2)*self.getEStat2(Ntt)
											 + (f1-1.)*    p2 *(f1-1.)*    p2 *self.getEStat2(Ntl)
											 +  f1    *(1.-p2)* f1    *(1.-p2)*self.getEStat2(Nlt)
											 +  f1    *    p2 * f1    *    p2 *self.getEStat2(Nll) )

	def getNfpEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return f1*p2/((f1-p1)*(f2-p2)) * math.sqrt( (1.-p1)*(f2-1.)*(1.-p1)*(f2-1.)*self.getEStat2(Ntt)
											 + (1.-p1)* f2    *(1.-p1)* f2    *self.getEStat2(Ntl)
											 +     p1 *(f2-1.)*    p1 *(f2-1.)*self.getEStat2(Nlt)
											 +     p1 * f2    *    p1 * f2    *self.getEStat2(Nll) )

	def getNffEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return f1*f2/((f1-p1)*(f2-p2)) * math.sqrt( (1.-p1)*(1.-p2)*(1.-p1)*(1.-p2)*self.getEStat2(Ntt)
											 + (1.-p1)*    p2 *(1.-p1)*    p2 *self.getEStat2(Ntl)
											 +     p1 *(1.-p2)*    p1 *(1.-p2)*self.getEStat2(Nlt)
											 +     p1 *    p2 *    p1 *    p2 *self.getEStat2(Nll) )

	def getNfpNpfNffSumEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return 1./((f1-p1)*(f2-p2)) * math.sqrt( (    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1) + f1*f2*(1-p1)*(1-p2))*(    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1) + f1*f2*(1-p1)*(1-p2))*self.getEStat2(Ntt)
										 +  (-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2    - f1*f2*(1-p1)*   p2 )*(-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2    - f1*f2*(1-p1)*   p2 )*self.getEStat2(Ntl)
										 +  (    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1) - f1*f2*   p1 *(1-p2))*(    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1) - f1*f2*   p1 *(1-p2))*self.getEStat2(Nlt)
										 +  (-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2    + f1*f2*   p1 *   p2 )*(-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2    + f1*f2*   p1 *   p2 )*self.getEStat2(Nll) )

	def getNfpNpfSumEStat(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2) :
		return 1./((f1-p1)*(f2-p2)) * math.sqrt( (    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1))*(    p1*f2*(f1-1)*(1-p2) + f1*p2*(1-p1)*(f2-1))*self.getEStat2(Ntt)
										 +  (-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2   )*(-1.*p1*f2*(f1-1)*   p2  + f1*p2*(1-p1)* f2   )*self.getEStat2(Ntl)
										 +  (    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1))*(    p1*f2* f1   *(1-p2) - f1*p2*   p1 *(f2-1))*self.getEStat2(Nlt)
										 +  (-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2   )*(-1.*p1*f2* f1   *   p2  - f1*p2*   p1 * f2   )*self.getEStat2(Nll) )


	###############################################################
	# Throw toy MC to determine systematics from errors on ratios #
	###############################################################

	def getESystFromToys2(self, Ntt, Ntl, Nlt, Nll, f1, f2, p1, p2, df1, df2, dp1, dp2, func) :
		'''
		Assume errors of p and f are uncorrelated
		Throw toys in a gaussian around f and p with df and dp as their sigmas
		Distributions for f and p are cut off at 0 and 1
		'''
		if self.fVerbose > 2 : print "getESystFromToys2 ..."
		rand = ROOT.TRandom3()
		rand.SetSeed()
		f_results = []
		i = 0
		while i < self.fNToyMCs : # vary f1
			f1_v = rand.Gaus(f1, df1)
			if f1_v > 1. or f1_v < 0. : continue # throw again if f<0 or f>1
			j = 0
			while j < self.fNToyMCs : # vary f2
				f2_v = rand.Gaus(f2, df2)
				if f2_v > 1. or f2_v < 0. : continue # throw again if f<0 or f>1
				result = func(Ntt, Ntl, Nlt, Nll, f1_v, f2_v, p1, p2)
				f_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_f = np.std(f_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_f
		if self.fVerbose > 2 : print ''
		p_results = []
		i = 0
		while i < self.fNToyMCs : # vary p1
			p1_v = rand.Gaus(p1, dp1)
			if p1_v > 1. or p1_v < 0. : continue # throw again if p<0 or p>1
			j = 0
			while j < self.fNToyMCs : # vary p2
				p2_v = rand.Gaus(p2, dp2)
				if p2_v > 1. or p2_v < 0. : continue # throw again if p<0 or p>1
				result = func(Ntt, Ntl, Nlt, Nll, f1, f2, p1_v, p2_v)
				p_results.append(result)
				if self.fVerbose > 2 : print result
				j += 1
			i += 1
		rms_p = np.std(p_results)
		if self.fVerbose > 2 : print " RMS = %f" % rms_p

		return rms_f*rms_f + rms_p*rms_p


	##########################
	# Event by event weights #
	##########################

	def getWpp(self, cat, f1, f2, p1, p2) :
		l = p1*p2/((f1-p1)*(f2-p2))
		if cat == 0 : return l*(f1-1.)*(f2-1.)
		if cat == 1 : return l*(f1-1.)* f2
		if cat == 2 : return l* f1*    (f2-1.)
		if cat == 3 : return l* f1*     f2
		return 0.

	def getWpf(self, cat, f1, f2, p1, p2) :
		l = p1*f2/((f1-p1)*(f2-p2))
		if cat == 0 : return l*(f1-1.)*(1.-p2)
		if cat == 1 : return l*(1.-f1)*    p2
		if cat == 2 : return l* f1*    (1.-p2)
		if cat == 3 : return l* f1*  (-1.)*p2
		return 0.

	def getWfp(self, cat, f1, f2, p1, p2) :
		l = f1*p2/((f1-p1)*(f2-p2))
		if cat == 0 : return l*(1.-p1)*(f2-1.)
		if cat == 1 : return l*(1.-p1)* f2
		if cat == 2 : return l*    p1 *(f2-1.)*(-1.)
		if cat == 3 : return l*    p1 * f2    *(-1.)
		return 0.

	def getWff(self, cat, f1, f2, p1, p2) :
		l = f1*f2/((f1-p1)*(f2-p2))
		if cat == 0 : return l*(1.-p1)*(1.-p2)
		if cat == 1 : return l*(1.-p1)*    p2 *(-1.)
		if cat == 2 : return l*    p1 *(1.-p2)*(-1.)
		if cat == 3 : return l*    p1 *    p2
		return 0.


	########################################################################
	# This is the central place to fix how to calculate statistical errors #
	########################################################################

	def getEStat2(self, N) :
		if self.fIsMC and self.fNGen > 0 :
			# If n passed of ngen generated, what is upper limit
			# on number of events passing?
			eff = ROOT.TEfficiency()
			upper = eff.ClopperPearson(self.fNGen, N, 0.68, true)
			delta = upper - float(N)/float(self.fNGen)
			err = delta * float(self.fNGen)
			return err * err

		if N > 3. : return N
		# Get these limits from TMath::ChisquaredQuantile() according to formulas
		# 32.51a and 32.51b in the PDG:
		# nulo(N) = 1/2 TMath::ChisquaredQuantile(1-0.6827, 2N)
		# nuup(N) = 1/2 TMath::ChisquaredQuantile(0.6827, 2(N+1))
		up = 0.5 * ROOT.TMath.ChisquareQuantile(0.6827, 2*(N+1))
		# lo = 0.5 * TMath::ChisquareQuantile(1.-0.6827, 2*N)

		# Always just return the upper one (will always be larger)
		return ((up-N)*(up-N))
