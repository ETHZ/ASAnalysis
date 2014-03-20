#! /usr/bin/python

import ROOT

class sample :
	'''stores informations and numbers about a sample'''


	def __init__(self, name, datamc, xsec, ngen) :
		self.name   = name
		self.datamc = datamc
		self.xsec   = xsec
		self.ngen   = ngen


	def getLumi(self) :
		'''returns luminosity'''
		if self.datamc is 0 : return -1.
		if self.ngen > 0 and self.xsec > 0 : return float(self.ngen)/self.xsec
		else : return -1.


	def getError(self, n) :
		## If n passed of ngen generated, what is upper limit
		## on number of events passing?
		if self.ngen <= 0 : return -1.
		eff = ROOT.TEfficiency()
		upper = eff.ClopperPearson(self.ngen, n, 0.68, True)
		delta = upper - float(n)/float(self.ngen)
		return delta * float(self.ngen)


	def getError2(self, n) :
		err = getError(n)
		return err*err


##		int getType(){ // -1: undef, 0: data, 1: QCD, 2: top, 3: EWK, 4: Rare SM, 5: diboson
##			if(datamc == 0) return 0;
##			if( (sname.Contains("QCD")) ||
##			    (sname) == "MuEnr15"    ||
##			    (sname) == "MuEnr10"    ||
##			    (sname) == "EMEnr20"    ||
##			    (sname) == "EMEnr30"    ) return 1;
##			if( (sname.Contains("SingleT")) ||
##			    (sname.Contains("TTJets"))  ) return 2;
##			if( (sname.Contains("DYJets")) ||
##			    (sname.Contains("GJets"))  ||
##			    (sname) == "WJets"  ||
##				(sname) == "WJets1" ||
##				(sname) == "WJets2" ||
##				(sname) == "WbbJets" )   return 3;
##			if( (sname) == "HWW"    ||
##			    (sname) == "HZZ"    ||
##			    (sname) == "HTauTau"   ||
##			    (sname) == "TTbarW"    ||
##			    (sname) == "TTbarZ"    ||
##			    (sname) == "TTbarG"    ||
##			    (sname) == "TbZ"       ||
##			    (sname) == "DPSWW"     ||
##			    (sname) == "WWZ"       ||
##			    (sname) == "WZZ"       ||
##			    (sname) == "WZZ"       ||
##			    (sname) == "WWG"       ||
##			    (sname) == "ZZZ"       ||
##			    (sname) == "WWW"       ||
##			    (sname) == "W+W+"      ||
##			    (sname) == "W-W-"      ||
##			    (sname) == "TTbarWW")     return 4;
##			if( (sname.Contains("GVJets"))    || 
##			    (sname.Contains("WGstarMu"))  || 
##			    (sname.Contains("WGstarTau")) ||
##			    (sname.Contains("WWTo2L2Nu")) ||
##			    (sname.Contains("WZTo3LNu"))  ||
##			    (sname.Contains("ZZTo4L")) ) return 5;
##			else {
##				cout << "SSDLDumper::Sample::getType() ==> ERROR: "<< sname << " has no defined type!" << endl;
##				return -1;
##			}
##		}
##		int getProc(){ // used for binned samples
##			if(datamc == 0)                                 return 0;
##			if(sname == "TTJets" )                          return 1;
##			if(sname.Contains("SingleT"))                   return 2;
##			if(sname == "WJets" )                           return 3;
##			if(sname.Contains("DYJets"))                    return 4;
##			if(sname.Contains("GJets"))                     return 5;
##			if(sname.Contains("WWTo2L2Nu"))                 return 6; 
##			if(sname.Contains("WZTo3LNu"))                  return 7; 
##			if(sname.Contains("ZZTo4L"))                    return 8; 
##			if(sname.Contains("GVJets") ||
##			   sname.Contains("WGstarMu") ||
##			   sname.Contains("WGstarTau"))                 return 9; 
##			if(sname == "HWW" || 
##			   sname == "HZZ" || 
##			   sname == "HTauTau")                          return 10;
##			if(sname == "TTbarW")                           return 11;
##			if(sname == "TTbarZ")                           return 12;
##			if(sname == "TTbarG")                           return 13;
##			if(sname == "W+W+" || sname == "W-W-")          return 14;
##			if(sname == "WWZ" ||                                     
##			   sname == "WZZ" ||                                     
##			   sname == "WWG" ||                                     
##			   sname == "ZZZ" ||                                     
##			   sname == "WWW" )                             return 15;
##			if(sname == "DPSWW")                            return 16;
##			if(sname.Contains("QCD") || 
##			   sname == "MuEnr10"    ||
##			   sname == "MuEnr15"    ||
##			   sname == "EMEnr20"    ||
##			   sname == "EMEnr30")                          return 17;
##			if(sname == "TbZ")                              return 18;
##			if(sname == "TTbarWW")                          return 19;
##			if(sname == "TTJets_matchingdown")				return 20;
##			if(sname == "TTJets_matchingup")				return 21;
##			if(sname == "TTJets_scaledown")					return 22;
##			if(sname == "TTJets_scaleup")					return 23;
##			else {
##				cout << "SSDLDumper::Sample::getProc() ==> ERROR: "<< sname << " has no defined process!" << endl;
##				return -1;
##			}
##		}
##		TString getProcName(int proc){ // used for binned samples
##			if(proc ==  0) return "Data";
##			if(proc ==  1) return "$t\\bar{t}$";
##			if(proc ==  2) return "Single t";
##			if(proc ==  3) return "W+jets";
##			if(proc ==  4) return "Z+jets";
##			if(proc ==  5) return "$\\gamma$+jets";
##			if(proc ==  6) return "WW";
##			if(proc ==  7) return "WZ";
##			if(proc ==  8) return "ZZ";
##			if(proc ==  9) return "V$\\gamma$+jets";
##			if(proc == 10) return "$t\\bar{t}$H";
##			if(proc == 11) return "$t\\bar{t}$W";
##			if(proc == 12) return "$t\\bar{t}$Z";
##			if(proc == 13) return "$t\\bar{t}\\gamma$";
##			if(proc == 14) return "W$^{\\pm}$W$^{\\pm}$";
##			if(proc == 15) return "Tri-Boson";
##			if(proc == 16) return "DPS (2$\\times$ W+jets)";
##			if(proc == 17) return "QCD";
##			if(proc == 18) return "TbZ";
##			if(proc == 19) return "$t\\bar{t}$matchingdown";
##			if(proc == 20) return "$t\\bar{t}$matchingup";
##			if(proc == 21) return "$t\\bar{t}$scaledown";
##			if(proc == 22) return "$t\\bar{t}$scaleup";
##			else {
##				cout << "SSDLDumper::Sample::getProcName() ==> ERROR: "<< proc << " has no defined process name!" << endl;
##				return "";
##			}
##		}
##		inline int getNProcs(){return 17;} // make sure this number corresponds to the number of
##		                                   // processes define in the previous method
