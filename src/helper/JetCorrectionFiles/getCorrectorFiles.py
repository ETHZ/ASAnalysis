import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# define your favorite global tag
## globalTag = 'START53_V7A'
#globalTag = 'START52_V9'
#globalTag = 'START53_V20'
globalTag = 'GR_P_V42_AN2'

process.GlobalTag.globaltag = globalTag+'::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PF    = cms.EDAnalyzer('JetCorrectorDBReader',  
      # below is the communication to the database 
      payloadName    = cms.untracked.string('AK5PF'),
      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
      # but it is recommended to use the GT name that you retrieved the files from.
      globalTag      = cms.untracked.string(globalTag),  
      printScreen    = cms.untracked.bool(False),
      createTextFile = cms.untracked.bool(True)
)
# process.readAK5Calo = process.readAK5PF.clone(payloadName = 'AK5Calo')
# process.readAK5JPT = process.readAK5PF.clone(payloadName = 'AK5JPT')
process.p = cms.Path(process.readAK5PF ) #* process.readAK5Calo * process.readAK5JPT)

