import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


varOptions = VarParsing('analysis')
varOptions.register("isMC", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.inputFiles = "file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root"
varOptions.parseArguments()
print varOptions

process = cms.Process("MakeHLTIndex")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.autoCond import autoCond

#process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
#process.GlobalTag.globaltag = '94X_dataRun2_v6'

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(varOptions.inputFiles))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

print "### isMC ", varOptions.isMC
print "### InFileName ", process.source.fileNames 

process.makeHLTIndex = cms.EDAnalyzer('MakeHLTIndex',
	triggerResults = cms.InputTag("TriggerResults","","HLT"),
	triggerLabelsName  = cms.vstring(),
	objects = cms.InputTag("selectedPatTrigger"),	

)

import os
HLTLabelNameFile = open(os.environ['CMSSW_BASE'] + '/src/MiniAnalyzer/MiniAnalyzer/data/HLTNames.txt', 'r')
process.makeHLTIndex.triggerLabelsName = HLTLabelNameFile.read().splitlines()
HLTLabelNameFile.close()


process.p = cms.Path(
process.makeHLTIndex
)


