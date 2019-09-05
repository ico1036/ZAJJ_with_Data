import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


varOptions = VarParsing('analysis')
varOptions.register("isMC", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
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

if (varOptions.isMC):
    #process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
    inputFileName = 'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'
else:
    #process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
    inputFileName = 'file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root'

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFileName))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

print "### isMC ", varOptions.isMC
print "### InFileName ", process.source.fileNames 

process.makeHLTIndex = cms.EDAnalyzer('MakeHLTIndex',
	triggerResults = cms.InputTag("TriggerResults","","HLT"),
	triggerLabelsName  = cms.vstring(),
	objects = cms.InputTag("selectedPatTrigger"),	

)

import os
#HLTLabelNameFile = open(os.environ['CMSSW_BASE'] + '/src/MiniAnalyzer/MiniAnalyzer/data/HLTNames.txt', 'r')
HLTLabelNameFile = open(os.environ['CMSSW_BASE'] + '/src/MiniAnalyzer/MiniAnalyzer/data/HLTNames.txt.MC', 'r')
process.makeHLTIndex.triggerLabelsName = HLTLabelNameFile.read().splitlines()
HLTLabelNameFile.close()


process.p = cms.Path(
process.makeHLTIndex
)


