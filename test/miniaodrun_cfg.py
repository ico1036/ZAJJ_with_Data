import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

### VarParcing
options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.parseArguments()



process = cms.Process("MiniAnalyzer")

### Standard modules Loading
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")



# Configuration options =======================

runOnData=False # or False (when you run on MC)

# Input, GT, Output
from Configuration.AlCa.GlobalTag import GlobalTag

if (varOptions.isMC):
	process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
	inputFileName = 'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'
	outFileName = "DYjet.root"
else:
	process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'  
	inputFileName = 'file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root' 
	outFileName = 'Data.root' 

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputFileName))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100  ## --How often you're updated on the progress
process.TFileService = cms.Service("TFileService",
      closeFileFast = cms.untracked.bool(True),
      fileName = cms.string(outFileName)   
)

## Input tag============================================

process.MiniAnalyzer = cms.EDAnalyzer("MiniAnalyzer",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
#    taus = cms.InputTag("slimmedTaus"),
#    photons = cms.InputTag("slimmedPhotons"),
#    jets = cms.InputTag("slimmedJets"),
#    fatjets = cms.InputTag("slimmedJetsAK8"),
#    mets = cms.InputTag("slimmedMETs"),

)

## Import Golden JSON file ===========================
from CondCore.CondDB.CondDB_cfi import *
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename='/hcp/data/data02/jwkim2/JSON/Final/test_DoubleEG_GT_Run2016BJSON.txt').getVLuminosityBlockRange()


## Run ===================================
process.p = cms.Path(process.MiniAnalyzer)
