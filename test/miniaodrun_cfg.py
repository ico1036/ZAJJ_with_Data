import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("MiniAnalyzer")

### VarParcing
options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.parseArguments()


### Standard modules Loading
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


#@#@ --Added------------------------------------------------
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#-------------------------------------------------------------------


# Input, GT, Output
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.autoCond import autoCond

if (varOptions.isMC):
	#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
	process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'
	inputFileName = 'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'
	outFileName = "DYjet.root"
else:
	#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'  
	process.GlobalTag.globaltag = '94X_dataRun2_v10'  
	inputFileName = 'file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root' 
	outFileName = 'Data.root' 



process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputFileName))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100  ## --How often you're updated on the progress
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))  

print "### isMC ", varOptions.isMC
print "### GlobalTag ", process.GlobalTag.globaltag
print "### InFileName ", process.source.fileNames
print "### OutFileName ", process.TFileService.fileName


#@#@ Added --------------------------------------------------------------------------

## --EGamma ID Sequence
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
						applyEnergyCorrections=False,
						applyVIDOnCorrectedEgamma=False,
						isMiniAOD=True,
						era='2016-Legacy',
						runVID=True,
						runEnergyCorrections=False,
						applyEPCombBug=False)
process.p = cms.Path( process.egammaPostRecoSeq )
# ----------------------------------------------------------------------------------------------------------------


## Input tag============================================

process.MiniAnalyzer = cms.EDAnalyzer("MiniAnalyzer",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
	eleCutTight      = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
	triggerResults = cms.InputTag("TriggerResults","","HLT"),
    #triggerObjects = cms.InputTag("slimmedPatTrigger"),# Low version
    triggerObjects = cms.InputTag("selectedPatTrigger"), # High version
    triggerPrescales = cms.InputTag("patTrigger"),
    triggerIdentifiers = cms.vstring(['HLT_Ele*','HLT_Mu*','HLT_TkMu*']),
    triggerLabelsName  = cms.vstring(),

#    taus = cms.InputTag("slimmedTaus"),
#    photons = cms.InputTag("slimmedPhotons"),
#    jets = cms.InputTag("slimmedJets"),
#    fatjets = cms.InputTag("slimmedJetsAK8"),
#    mets = cms.InputTag("slimmedMETs"),
)

import os
HLTLabelNameFile = open(os.environ['CMSSW_BASE'] + '/src/MiniAnalyzer/MiniAnalyzer/data/HLTNames.txt', 'r')
process.MiniAnalyzer.triggerLabelsName = HLTLabelNameFile.read().splitlines()
HLTLabelNameFile.close()

## Run ===================================
process.p = cms.Path(process.MiniAnalyzer)

### Import Golden JSON file ===========================
#from CondCore.CondDB.CondDB_cfi import *
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename='/hcp/data/data02/jwkim2/JSON/Final/test_DoubleEG_GT_Run2016BJSON.txt').getVLuminosityBlockRange()


