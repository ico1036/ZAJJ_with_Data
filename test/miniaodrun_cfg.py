import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("MiniAnalyzer")

### VarParcing
options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.register("isSig", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isSig" )


varOptions.parseArguments()


### Standard modules Loading
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


# Input, GT, Output
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.autoCond import autoCond

if (varOptions.isMC):
	#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
	process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'
	
	if (varOptions.isSig):
		inputFileName = 'file:/xrootd/store/mc/RunIISummer16MiniAODv3/LLAJJ_EWK_MLL-50_MJJ-120_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/110000/FC230F3A-52D2-E811-A2AA-0242AC130002.root'
		outFileName = "LLAJJ_EWK.root"
	
#	else:
#		#inputFileName = 'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'
#		inputFileName = ['file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/46736EB0-5D08-E811-A92F-02163E019B3C.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/B0A89FD3-A508-E811-A443-EC0D9A822626.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/FCD41A19-9F08-E811-809D-0CC47A4D75EE.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/A2B59919-A808-E811-AB77-0CC47A4D76B6.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/E2CA66CB-7309-E811-89D3-0CC47A4C8E5E.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/B00660CD-9609-E811-B80C-1866DA87A7E7.root'
#		,'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/QCDLLAJJ/F2AF1884-9609-E811-B499-0CC47A7C340E.root'
#		]
#		outFileName = "LLAJJ_QCD.root"

	else:
		inputFileName = [#'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/00B01967-E5EC-E611-B179-A0000220FE80.root'
	#	,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/021F3A39-FFEC-E611-8C5C-9CB65404F364.root'
		'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/023F16F7-71ED-E611-9F9F-02163E01A501.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/044E5A6B-E5EC-E611-A7D2-001E67F8F727.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/04BF96BC-DDEB-E611-B2BB-02163E01A843.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/0664AD0A-5AEB-E611-A94C-38EAA7A6D65C.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/0893E9F1-71ED-E611-8AD7-02163E019DE9.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/0C4FC5FE-61EB-E611-81CF-9CB654AD7670.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/14E0B6FB-61EB-E611-B929-001E6739B871.root'
#		,'file:/xrootd/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/16B7270A-5AEB-E611-BEDB-001E67F8FA15.root'
		]
		outFileName="DYjets_check.root"


else:
	#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'  
	process.GlobalTag.globaltag = '94X_dataRun2_v10'  

	inputFileName = 'file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root' 
	outFileName = 'Data.root' 




process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputFileName))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))  # If you want to limit change max evts number
process.MessageLogger.cerr.FwkReport.reportEvery = 100  ## --How often you're updated on the progress
process.TFileService = cms.Service("TFileService",fileName = cms.string(outFileName))  

print "### isMC ", varOptions.isMC
print "### GlobalTag ", process.GlobalTag.globaltag
print "### InFileName ", process.source.fileNames
print "### OutFileName ", process.TFileService.fileName



## --EGamma ID Sequence ------------------------------------------------------
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
						applyEnergyCorrections=False,
						applyVIDOnCorrectedEgamma=False,
						isMiniAOD=True,
						era='2016-Legacy',
						runVID=True,
						runEnergyCorrections=True,
						applyEPCombBug=False)
process.p = cms.Path( process.egammaPostRecoSeq )
## ----------------------------------------------------------------------------



jetCorrectionsType = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')
if (varOptions.isMC):
   jetCorrectionsType = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsType,
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)



## Input tag============================================
if(varOptions.isSig):
	triggerObj_tag_name = 'slimmedPatTrigger'  #high version
else:
	triggerObj_tag_name = 'selectedPatTrigger' #old version


process.MiniAnalyzer = cms.EDAnalyzer("MiniAnalyzer",
#    muons = cms.InputTag("slimmedMuons"),
	triggerResults = cms.InputTag("TriggerResults","","HLT"),
	triggerObjects = cms.InputTag(triggerObj_tag_name), 
	triggerPrescales = cms.InputTag("patTrigger"),
	triggerIdentifiers = cms.vstring(['HLT_Ele*','HLT_Mu*','HLT_TkMu*']),
	triggerLabelsName  = cms.vstring(),
	vertex          = cms.InputTag("offlineSlimmedPrimaryVertices"),
	pileup = cms.InputTag("slimmedAddPileupInfo"),
	electrons = cms.InputTag("slimmedElectrons"),
	photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    jets_update = cms.InputTag("updatedPatJetsUpdatedJEC"),

#    taus = cms.InputTag("slimmedTaus"),
#    fatjets = cms.InputTag("slimmedJetsAK8"),
#    mets = cms.InputTag("slimmedMETs"),
)

import os
HLTLabelNameFile = open(os.environ['CMSSW_BASE'] + '/src/MiniAnalyzer/MiniAnalyzer/data/HLTNames.txt', 'r')
process.MiniAnalyzer.triggerLabelsName = HLTLabelNameFile.read().splitlines()
HLTLabelNameFile.close()

## Run ===================================
process.p = cms.Path(process.egammaPostRecoSeq*
					process.jecSequence*
					process.MiniAnalyzer)

## Do it after crab
### Import Golden JSON file ===========================
#from CondCore.CondDB.CondDB_cfi import *
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename='/hcp/data/data02/jwkim2/JSON/Final/test_DoubleEG_GT_Run2016BJSON.txt').getVLuminosityBlockRange()


