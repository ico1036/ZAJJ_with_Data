import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8') # for mc
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7') # for data


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100  ## --How often you're updated on the progress


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'  
		#'file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root'

	)
)


process.TFileService = cms.Service("TFileService",
      closeFileFast = cms.untracked.bool(True),
      fileName = cms.string("test_GT_DYNt.root")  # for mc
      #fileName = cms.string("test_DoubleEG_GT_Run2016B.root")   # for data


)


process.MiniAnalyzer = cms.EDAnalyzer("MiniAnalyzer",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons")
#    taus = cms.InputTag("slimmedTaus"),
#    photons = cms.InputTag("slimmedPhotons"),
#    jets = cms.InputTag("slimmedJets"),
#    fatjets = cms.InputTag("slimmedJetsAK8"),
#    mets = cms.InputTag("slimmedMETs"),
)

process.p = cms.Path(process.MiniAnalyzer)
