import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100  ## --How often you're updated on the progress

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         'file:/hcp/data/data02/jwkim2/store/mc/RunIISummer16MiniAODv2/701979BD-51C4-E611-9FC9-C4346BC84780.root'    )
)

process.TFileService = cms.Service("TFileService",
      closeFileFast = cms.untracked.bool(True),
      fileName = cms.string("test_DYNt.root") 
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
