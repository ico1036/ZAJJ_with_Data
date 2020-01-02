import FWCore.ParameterSet.Config as cms
import os
import sys
# set up process
process = cms.Process("EGAMMA")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis') 
options.register('isMiniAOD',True,options.multiplicity.singleton,options.varType.bool," whether we are running on miniAOD or not")
options.register('runVID',True,options.multiplicity.singleton,options.varType.bool," ")

#@#@ Modified (Default = All False) -----------------------------------------------------
options.register('runEnergyCorrections',True,options.multiplicity.singleton,options.varType.bool," ")
# --------------------------------------------------------------------

options.register('applyEnergyCorrections',False,options.multiplicity.singleton,options.varType.bool," ")
options.register('applyVIDOnCorrectedEgamma',False,options.multiplicity.singleton,options.varType.bool," ")
options.register('applyEPCombBug',False,options.multiplicity.singleton,options.varType.bool," ")
options.register('era','2017-Nov17ReReco',options.multiplicity.singleton,options.varType.string," ")

#@#@ Modified (Default = False) -----------------------------------------------------
options.register('isMC',True,options.multiplicity.singleton,options.varType.bool," ")
# --------------------------------------------------------------------


options.register('unscheduled',False,options.multiplicity.singleton,options.varType.bool," ")
options.parseArguments()

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(250),
    limit = cms.untracked.int32(10000000)
)


### --set the number of events
#For Localm turn on
#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(25000)
#)




### -- Input data
## --For local

#inputfile='file:/hcp/data/data02/jwkim2/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver3-v1/50000/B064E00E-AAEB-E611-9D8E-0CC47A7E018E.root'
#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(options.inputFiles),  
#                            fileNames = cms.untracked.vstring(inputFile),  
#                          )

## --For Crab3
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())



era='2016-Legacy'
isMC = False
#outfile = 'DataEDM.root'

def getGlobalTagName(isMC,era):
    if era=='2018-Prompt':
        if isMC: return '102X_upgrade2018_realistic_v15'
        else: return '102X_dataRun2_Prompt_v11'
    elif era=='2017-Nov17ReReco':
        if isMC: return '94X_mc2017_realistic_v10'
        else: return '94X_dataRun2_ReReco_EOY17_v2'
    elif era=='2016-Legacy':
        if isMC: return '94X_mcRun2_asymptotic_v3'
        #if isMC: return '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
        else: return '94X_dataRun2_v10'
        #else: return '80X_dataRun2_2016SeptRepro_v7'
    else:
        raise RuntimeError('Error in runPostRecoEgammaTools, era {} not currently implimented. Allowed eras are 2018-Prompt 2017-Nov17ReReco 2016-Legacy'.format(era)) 
    

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag


#process.GlobalTag = GlobalTag(process.GlobalTag, getGlobalTagName(isMC=options.isMC,era=options.era), '')
process.GlobalTag = GlobalTag(process.GlobalTag, getGlobalTagName(isMC=isMC,era=era), '')


from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       applyEnergyCorrections=options.applyEnergyCorrections,
                       applyVIDOnCorrectedEgamma=options.applyVIDOnCorrectedEgamma,
                       isMiniAOD=options.isMiniAOD, 
                       era=era,                    
                       runVID=options.runVID,
                       runEnergyCorrections=options.runEnergyCorrections,
                       applyEPCombBug=options.applyEPCombBug)


process.p = cms.Path( process.egammaPostRecoSeq )


process.egammaOutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#    fileName = cms.untracked.string(outfile),
    outputCommands = cms.untracked.vstring('drop *',
                                           "keep *_*_*_RECO",
                                           "keep *_*_*_PAT",
                                           'keep *_*_*_HLT',
                                           'keep *_slimmedElectrons*_*_*',
                                           'keep *_slimmedPhotons*_*_*')
                                        )
if not options.isMiniAOD:
    process.egammaOutput.outputCommands = cms.untracked.vstring('drop *',
                                                                'keep *_gedGsfElectrons_*_*',
                                                                'keep *_gedPhotons_*_*',
                                                                'keep *_calibratedElectrons_*_*',
                                                                'keep *_calibratedPhotons_*_*',
                                                                'keep *_egmGsfElectronIDs_*_*',
                                                                'keep *_egmPhotonIDs_*_*')
    

process.outPath = cms.EndPath(process.egammaOutput)

residualCorrFileName = None
if options.isMiniAOD:
    try: 
        residualCorrFileName = process.calibratedPatElectrons.correctionFile.value()
    except AttributeError:
        pass
else:
    try:
        residualCorrFileName = process.calibratedElectrons.correctionFile.value()
    except AttributeError:
        pass

msgStr='''EgammaPostRecoTools:
  running with GT: {}
  running residual E corr: {}'''
print msgStr.format(process.GlobalTag.globaltag.value(),residualCorrFileName)

if options.unscheduled:
    print "  converting to unscheduled"
    from FWCore.ParameterSet.Utilities import convertToUnscheduled
    process=convertToUnscheduled(process)