from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Crab_EGamma_data'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runEgammaPostRecoTools.py'

config.Data.inputDataset = '/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.JobType.outputFiles = ['DATA_EDM.root']
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 60
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.runRange = '274441-275376' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/jiwoong/'
config.Data.publication = False # wheather upload this to CMSDAS or not 
config.Data.outputDatasetTag = 'Crab_EGamma_data'

config.Site.storageSite = 'T3_KR_KISTI'

