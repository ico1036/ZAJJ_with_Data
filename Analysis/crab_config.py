from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Crab_ZAJJ_data'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniaodrun_cfg.py'

config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.JobType.outputFiles = ['DATA_ntuple.root']
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 1
config.Data.lumiMask = '/hcp/data/data02/jwkim2/JSON/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.runRange = '273158-284044' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/jiwoong/'
config.Data.publication = False # wheather upload this to CMSDAS or not 
config.Data.outputDatasetTag = 'Crab_ZAJJ_data'

config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T2_KR_KISTI'