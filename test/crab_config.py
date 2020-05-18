
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Crab_ZAJJ_data'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniaodrun_cfg.py'


config.section_("Data")
config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 1
config.Data.lumiMask = '/hcp/data/data02/jwkim2/JSON/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.runRange = '273158-284044' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/jiwoong/'
config.Data.publication = False # wheather upload this to CMSDAS or not 
config.Data.outputDatasetTag = 'Crab_ZAJJ_data'

config.section_('User')
config.section_("Site")
config.Site.storageSite = 'T3_KR_KISTI'
