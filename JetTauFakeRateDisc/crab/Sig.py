from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'Tau_GGHTauTau_FligtDistChi2'
config.General.workArea    = 'Tau_GGHTauTau_FligtDistChi2'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = '/afs/cern.ch/work/f/fromeo/CMSSW_7_2_4/src/TauRunII/JetTauFakeRateDisc/python/ConfFile_cfg_sig.py'

config.section_('Data')
config.Data.inputDataset  = '/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 1
config.Data.totalUnits    = 200
config.Data.outLFNDirBase = '/store/user/fromeo/'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Pisa'
#'T2_IT_Bari' 
#'T2_IT_Pisa' 
