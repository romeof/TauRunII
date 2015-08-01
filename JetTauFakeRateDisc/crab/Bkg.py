from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Tau_QCD_AEIP1D_3'
config.General.workArea    = 'Tau_QCD_AEIP1D_3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = '/afs/cern.ch/work/f/fromeo/CMSSW_7_2_4/src/TauRunII/JetTauFakeRateDisc/python/ConfFile_cfg_bkg.py'

config.section_("Data")
config.Data.inputDataset  = '/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 1
config.Data.totalUnits    = 300
config.Data.outLFNDirBase = '/store/user/fromeo/'

config.section_("Site")
config.Site.storageSite = 'T2_IT_Pisa' 
