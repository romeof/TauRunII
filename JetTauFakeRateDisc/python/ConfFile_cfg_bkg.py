import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
#####
##   Modules for the analysis
#####
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#####
##   Interaction with the analyser
#####
process.demo = cms.EDAnalyzer('JetTauFakeRateDisc',
 tau_min_pt  = cms.untracked.double(15),
 tau_max_eta = cms.untracked.double(2.5),
 tau_gentauh_match = cms.bool(False),
 tau_genjet_match  = cms.bool(True)
)
#####
##   Do the Tau
#####
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:/afs/cern.ch/work/f/fromeo/CMSSW_7_2_3/src/TauRunII/JetTauFakeRateDisc/rootfiles/GluGluToHToTauTau_M125_D08.root'
    #'/store/mc/Phys14DR/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/1813E94A-D36E-E411-8EDC-3417EBE34D08.root'
    #'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/125A6B71-C56A-E411-9D2B-0025907609BE.root'
    '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/74A26ADE-806B-E411-9E7B-002590D0B030.root'
    ),
 skipEvents = cms.untracked.uint32(780) #Skip the first n evt, or comment this line if you do not want to skip evt
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) ) #Num of evt to be analysed (whatever is the starting evt)
#####
##   Output file
#####
process.TFileService = cms.Service("TFileService",
 fileName = cms.string('jettaufakeratedisc.root')
)
#####
##   Analysis chain
#####
process.p = cms.Path(process.PFTau*process.demo)
