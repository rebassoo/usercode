import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
#            'file:/uscms_data/d1/rebassoo/Root_files/GamGamWW_AQGC_ScanPoint-a0W-0point0002aCW-0point0008_GEN-FASTSIM-GRun_RECO_withPU.root'
#            'file:/uscms_data/d1/rebassoo/Root_files/GamGamWW_AQGC_ScanPoint-a0W-0point0002aCW-0point0008_START44_Fall11PU_RECO_WithCondor.root'
#'file:/uscms_data/d1/rebassoo/Root_files/GamGamWW_AQGC_ScanPoint-SM_a0W0_aCW0_GEN-FASTSIM-GRun_RECO_withPU.root'
#    'file:/uscms_data/d1/rebassoo/RootFiles/GamGamWW_AnomalousPoint2_Gustavo_START44_Fall11PU_RECO.root'

#'dcache:/pnfs/cms/WAX/11/store/mc/Summer12_DR53X/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/EE7FCC6E-FBE0-E111-97F9-00215E2217BE.root'
#'file:/uscms_data/d2/rebassoo/2012_12_11_MakingFastSim8TeV/FastSimToProduceAODFile/lhe_GEN_FASTSIM_HLT_PU_withWDecays_withTaus.root'
    'file:/uscms_data/d2/rebassoo/2013_2_20_MakingFullSim8TeV/GamGamWW_SM_RECO.root'
    )
)


# particle flow isolation
#
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

#
# rho value for isolation
#
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond

#IS sample MC or Data
ISMC = False

if ISMC:
    process.GlobalTag.globaltag = autoCond['startup']
else:
    process.GlobalTag.globaltag = 'FT_R_53_V6_AN2::All'

process.TFileService = cms.Service("TFileService", fileName = cms.string("SM_WW_8TeV.root") )

#process.load("SlimmedNtuple.Ntupler.ntuplerDATA_cff")

#process.demo = cms.EDAnalyzer('Ntupler',
#          ismc=cms.bool(ISMC)
#)

process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.demo)
#process.p = cms.Path(process.demo)
