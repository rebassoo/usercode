import FWCore.ParameterSet.Config as cms

#process = cms.Process("Demo")

#process.TFileService = cms.Service("TFileService", fileName = cms.string("TestingTrigger.root") )
TFileService = cms.Service("TFileService", fileName = cms.string("SlimmedNtuple.root") )


#IS sample MC or Data
ISMC = False
#process.demo = cms.EDAnalyzer('Ntupler',
demo = cms.EDAnalyzer('Ntupler',
          ismc=cms.bool(ISMC),
          isoValInputTags=cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                        cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                        cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'))
)

#process.runNtupler = cms.Path(process.demo)
