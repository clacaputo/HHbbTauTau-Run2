import FWCore.ParameterSet.Config as cms

triggerObjectBlock = cms.EDAnalyzer("TriggerObjectBlock",
  verbosity = cms.int32(0),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  triggerEventTag = cms.InputTag('patTriggerEvent'),
  hltPathsOfInterest = cms.vstring ("HLT_DoubleMu",
                                    "HLT_Mu",
                                    "HLT_IsoMu",
                                    "HLT_TripleMu",
                                    "IsoPFTau",
                                    "TrkIsoT",
                                    "HLT_Ele"),
  May10ReRecoData = cms.bool(False)
)
