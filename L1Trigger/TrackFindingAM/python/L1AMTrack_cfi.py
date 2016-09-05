import FWCore.ParameterSet.Config as cms

# AM-based pattern recognition default sequence
TTPatternsFromStub = cms.EDProducer("TrackFindingAMProducer",
   TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTPatternName      = cms.string("AML1Patterns"),
   inputBankFile      = cms.string('/afs/cern.ch/work/s/sviret/testarea/PatternBanks/BE_5D/Eta7_Phi8/ss32_cov40/612_SLHC6_MUBANK_lowmidhig_sec37_ss32_cov40.pbk'),
   threshold          = cms.int32(5),
   nbMissingHits      = cms.int32(-1)
)

## Trackfit default sequence


TTTCsFromPattern = ( cms.EDProducer("TrackFitTCProducer",
                                    TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                    TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                    TTTrackName        = cms.string("AML1TCs"),
                                    TTTrackBinaryName  = cms.string("AML1BinTCs")
                                    )
                     )

TTCBsFromPattern = ( cms.EDProducer("TrackFitCBProducer",
                                    TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                    TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                    TTTrackName        = cms.string("AML1CBs"),
                                    AdvancedCombinationBuilder = cms.bool(False)
                                    )
                     )

TTTracksINFNFromTC = ( cms.EDProducer("TrackFitPCAProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("MergeTCOutput", "AML1TCs"),
                                  TTTrackName        = cms.string("AML1Tracks"),
                                  verboseLevel       = cms.untracked.int32(1),
                                  fitPerTriggerTower = cms.untracked.bool(False),
                                  removeDuplicates   = cms.untracked.int32(1)
                                  )
                   )

TTTracksTAMUFromTC = ( cms.EDProducer("AMTrackProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("MergeTCOutput", "AML1TCs"),
                                  TTTrackName        = cms.string("AML1Tracks"),
                                  CutOnPrincipals    = cms.bool(False)
                                  )
                   )

TTTracksTAMUFromCB = ( cms.EDProducer("AMTrackProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("MergeCBOutput", "AML1CBs"),
                                  TTTrackName        = cms.string("AML1Tracks"),
                                  CutOnPrincipals    = cms.bool(False)
                                  )
                   )

DuplicateRemovalTAMU = ( cms.EDProducer("TrackFitDRProducer",
                                        TTTrackName        = cms.InputTag("MergeFITCBOutput", "AML1Tracks"),
                                        CleanedTTTrackName = cms.string("CleanedAML1Tracks"),
                                        ParameterBased     = cms.bool(False),
                                        MaxCommonStubs     = cms.uint32(2)
                                        )
                     )


# AM output merging sequence
MergePROutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTPatternsFromStub", "AML1Patterns")),                               
   TTFiltClustersName  = cms.string("ClusInPattern"),
   TTFiltStubsName     = cms.string("StubInPattern"),
   TTPatternsName      = cms.string("AML1Patterns")                         
)

MergeTCOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1TCs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1TCs")                         
)

MergeCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTCBsFromPattern", "AML1CBs")),                               
   TTFiltClustersName  = cms.string("ClusInCB"),
   TTFiltStubsName     = cms.string("StubInCB"),
   TTPatternsName      = cms.string("AML1CBs")                         
)

MergeTCOutputb = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1BinTCs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1BinTCs")                         
)

MergeFITOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksTAMUFromTC", "AML1Tracks")),                               
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("AML1Tracks")                         
)

MergeFITCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksTAMUFromCB", "AML1Tracks")),                               
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("AML1Tracks")                         
)

MergeDROutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPixelDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("DuplicateRemovalTAMU", "CleanedAML1Tracks")),
   TTFiltClustersName  = cms.string("ClusInTrack"),
   TTFiltStubsName     = cms.string("StubInTrack"),
   TTPatternsName      = cms.string("CleanedAML1Tracks")
)
