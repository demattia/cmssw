############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")
 
 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
Source_Files = cms.untracked.vstring(
    ## ttbar PU=140
    # '/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/TTbar/CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/FC5CDAC1-4E2F-E511-8085-0026189438D9.root'
    'file:/home/demattia/Seb/CMSSW_6_2_0_SLHC28_patch1/src/L1Trigger/TrackFindingAM/test/AMFIT_output.root'
    ## single muons PU=140
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt10_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/067B7C51-3B2F-E511-B41F-0025905A607E.root'
	)
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('TTbar_PU140.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# Path definitions & schedule
############################################################

#run the tracking (example for tracklet method)
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
#process.TTTracksFromPixelDigis.phiWindowSF = cms.untracked.double(2.0)  ## uncomment this to run with wider projection windows (for electrons)
# process.TT_step = cms.Path(process.TrackTriggerTTTracks)
# process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# Valid options are:
#      single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets from ttbar = 6
#      pions from taus = 15
#      inclusively, store all TPs (also those not from primary interaction, if available in samples you are running on) = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(True),       # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(5),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       # L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks"),               # TTTrack input
                                       # MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input 
                                       # TCs
                                       # L1TrackInputTag = cms.InputTag("MergeTCOutput", "AML1TCs"),
                                       # MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "AML1TCs"),
                                       # Tracks
                                       L1TrackInputTag = cms.InputTag("MergeFITOutput", "AML1Tracks"),
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "AML1Tracks"),
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

############################################################
# output module (can use this to store edm-file instead of root-ntuple)
############################################################

#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("FileOut.root"),
#                                fastCloning = cms.untracked.bool( False ),
#                                outputCommands = cms.untracked.vstring('drop *',
#                                                                       'keep *_*_Level1TTTracks_*',
#                                                                       'keep *_*_StubAccepted_*',
#                                                                       'keep *_*_ClusterAccepted_*',
#                                                                       'keep *_*_MergedTrackTruth_*')
#)
#process.FEVToutput_step = cms.EndPath(process.out)


############################################################
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
############################################################

# OLD, FOR TTI SAMPLES
#from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
#process = cust_2023TTI(process)

# NEW, FOR SCOPE DOCUMENT REL VAL SAMPLES
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D
process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

############################################################
# process schedule
############################################################

# process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.ana)
process.schedule = cms.Schedule(process.ana)
