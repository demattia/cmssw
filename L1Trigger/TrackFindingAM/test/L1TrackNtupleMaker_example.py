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
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring(
        'file:AMFIT_output.root'
	)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('TrkEff_output.root'), closeFileFast = cms.untracked.bool(True))

############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save all L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(5),         # use 4 or 5-parameter L1 track fit ??
                                       TP_minNStub = cms.int32(4),       # require TP to have >= minNstub (efficiency denominator)
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("MergeFITOutput", "AML1Tracks"),               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "AML1Tracks") ## MCTruth input 
                                       )
process.ana = cms.Path(process.L1TrackNtuple)


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
process = cust_2023TTI(process)

process.schedule = cms.Schedule(process.ana)
