import FWCore.ParameterSet.Config as cms

process = cms.Process('AMTCBASE')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('L1Trigger.TrackFindingAM.L1AMTrack_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
#
# You can use as input file the result of the script AMPR_test.py of part 5.2.2 of the tutorial
#
# Any other EDM file containing patterns and produced with CMSSW 620_SLHC13 should also work
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('INPUTFILENAME'),     
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)


process.TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergeTCOutput", "StubInTC"))
process.TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))
process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag("MergeTCOutput", "AML1TCs"),cms.InputTag("MergeTCOutputb", "AML1BinTCs"))

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('OUTPUTFILENAME'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

# Keep the PR output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR*')

# Keep the TC output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMTC*')
process.RAWSIMoutput.outputCommands.append('drop *_TTTCsFromPattern_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

# Path and EndPath definitions
process.L1AMTC_step         = cms.Path(process.TTTCsFromPatternswStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1AMTC_step,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
