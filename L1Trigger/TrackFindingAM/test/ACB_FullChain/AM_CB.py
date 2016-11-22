#########################
#
# Configuration file for L1 PCA fit
# using a file with AMTC output 
#
# This script works on any official production sample
# (assuming that this sample contains a container of TTStubs,
# a container of TTClusters, and a container of TrackingParticles)
#
# And of course, a container of TCs.... (TTTracks) 
#
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 04/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('AMCBBASE')

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
#
# You can use as input file the result of the script AMTC_test.py of part 6.1.2 of the tutorial
#
# Any other EDM file containing TCs and produced with CMSSW 620_SLHC27 should also work
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(INPUT_FILE_NAME),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Select which combination builder (default False to use the simple combination builder)
process.TTCBsFromPattern.AdvancedCombinationBuilder   = cms.bool(True)

# Location of the coefficient files
process.TTTracksTAMUFromCB.ConstantsDir               = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")
# Apply cut on the chi2 components
process.TTTracksTAMUFromCB.CutOnPrincipals            = cms.bool(True)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('CB_output_INDEX.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

# Keep the PR output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR')
# process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMTC')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMCBBASE')

# Keep the FIT output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMFIT2')
process.RAWSIMoutput.outputCommands.append('drop *_TTTracks*FromCB_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

# Keep everything
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_*')

# Path and EndPath definitions
process.L1AMCB_FIT_DR_steps = cms.Path(process.TT_CB_FIT_DR_FromPatternswStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1AMCB_FIT_DR_steps,process.endjob_step,process.RAWSIMoutput_step)

# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
