#########################
#
# Configuration file for L1 pattern recognition
# using a pattern bank
#
# This script works on any official production sample
# (assuming that this sample contains a container of TTStubs,
# a container of TTClusters, and a container of TrackingParticles)
#
# Instruction to run this script are provided on this page:
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP V
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 17/02/2014
# Upd.        : 10/02/2016
#
# Script tested with release CMSSW_6_2_X_SLHC
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('AMPR')

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
    input = cms.untracked.int32(20)
)

# Input source
#
# You can use as input file the result of the script SLHC_PGUN_off.py of part 2.2 of the tutorial
#
# Any other EDM file containing stubs and produced with CMSSW 620_SLHC13 and later should also work
# 
# Below some examples are given
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/TTbar/CMSSW_6_2_0_SLHC26-PU_DES23_62_V1_LHCCRefPU200-v1/F24BAA21-8818-E511-9030-0025905A6134.root'), 
                            skipEvents=cms.untracked.uint32(0),
			    duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Some pattern recognition options
process.TTPatternsFromStub.inputBankFile = cms.string('/afs/cern.ch/work/s/sviret/public/PatternBanks/TP2014_extended/620_SLHC7_MUBANK_lowmidhig_sec30_varDC_1000000.pbk')
process.TTPatternsFromStub.threshold     = cms.int32(5)
process.TTPatternsFromStub.nbMissingHits = cms.int32(1)
process.TTPatternsFromStub.debugMode     = cms.int32(0)

# The name of the stub container over which the association is done, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

process.TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergePROutput", "StubInPattern"))
process.TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('AM_output.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)


# Keep the PR output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR')
process.RAWSIMoutput.outputCommands.append('drop *_TTPatternsFromStub_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

# Path and EndPath definitions
process.L1AMPR_step          = cms.Path(process.TTPatternsFromStubswStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1AMPR_step,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions	

