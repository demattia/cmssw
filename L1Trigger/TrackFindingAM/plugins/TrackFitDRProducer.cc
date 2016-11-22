/*! \class   TrackFitDRProducer
 *
 *  \author Marco De Mattia
 *  \date   2016, Jul 25
 *
 */

#ifndef NTupleTools_TrackFitDRProducer_h_

#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "L1Trigger/TrackFindingAM/interface/LinearizedTrackFitter.h"
#include "L1Trigger/TrackFindingAM/interface/DuplicateRemoval.h"
#include "L1Trigger/TrackFindingAM/interface/ParameterDuplicateRemoval.h"
#include<map>


class TrackFitDRProducer : public edm::EDProducer
{
 public:
  explicit TrackFitDRProducer(const edm::ParameterSet&);

 private:

  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

  edm::InputTag TracksTag_;
  std::string CleanedTracksTag_;
  bool parameterBased_;
  slhcl1tt::ParameterDuplicateRemoval parameterDuplicateRemoval_;
  slhcl1tt::DuplicateRemoval duplicateRemoval_;
  unsigned int maxCommonStubs_;
};


/// Begin run
void TrackFitDRProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}


/// End run
void TrackFitDRProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}


TrackFitDRProducer::TrackFitDRProducer(const edm::ParameterSet& iConfig)
{
  TracksTag_ = iConfig.getParameter<edm::InputTag>("TTTrackName");
  CleanedTracksTag_= iConfig.getParameter<std::string>("CleanedTTTrackName");
  parameterBased_ = iConfig.getParameter<bool>("ParameterBased");
  maxCommonStubs_ = iConfig.getParameter<unsigned int>("MaxCommonStubs");
  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( CleanedTracksTag_ );
}


void TrackFitDRProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<TTTrack<Ref_PixelDigi_ > > > TTTrackHandle;
  iEvent.getByLabel(TracksTag_, TTTrackHandle);
  if (parameterBased_) {
    std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > L1TkTracksForOutput(new std::vector< TTTrack< Ref_PixelDigi_ > >(parameterDuplicateRemoval_.ReduceTracks(*TTTrackHandle)));
    iEvent.put(L1TkTracksForOutput, CleanedTracksTag_);
  }
  else {
    std::auto_ptr< std::vector< TTTrack< Ref_PixelDigi_ > > > L1TkTracksForOutput(new std::vector< TTTrack< Ref_PixelDigi_ > >(duplicateRemoval_.CheckTracks(*TTTrackHandle, maxCommonStubs_)));
    iEvent.put(L1TkTracksForOutput, CleanedTracksTag_);
  }
}


DEFINE_FWK_MODULE(TrackFitDRProducer);

#endif


