#ifndef NTupleTools_AMTrackProducer_h_

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
#include<map>


class AMTrackProducer : public edm::EDProducer
{
 public:
  explicit AMTrackProducer(const edm::ParameterSet &);

 private:

  virtual void beginRun(const edm::Run &run, const edm::EventSetup &iSetup);

  virtual void endRun(const edm::Run &run, const edm::EventSetup &iSetup);

  virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup);

  double mMagneticField;
  const StackedTrackerGeometry *theStackedTracker;
  edm::InputTag StubsTag_;
  edm::InputTag RoadsTag_;
  std::string TracksTag_;
  bool cutOnPrincipals_;
  std::string constantsDir_;
  std::shared_ptr<LinearizedTrackFitter> linearizedTrackFitter_;
};


/// Begin run
void AMTrackProducer::beginRun(const edm::Run &run, const edm::EventSetup &iSetup)
{
  /// Get the geometry references
  edm::ESHandle <StackedTrackerGeometry> StackedTrackerGeomHandle;
  iSetup.get<StackedTrackerGeometryRecord>().get(StackedTrackerGeomHandle);
  theStackedTracker = StackedTrackerGeomHandle.product();

  /// Get magnetic field
  edm::ESHandle <MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField *theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0, 0, 0)).z();
  mMagneticField = (floor(mMagneticFieldStrength * 10.0 + 0.5)) / 10.0;
}


/// End run
void AMTrackProducer::endRun(const edm::Run &run, const edm::EventSetup &iSetup) {}


AMTrackProducer::AMTrackProducer(const edm::ParameterSet &iConfig) :
    StubsTag_(iConfig.getParameter<edm::InputTag>("TTInputStubs")),
    RoadsTag_(iConfig.getParameter<edm::InputTag>("TTInputPatterns")),
    TracksTag_(iConfig.getParameter<std::string>("TTTrackName")),
    cutOnPrincipals_(iConfig.getParameter<bool>("CutOnPrincipals"))
{
  produces<std::vector<TTTrack<Ref_PixelDigi_> > >(TracksTag_);

  edm::FileInPath fp = iConfig.getParameter<edm::FileInPath>("ConstantsDir");
  constantsDir_ = fp.fullPath();
  // Remove the file part to get the constants dir
  constantsDir_ = constantsDir_.substr(0, constantsDir_.rfind("PreEstimate_Transverse/matrixVD_2016.txt"));
  linearizedTrackFitter_ = (std::make_shared<LinearizedTrackFitter>(constantsDir_.c_str(), true, 0, true, 14,
                                                                    cutOnPrincipals_));
}


void AMTrackProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  std::auto_ptr<std::vector<TTTrack<Ref_PixelDigi_> > > L1TkTracksForOutput(new std::vector<TTTrack<Ref_PixelDigi_> >);
  edm::Handle <edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >> TTStubs;
  edm::Handle <std::vector<TTTrack<Ref_PixelDigi_> >> TTRoadHandle;

  iEvent.getByLabel(StubsTag_, TTStubs);
  iEvent.getByLabel(RoadsTag_, TTRoadHandle);

  L1TkTracksForOutput->clear();


  if (TTRoadHandle->size() > 0) {
    unsigned int tkCnt = 0;
    std::vector<TTTrack<Ref_PixelDigi_> >::const_iterator iterTTTrack;

    for (iterTTTrack = TTRoadHandle->begin();
         iterTTTrack != TTRoadHandle->end();
         ++iterTTTrack) {
      edm::Ptr <TTTrack<Ref_PixelDigi_>> tempTrackPtr(TTRoadHandle, tkCnt++);

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();

      std::vector<edm::Ref < edmNew::DetSetVector < TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > > trackStubs = tempTrackPtr->getStubRefs();

      std::vector<double> vars;
      std::vector<int> layers;

      for (unsigned int i = 0; i < trackStubs.size(); i++) {
        edm::Ref <edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_>> tempStubRef = trackStubs.at(i);
        StackedTrackerDetId detIdStub(tempStubRef->getDetId());
        int layer = 0;
        if (detIdStub.isBarrel()) layer = detIdStub.iLayer() + 4;
          // else  continue; //layer = 10+detIdStub.iZ()+abs((int)(detIdStub.iSide())-2)*7;
        else {
          // std::cout << "disk = " << detIdStub.iZ() << std::endl;
          // std::cout << "side = " << detIdStub.iSide() << std::endl;
          // Side is 1 for negative and 2 for positive.
          layer = 10 + detIdStub.iZ() + (2 - detIdStub.iSide()) * 5;
        }
        layers.push_back(layer);
        GlobalPoint posStub = theStackedTracker->findGlobalPosition(&(*tempStubRef));
        vars.push_back(posStub.phi());
        vars.push_back(posStub.perp());
        vars.push_back(posStub.z());
        // std::cout << "z = " << posStub.z() << std::endl;
      }
      //which layers are hit 
      if (layers.size() < 4)continue;
      // int bits=0;
      // if(layers[0]!=5)bits=1;
      // if(layers[0]==5 && layers[1]==7)bits=2;
      // if(layers[1]==6 && layers[2]==8)bits=3;
      // if(layers[3]==9 && layers[2]==7)bits=4;
      // if(layers[3]==8 && layers[4]==10)bits=5;
      // if(layers[4]!=10 && layers.size()<6)bits=6;
      // double normChi2 = linearizedTrackFitter_->fit(vars, bits);
      // std::cout << "layers = " << std::endl;
      // for (auto l : layers) std::cout << l << ", ";
      // std::cout << std::endl;
      double normChi2 = linearizedTrackFitter_->fit(vars, layers);
      // if (layers.size() == 4) {
      // 	if (normChi2 != -1) {
      // 	  std::cout << std::cout << "4-layers track: ";
      // 	  for (auto l : layers) std::cout << l << " ";
      // 	}
      // 	std::cout << std::endl;
      // }
      // chi2/ndf = -1 means the fit did not run on this combination. Either because it has < 5 stubs
      // or beacuse it is a combination for which coefficients are not available.
      if (normChi2 == -1) {
        // std::cout << "layers =" << std::endl;
        // for (auto l : layers) std::cout << l << " ";
        // std::cout << std::endl;
        // std::cout << "R = " << std::endl;
        // for (size_t i=0; i<vars.size()/3; ++i) std::cout << vars[i*3+1] << " ";
        // std::cout << std::endl;
        continue;
      }
      const std::vector<double> &pars = linearizedTrackFitter_->estimatedPars();
      float pt = 1.0 / fabs(pars[0]);
      float px = pt * cos(pars[1]);
      float py = pt * sin(pars[1]);
      float pz = pt * pars[2];
      GlobalVector p3(px, py, pz);

      TTTrack<Ref_PixelDigi_> aTrack(trackStubs);

      aTrack.setSector(seedSector);
      aTrack.setMomentum(p3);
      aTrack.setRInv(0.003 * 3.8114 * pars[0]);
      aTrack.setChi2(normChi2);
      GlobalPoint POCA(0, 0, pars[3]);
      aTrack.setPOCA(POCA);
      aTrack.setRoadRef(tempTrackPtr->roadRef());
      aTrack.setPrincipals(linearizedTrackFitter_->normalizedPrincipalComponents());
      aTrack.setPreEstimatePt(linearizedTrackFitter_->preEstimatedPt());
      aTrack.setLayers(layers);
      L1TkTracksForOutput->push_back(aTrack);
    }
  }

  iEvent.put(L1TkTracksForOutput, TracksTag_);
}


DEFINE_FWK_MODULE(AMTrackProducer);

#endif
