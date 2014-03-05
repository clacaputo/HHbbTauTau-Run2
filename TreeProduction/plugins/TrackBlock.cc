#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "Math/GenVector/VectorUtil.h"

#include "HHbbTauTau/TreeProduction/plugins/TrackBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

TrackBlock::TrackBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
  _beamSpot(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))
{}
void TrackBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneTrack = new TClonesArray("vhtm::Track");
  tree->Branch("Track", &cloneTrack);
  tree->Branch("nTrack", &fnTrack,  "fnTrack/I");
}
void TrackBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTrack->Clear();
  fnTrack = 0;

  // read the beam spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(_beamSpot, beamSpot);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(_inputTag, tracks);

  if (tracks.isValid()) {
    edm::LogInfo("TrackBlock") << "Total # of Tracks: " << tracks->size();
    reco::Track::TrackQuality quality = reco::Track::qualityByName("loose");
    for (reco::TrackCollection::const_iterator  it  = tracks->begin(); 
                                                it != tracks->end(); ++it) {
      const reco::Track &track = *it;
      if (!track.quality(quality)) continue;

      if (fnTrack == kMaxTrack) {
	edm::LogInfo("TrackBlock") << "Too many Tracks in the event, fnTrack = " << fnTrack; 
	break;
      }
      trackB = new ((*cloneTrack)[fnTrack++]) vhtm::Track();
      trackB->eta         = track.eta();
      trackB->etaError    = track.etaError();
      trackB->theta       = track.theta();
      trackB->thetaError  = track.thetaError();
      trackB->phi         = track.phi();
      trackB->phiError    = track.phiError();
      trackB->p           = track.p();
      trackB->pt          = track.pt();
      trackB->ptError     = track.ptError();
      trackB->qoverp      = track.qoverp();
      trackB->qoverpError = track.qoverpError(); 
      trackB->charge      = track.charge();

      trackB->nValidHits  = track.numberOfValidHits();
      trackB->nLostHits   = track.numberOfLostHits();
      trackB->validFraction = track.validFraction();

      const reco::HitPattern& hitp = track.hitPattern();
      trackB->nValidTrackerHits  = hitp.numberOfValidTrackerHits();      
      trackB->nValidPixelHits    = hitp.numberOfValidPixelHits();
      trackB->nValidStripHits    = hitp.numberOfValidStripHits();
      trackB->trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement(); 
      trackB->pixelLayersWithMeasurement   = hitp.pixelLayersWithMeasurement();
      trackB->stripLayersWithMeasurement   = hitp.stripLayersWithMeasurement();
 
      trackB->dxy       = track.dxy(beamSpot->position());
      trackB->dxyError  = track.dxyError();
      trackB->dz        = track.dz(beamSpot->position());
      trackB->dzError   = track.dzError();
     
      trackB->chi2      = track.chi2();
      trackB->ndof      = track.ndof();
      trackB->vx        = track.vx();
      trackB->vy        = track.vy();
      trackB->vz        = track.vz();
    }
  } 
  else {
    edm::LogError("TrackBlock") << "Error! Failed to get reco::Track collection, " 
                                << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackBlock);
