#ifndef __VHTauTau_TreeMaker_TrackBlock_h
#define __VHTauTau_TreeMaker_TrackBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class Track;

class TrackBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}

public:
  explicit TrackBlock(const edm::ParameterSet& iConfig);
  virtual ~TrackBlock() {}

  enum {
    kMaxTrack = 200
  };

private:
  TClonesArray* cloneTrack; 
  int  fnTrack;

  int _verbosity;

  edm::InputTag _inputTag;
  edm::InputTag _beamSpot;

  vhtm::Track* trackB;
};
#endif
