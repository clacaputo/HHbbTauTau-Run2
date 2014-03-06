#ifndef __VHTauTau_TreeMaker_TauBlock_h
#define __VHTauTau_TreeMaker_TauBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class Tau;

class TauBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  const TransientTrackBuilder* trackBuilder_;
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TauBlock(const edm::ParameterSet& iConfig);

  static const reco::PFJetRef& getJetRef(const reco::PFTau& tau);

private:
  TClonesArray* cloneTau; 
  int  fnTau;
  int _verbosity;
  edm::InputTag _inputTag;
  edm::InputTag _vtxInputTag;
};
#endif
