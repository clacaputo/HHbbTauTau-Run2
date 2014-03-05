#ifndef __VHTauTau_TreeMaker_MuonBlock_h
#define __VHTauTau_TreeMaker_MuonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class Muon;

class MuonBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}

public:
  explicit MuonBlock(const edm::ParameterSet& iConfig);
  virtual ~MuonBlock();


private:
  TClonesArray* cloneMuon; 
  int  fnMuon;

  int _verbosity;

  edm::InputTag _muonInputTag;
  edm::InputTag _vtxInputTag;
  edm::InputTag _beamSpotInputTag;
  bool _beamSpotCorr;
  std::string _muonID;

  vhtm::Muon* muonB;
};
#endif
