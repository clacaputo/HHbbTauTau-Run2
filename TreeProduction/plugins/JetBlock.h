#ifndef __VHTauTau_TreeMaker_JetBlock_h
#define __VHTauTau_TreeMaker_JetBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class Jet;

class JetBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit JetBlock(const edm::ParameterSet& iConfig);
  virtual ~JetBlock() {}

  enum {
    kMaxJet = 100
  };
private:
  TClonesArray* cloneJet; 
  int fnJet;

  int _verbosity;
  edm::InputTag _inputTag;
  std::string _jecUncPath;
  bool _applyResJEC;
  std::string _resJEC;

  vhtm::Jet* jetB;
};
#endif
