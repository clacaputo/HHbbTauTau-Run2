#ifndef __VHTauTau_TreeMaker_GenEventBlock_h
#define __VHTauTau_TreeMaker_GenEventBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class GenEvent;

class GenEventBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob() {}

public:
  explicit GenEventBlock(const edm::ParameterSet& iConfig);
  virtual ~GenEventBlock();

private:
  TClonesArray* cloneGenEvent; 
  int _verbosity;

  const edm::InputTag _genEvtInfoInputTag;
  const bool          _storePDFWeights;
  const edm::InputTag _pdfWeightsInputTag;

  std::vector<double> *_pdfWeights;

  vhtm::GenEvent* genEventB;
};
#endif
