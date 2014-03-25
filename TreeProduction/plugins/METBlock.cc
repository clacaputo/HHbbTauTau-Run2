#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HHbbTauTau/TreeProduction/interface/MET.h"

class METBlock : public edm::EDAnalyzer {
public:

    explicit METBlock(const edm::ParameterSet& iConfig) :
      _verbosity(iConfig.getParameter<int>("verbosity")),
      _pfinputTag(iConfig.getParameter<edm::InputTag>("metSrc")),
      _corrinputTag(iConfig.getParameter<edm::InputTag>("corrmetSrc")),
      _mvainputTag(iConfig.getParameter<edm::InputTag>("mvametSrc"))
    {}

private:
    virtual void endJob() { METTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
  int _verbosity;
  edm::InputTag _pfinputTag;
  edm::InputTag _corrinputTag;
  edm::InputTag _mvainputTag;
  ntuple::METTree METTree;
};

void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    METTree.EventId() = iEvent.id().event();

    edm::Handle<std::vector<pat::MET> > mets;
    iEvent.getByLabel(_pfinputTag, mets);

    if (!mets.isValid()) {
        edm::LogError("METBlock") << "Error >>  Failed to get METCollection for label: "
                                     << _pfinputTag;
        throw std::runtime_error("Failed to get METCollection.");
    }
      edm::LogInfo("METBlock") << "Total # METs: " << mets->size();
      for (const pat::MET& MET : *mets) {
        METTree.met()   = MET.pt();
        METTree.metphi() = MET.phi();
        METTree.sumet()  = MET.sumEt();
        METTree.metuncorr() = MET.uncorrectedPt(pat::MET::uncorrALL);
        METTree.metphiuncorr() = MET.uncorrectedPhi(pat::MET::uncorrALL);
        METTree.sumetuncorr() = MET.sumEt() - MET.corSumEt(pat::MET::uncorrALL);

        METTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
