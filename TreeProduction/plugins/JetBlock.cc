#include <iostream>
#include <algorithm>
#include <boost/shared_ptr.hpp>

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

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
 
#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/Jet.h"

#define SIMPLE_VAR(type, name, default_value) jetTree.name() = jet.bDiscriminator(#name);

namespace {
PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();
}

class JetBlock : public edm::EDAnalyzer {
public:
    explicit JetBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("jetSrc")),
        _jecUncPath(iConfig.getParameter<std::string>("jecUncertainty")),
        _applyResJEC (iConfig.getParameter<bool>     ("applyResidualJEC")),
        _resJEC (iConfig.getParameter<std::string>   ("residualJEC")) {}

private:
    virtual void endJob() { jetTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    std::string _jecUncPath;
    bool _applyResJEC;
    std::string _resJEC;

    ntuple::JetTree jetTree;
};

void JetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    jetTree.EventId() = iEvent.id().event();

    boost::shared_ptr<JetCorrectionUncertainty> jecUnc;
    boost::shared_ptr<JetCorrectorParameters> ResJetCorPar;
    boost::shared_ptr<FactorizedJetCorrector> JEC;
    if (_applyResJEC) {
        edm::FileInPath fipUnc(_jecUncPath);
        jecUnc = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(fipUnc.fullPath()));

        edm::FileInPath fipRes(_resJEC);
        ResJetCorPar = boost::shared_ptr<JetCorrectorParameters>(new JetCorrectorParameters(fipRes.fullPath()));
        std::vector<JetCorrectorParameters> vParam;
        vParam.push_back(*ResJetCorPar);
        JEC = boost::shared_ptr<FactorizedJetCorrector>(new FactorizedJetCorrector(vParam));
    }
    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(_inputTag, jets);

    if (!jets.isValid()) {
        edm::LogError("JetBlock") << "Error >> Failed to get pat::Jet collection for label: "
                                  << _inputTag;
        throw std::runtime_error("Failed to get pat::Jet collection");
    }

    unsigned int njets = jets->size();
    edm::LogInfo("JetBlock") << "Total # PAT Jets: " << njets;
    for (size_t i = 0; i < njets; ++i) {
      const pat::Jet& jet = jets->at(i);

      retpf.set(false);
      jetTree.passLooseID() = pfjetIDLoose(jet, retpf);

      retpf.set(false);
      jetTree.passTightID() = pfjetIDTight(jet, retpf);

      double corr = 1.;
      if (_applyResJEC && iEvent.isRealData() ) {
        JEC->setJetEta(jet.eta());
        JEC->setJetPt(jet.pt()); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      if (jecUnc) {
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt(jet.pt()*corr); // the uncertainty is a function of the corrected pt
      }

      // fill in all the vectors
      jetTree.eta()        = jet.eta();
      jetTree.phi()        = jet.phi();
      jetTree.pt()         = jet.pt()*corr;
      jetTree.pt_raw()     = jet.correctedJet("Uncorrected").pt();
      jetTree.energy()     = jet.energy()*corr;
      jetTree.energy_raw() = jet.correctedJet("Uncorrected").energy();
      jetTree.jecUnc()     = (jecUnc) ? jecUnc->getUncertainty(true) : -1;
      jetTree.resJEC()     = corr;
      jetTree.partonFlavour() = jet.partonFlavour();

      // Jet identification in high pile-up environment
      float mvaValue = jet.userFloat("pileupJetIdProducer:fullDiscriminant");
      jetTree.puIdMVA() = mvaValue;
      jetTree.puIdBits() = jet.userInt("pileupJetIdProducer:fullId"); // Bits: 0:Tight,1:Medium,2:Loose

      jetTree.chargedEmEnergyFraction()     = jet.chargedEmEnergyFraction();
      jetTree.chargedHadronEnergyFraction() = jet.chargedHadronEnergyFraction();
      jetTree.chargedMuEnergyFraction()     = jet.chargedMuEnergyFraction();
      jetTree.electronEnergyFraction()      = jet.electronEnergy()/jet.energy();
      jetTree.muonEnergyFraction()          = jet.muonEnergyFraction();
      jetTree.neutralEmEnergyFraction()     = jet.neutralEmEnergyFraction();
      jetTree.neutralHadronEnergyFraction() = jet.neutralHadronEnergyFraction();
      jetTree.photonEnergyFraction()        = jet.photonEnergyFraction();

      jetTree.chargedHadronMultiplicity()   = jet.chargedHadronMultiplicity();
      jetTree.chargedMultiplicity()         = jet.chargedMultiplicity();
      jetTree.electronMultiplicity()        = jet.electronMultiplicity();
      jetTree.muonMultiplicity()            = jet.muonMultiplicity();
      jetTree.neutralHadronMultiplicity()   = jet.neutralHadronMultiplicity();
      jetTree.neutralMultiplicity()         = jet.neutralMultiplicity();
      jetTree.photonMultiplicity()          = jet.photonMultiplicity();

      jetTree.nConstituents()               = jet.numberOfDaughters();

      B_TAG_DATA()

      jetTree.Fill();
    }
}

#undef SIMPLE_VAR
#undef B_TAG_DATA

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetBlock);
