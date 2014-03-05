#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
 
#include "TTree.h"
#include "TClonesArray.h"

#include "HHbbTauTau/TreeProduction/plugins/JetBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

// Constructor
JetBlock::JetBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("jetSrc")),
  _jecUncPath(iConfig.getParameter<std::string>("jecUncertainty")),
  _applyResJEC (iConfig.getParameter<bool>     ("applyResidualJEC")),
  _resJEC (iConfig.getParameter<std::string>   ("residualJEC"))
{}
void JetBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneJet = new TClonesArray("vhtm::Jet");
  tree->Branch("Jet", &cloneJet, 32000, 2);
  tree->Branch("nJet", &fnJet, "fnJet/I");
}
void JetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneJet->Clear();
  fnJet = 0;

  JetCorrectionUncertainty *jecUnc = 0;
  JetCorrectorParameters *ResJetCorPar = 0;
  FactorizedJetCorrector *JEC = 0;
  bool applyResJECLocal = _applyResJEC;
  if (_applyResJEC) {
    try {
      edm::FileInPath fipUnc(_jecUncPath);;
      jecUnc = new JetCorrectionUncertainty(fipUnc.fullPath());

      edm::FileInPath fipRes(_resJEC);
      ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
      std::vector<JetCorrectorParameters> vParam;
      vParam.push_back(*ResJetCorPar);
      JEC = new FactorizedJetCorrector(vParam);
    }
    catch (std::exception& ex) {
      edm::LogInfo("JetBlock") << "The following exception occurred:" 
                               << std::endl 
                               << ex.what();
      applyResJECLocal = false; 
    } 
  }
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(_inputTag, jets);

  if (jets.isValid()) {
    unsigned int njets = jets->size();
    edm::LogInfo("JetBlock") << "Total # PAT Jets: " << njets;
    for (size_t i = 0; i < njets; ++i) {
      if (fnJet == kMaxJet) {
	edm::LogInfo("JetBlock") << "Too many PAT Jets, fnJet = " << fnJet; 
	break;
      }
      const pat::Jet& jet = jets->at(i);
      retpf.set(false);
      int passjetLoose = (pfjetIDLoose(jet, retpf)) ? 1 : 0;

      retpf.set(false);
      int passjetTight = (pfjetIDTight(jet, retpf)) ? 1 : 0;

      double corr = 1.;
      if (applyResJECLocal && iEvent.isRealData() ) {
        JEC->setJetEta(jet.eta());
        JEC->setJetPt(jet.pt()); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      if (jecUnc) {
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt(jet.pt()*corr); // the uncertainty is a function of the corrected pt
      }
      jetB = new ((*cloneJet)[fnJet++]) vhtm::Jet();

      // fill in all the vectors
      jetB->eta        = jet.eta();
      jetB->phi        = jet.phi();
      jetB->pt         = jet.pt()*corr;
      jetB->pt_raw     = jet.correctedJet("Uncorrected").pt();
      jetB->energy     = jet.energy()*corr;
      jetB->energy_raw = jet.correctedJet("Uncorrected").energy();
      jetB->jecUnc     = (jecUnc) ? jecUnc->getUncertainty(true) : -1;
      jetB->resJEC     = corr;
      jetB->partonFlavour = jet.partonFlavour();

      // Jet identification in high pile-up environment
      float mvaValue = jet.userFloat("pileupJetIdProducer:fullDiscriminant");
      jetB->puIdMVA = mvaValue;
      jetB->puIdBits = jet.userInt("pileupJetIdProducer:fullId"); // Bits: 0:Tight,1:Medium,2:Loose

      jetB->chargedEmEnergyFraction     = jet.chargedEmEnergyFraction();
      jetB->chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
      jetB->chargedMuEnergyFraction     = jet.chargedMuEnergyFraction();
      jetB->electronEnergyFraction      = jet.electronEnergy()/jet.energy();
      jetB->muonEnergyFraction          = jet.muonEnergyFraction();
      jetB->neutralEmEnergyFraction     = jet.neutralEmEnergyFraction();
      jetB->neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
      jetB->photonEnergyFraction        = jet.photonEnergyFraction();

      jetB->chargedHadronMultiplicity   = jet.chargedHadronMultiplicity();
      jetB->chargedMultiplicity         = jet.chargedMultiplicity();
      jetB->electronMultiplicity        = jet.electronMultiplicity();
      jetB->muonMultiplicity            = jet.muonMultiplicity();
      jetB->neutralHadronMultiplicity   = jet.neutralHadronMultiplicity();
      jetB->neutralMultiplicity         = jet.neutralMultiplicity();
      jetB->photonMultiplicity          = jet.photonMultiplicity();

      jetB->nConstituents               = jet.numberOfDaughters();
      jetB->trackCountingHighEffBTag         = jet.bDiscriminator("trackCountingHighEffBJetTags");
      jetB->trackCountingHighPurBTag         = jet.bDiscriminator("trackCountingHighPurBJetTags");
      jetB->simpleSecondaryVertexHighEffBTag = jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jetB->simpleSecondaryVertexHighPurBTag = jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      jetB->jetProbabilityBTag               = jet.bDiscriminator("jetProbabilityBJetTags");
      jetB->jetBProbabilityBTag              = jet.bDiscriminator("jetBProbabilityBJetTags");
      jetB->combinedSecondaryVertexBTag      = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
      jetB->combinedSecondaryVertexMVABTag   = jet.bDiscriminator("combinedSecondaryVertexMVABJetTags");
      jetB->combinedInclusiveSecondaryVertexBTag = jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
      jetB->combinedMVABTag                  = jet.bDiscriminator("combinedMVABJetTags");
      jetB->passLooseID = passjetLoose;
      jetB->passTightID = passjetTight;
      if (_verbosity > 0) 
      std::cout << "JetBlock: trackCountingHighEffBJetTag = " << jetB->trackCountingHighEffBTag
                << ", trackCountingHighPurBJetTag = " << jetB->trackCountingHighPurBTag
                << ", jetProbabilityBTag = " << jetB->jetProbabilityBTag  
                << ", jetBProbabilityBTag = " << jetB->jetBProbabilityBTag  
                << std::endl;
    }
  } 
  else {
    edm::LogError("JetBlock") << "Error >> Failed to get pat::Jet collection for label: " 
                              << _inputTag;
  }
  if (jecUnc) delete jecUnc;
  if (ResJetCorPar) delete ResJetCorPar;
  if (JEC) delete JEC;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetBlock);
