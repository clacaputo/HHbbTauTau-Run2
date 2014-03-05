#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HHbbTauTau/TreeProduction/plugins/GenParticleBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"
#include <limits>

namespace {
size_t FindParticleIndex(const reco::GenParticleCollection& particles, const reco::Candidate* candidate)
    {
        const auto genCandidate = dynamic_cast<const reco::GenParticle*>(candidate);
        if(genCandidate)
        {
            for(size_t index = 0; index < particles.size(); ++index)
            {
                if(&particles[index] == genCandidate)
                    return index;
            }
        }
        throw std::logic_error("bad candidate");
    }
}

GenParticleBlock::GenParticleBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genParticleSrc"))
{
}
GenParticleBlock::~GenParticleBlock() {
}
void GenParticleBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneGenParticle = new TClonesArray("vhtm::GenParticle");
  tree->Branch("GenParticle", &cloneGenParticle, 32000, 2);
  tree->Branch("nGenParticle", &fnGenParticle,  "fnGenParticle/I");
}
void GenParticleBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenParticle->Clear();
  fnGenParticle = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles_handle;
    iEvent.getByLabel(_inputTag, genParticles_handle);
    
    if (genParticles_handle.isValid()) {
	const reco::GenParticleCollection& genParticles = *genParticles_handle.product();
      edm::LogInfo("GenParticleBlock") << "Total # GenParticles: " << genParticles.size();

      unsigned n = 0;
      for (reco::GenParticleCollection::const_iterator it = genParticles.begin(); 
                                                      it != genParticles.end(); ++it, ++n ) {
        /*if (fnGenParticle == kMaxGenParticle) {
	  edm::LogInfo("GenParticleBlock") << "Too many GenParticles, fnGenParticle = " 
                                           << fnGenParticle;
	  break;
        }*/

        // Doat store low energy gluons
        int pdgid     = it->pdgId(); 
        double pt     = it->pt();  
        unsigned numMothers = it->numberOfMothers();
        unsigned numDaughters = it->numberOfDaughters();

        if (numMothers > 2) continue;

        genParticleB = new ((*cloneGenParticle)[fnGenParticle++]) vhtm::GenParticle();

        // fill in all the vectors
        genParticleB->index     = n;
        genParticleB->eta       = it->eta();
        genParticleB->phi       = it->phi();
        genParticleB->p         = it->p();
        genParticleB->px        = it->px();
        genParticleB->py        = it->py();
        genParticleB->pz        = it->pz();
        genParticleB->pt        = pt;
        genParticleB->energy    = it->energy();
        genParticleB->pdgId     = pdgid;
        genParticleB->vx        = it->vx();
        genParticleB->vy        = it->vy();
        genParticleB->vz        = it->vz();
        genParticleB->status    = it->status();
        genParticleB->charge    = it->charge();
        genParticleB->numDaught = numDaughters;
        genParticleB->numMother = numMothers;

        if (numMothers == 0){
            genParticleB->motherIndex_1 = std::numeric_limits<unsigned>::max();
            genParticleB->motherIndex_2 = std::numeric_limits<unsigned>::max();
        }

        if (numMothers == 1){
            genParticleB->motherIndex_1 = FindParticleIndex(genParticles, it->mother(0));;
            genParticleB->motherIndex_2 = std::numeric_limits<unsigned>::max();
        }

        if (numMothers == 2){
            genParticleB->motherIndex_1 = FindParticleIndex(genParticles, it->mother(0));;
            genParticleB->motherIndex_2 = FindParticleIndex(genParticles, it->mother(1));;
        }

	if (_verbosity > 0) std::cout << "pdgId=" << it->pdgId() << std::endl;
        /*for (size_t j = 0; j < it->numberOfMothers(); ++j) {
	  const reco::Candidate* m = it->mother(j);
          for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                          mit != genParticles->end(); ++mit) {
            if (m == &(*mit) ) { // -- keep all the entries && it->pdgId() != mit->pdgId() ) {
	      int idx = std::distance(genParticles->begin(), mit);
	      if (_verbosity > 0) 
                std::cout << "mother index/pdgId: " << idx << "/" << mit->pdgId() << std::endl;
              genParticleB->motherIndices.push_back(idx);
              break;
            }
          } 
        }


      for (size_t j = 0; j < it->numberOfDaughters(); ++j) {
	  const reco::Candidate* d = it->daughter(j);
          for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                          mit != genParticles->end(); ++mit) {
            if (d == &(*mit) ) { // -- keep all the entries  && it->pdgId() != mit->pdgId() ) {
	      int idx = std::distance(genParticles->begin(), mit);
	      if (_verbosity > 0) 
                std::cout << "daughter index/pdgId: " << idx << "/" << mit->pdgId() << std::endl;
              genParticleB->daughtIndices.push_back(idx);
              break;
            }
          } 
        }*/


      }
    } 
    else {
      edm::LogError("GenParticleBlock") << "Error >> Failed to get GenParticleCollection for label: " 
                                        << _inputTag;
    }
  }
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleBlock);
