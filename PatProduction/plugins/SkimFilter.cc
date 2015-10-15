/*!
 * \file SkimFilter.cc
 * \brief Definition of SkimFilter class which applies skim for X->HH->bbTauTau PAT event selection.
 * \author Rosamaria Venditti (INFN Bari, Bari University)
 * \author Contributing author: Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Contributing author: Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-30 created
 *
 * Copyright 2014 Rosamaria Venditti,
 *                Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "HHbbTauTau/PatProduction/interface/PatVertex.h"

class SkimFilter : public edm::EDFilter {
public:
    explicit SkimFilter(const edm::ParameterSet& iConfig)
        : vertexTag(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")),
          muonTag(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
          electronTag(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
          tauTag(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc")) {}

private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);

private:
    edm::InputTag vertexTag;
    edm::InputTag muonTag;
    edm::InputTag electronTag;
    edm::InputTag tauTag;
};

bool SkimFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::VertexCollection> hVertices;
    iEvent.getByLabel(vertexTag, hVertices);
    const pat::VertexCollection& vertices = *hVertices.product();

    edm::Handle<pat::MuonCollection> hMuons;
    iEvent.getByLabel(muonTag, hMuons);
    const pat::MuonCollection& muons = *hMuons.product();

    edm::Handle<pat::ElectronCollection> hElectrons;
    iEvent.getByLabel(electronTag, hElectrons);
    const pat::ElectronCollection& electrons = *hElectrons.product();

    edm::Handle<pat::TauCollection> hTaus;
    iEvent.getByLabel(tauTag, hTaus);
    const pat::TauCollection& taus = *hTaus.product();

    bool haveGoodVertex = false;
    unsigned nGoodTaus = 0;

    for (const pat::Vertex& vertex : vertices){
        if(vertex.ndof() > 4 && std::abs(vertex.z()) < 24 && vertex.position().rho() < 2){
            haveGoodVertex = true;
            break;
        }
    }
    if (!haveGoodVertex) return false;

    for(const pat::Tau& tau : taus) {
        if(tau.pt() > 15.0 && std::abs(tau.eta()) < 2.3  && tau.tauID("decayModeFinding") > 0.5)
            ++nGoodTaus;
        if(nGoodTaus >= 2) return true;
    }
    if(!nGoodTaus) return false;

    for(const pat::Electron& electron : electrons) {
        if(electron.pt() > 15.0 && std::abs(electron.eta()) < 2.5) return true;
    }

    for(const pat::Muon& muon : muons) {
        if(muon.pt() > 15.0 && std::abs(muon.eta())<2.4 && muon.isGlobalMuon()) return true;
    }

    return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimFilter);
