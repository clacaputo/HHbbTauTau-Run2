/*!
 * \file PatVertexProducer.cc
 * \brief Definition of PatVertexProducer class which produces analysis-level pat::Vertex objects from reco::Vertexes.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-27 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
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

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "HHbbTauTau/PatProduction/interface/PatVertex.h"

class PatVertexProducer : public edm::EDProducer {
public:
    explicit PatVertexProducer(const edm::ParameterSet&);

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::InputTag _inputTag;
};

PatVertexProducer::PatVertexProducer(const edm::ParameterSet& iConfig)
{
    _inputTag = iConfig.getParameter<edm::InputTag>("inputTag");
    produces<pat::VertexCollection>("");
}

void PatVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByLabel(_inputTag, vertexHandle);
    if (!vertexHandle.isValid()) {
        edm::LogError("PatVertexProducer") << "Error: Failed to get VertexCollection for label: " << _inputTag;
        throw std::runtime_error("Failed to get VertexCollection");
    }

    const reco::VertexCollection& recoVertices = *vertexHandle.product();
    std::auto_ptr<pat::VertexCollection> patVertices(new pat::VertexCollection());
    for(const reco::Vertex& recoVertex : recoVertices) {
        const pat::Vertex patVertex(recoVertex);
        patVertices->push_back(patVertex);
    }

    iEvent.put(patVertices);
}

DEFINE_FWK_MODULE(PatVertexProducer);
