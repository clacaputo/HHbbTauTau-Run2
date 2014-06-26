/*!
 * \file VertexBlock.cc
 * \author Original author: Subir Sarkar
 * \author Contributing author: Konstantin Androsov (Siena University, INFN Pisa)
 * \author Contributing author: Maria Teresa Grippo (Siena University, INFN Pisa)
 *
 * Copyright 2011 Subir Sarkar
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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HHbbTauTau/PatProduction/interface/PatVertex.h"

#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/Vertex.h"

class VertexBlock : public edm::EDAnalyzer {
public:
    explicit VertexBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")) {}

private:
    virtual void endJob()
    {
        vertexTree.Write();
    }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    ntuple::VertexTree vertexTree;
};

void VertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    vertexTree.RunId() = iEvent.id().run();
    vertexTree.LumiBlock() = iEvent.id().luminosityBlock();
    vertexTree.EventId() = iEvent.id().event();

    edm::Handle<pat::VertexCollection> primaryVertices;
    iEvent.getByLabel(_inputTag, primaryVertices);

    if (!primaryVertices.isValid()) {
        edm::LogError("VertexBlock") << "Error >> Failed to get VertexCollection for label: "
                                     << _inputTag;
        throw std::runtime_error("Failed to get VertexCollection");
    }
    edm::LogInfo("VertexBlock") << "Total # Primary Vertices: " << primaryVertices->size();

    for (const pat::Vertex& vertex : *primaryVertices) {
        vertexTree.x() = vertex.x();
        vertexTree.y() = vertex.y();
        vertexTree.z() = vertex.z();
        vertexTree.xErr() = vertex.xError();
        vertexTree.yErr() = vertex.yError();
        vertexTree.zErr() = vertex.zError();
        vertexTree.rho() =  vertex.position().rho();
        vertexTree.chi2() = vertex.chi2();
        vertexTree.ndf() = vertex.ndof();
        vertexTree.ntracks() = vertex.tracksSize();
        vertexTree.ntracksw05() = vertex.nTracks(0.5); // number of tracks in the vertex with weight above 0.5
        vertexTree.isfake() = vertex.isFake();
        vertexTree.isvalid() = vertex.isValid();
        vertexTree.sumPt() = vertex.sumPt();
        vertexTree.sumPtSquared() = vertex.sumPtSquared();

//        edm::LogError("VertexBlock") << vertex.hasRefittedTracks() << " " << vertex.nTracks(0.5) << "\n"
//                                     << vertex.sumPt() << " " << vertex.sumPtSquared() << "\n"
//                                     << vertex.numberOfTracks() << " " << vertex.numberOfWeights();

        vertexTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexBlock);
