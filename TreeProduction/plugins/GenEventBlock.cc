/*!
 * \file GenEventBlock.cc
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

#include <algorithm>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/GenEvent.h"

class GenEventBlock : public edm::EDAnalyzer {
public:
    explicit GenEventBlock(const edm::ParameterSet& iConfig) :
        _genEventSrc(iConfig.getParameter<edm::InputTag>("genEventSrc")),
        _lheProductSrc(iConfig.getParameter<edm::InputTag>("lheProductSrc")) {}

private:
    virtual void endJob() { genEventTree.Write(); }
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);

private:
    const edm::InputTag _genEventSrc;
    const edm::InputTag _lheProductSrc;


    ntuple::GenEventTree genEventTree;
};

void GenEventBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    genEventTree.RunId() = iEvent.id().run();
    genEventTree.LumiBlock() = iEvent.id().luminosityBlock();
    genEventTree.EventId() = iEvent.id().event();

    //Embedded part
    edm::Handle<GenFilterInfo> genInfoEmbedded_Handle;
    iEvent.getByLabel(_genEventSrc,genInfoEmbedded_Handle);


    edm::Handle< LHEEventProduct > LHEEvent_Handle;
    iEvent.getByLabel(_lheProductSrc,LHEEvent_Handle);

    if (!genInfoEmbedded_Handle.isValid()) {
        edm::LogError("GenEventBlock") << "Error >> Failed to get GenFilterInfor for label: "
                                     << _genEventSrc;
        genEventTree.embeddedWeight() = std::numeric_limits<float>::lowest();
    }
    else {
        const GenFilterInfo* genFilterInfo = genInfoEmbedded_Handle.product();
        genEventTree.embeddedWeight() = genFilterInfo->filterEfficiency();
    }

    if(!LHEEvent_Handle.isValid()){
        genEventTree.nup() = std::numeric_limits<unsigned>::lowest();
    }
    else {
        const LHEEventProduct* lheEvent = LHEEvent_Handle.product();
        genEventTree.nup() = lheEvent->hepeup().NUP;
    }

    genEventTree.Fill();
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
