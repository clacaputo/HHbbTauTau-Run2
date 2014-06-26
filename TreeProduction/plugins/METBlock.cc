/*!
 * \file METBlock.cc
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


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#define SMART_TREE_FOR_CMSSW

#include "HHbbTauTau/TreeProduction/interface/MET.h"

class METBlock : public edm::EDAnalyzer {
public:
    typedef std::shared_ptr<ntuple::METTree> TreePtr;
    typedef std::pair<edm::InputTag, TreePtr> TreeDescriptor;
    typedef std::vector<TreeDescriptor> TreeDescriptorVector;

    explicit METBlock(const edm::ParameterSet& iConfig)
    {
        const edm::ParameterSet& sources = iConfig.getParameterSet("metSrc");
        const auto names = sources.getParameterNames();
        for(const std::string& name : names) {
            const edm::InputTag tag = sources.getParameter<edm::InputTag>(name);
            const TreeDescriptor descriptor(tag, TreePtr(new ntuple::METTree(name)));
            descriptors.push_back(descriptor);
        }
    }

private:
    virtual void endJob()
    {
        for(const auto& descriptor : descriptors)
            descriptor.second->Write();
    }

    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    TreeDescriptorVector descriptors;
};

void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    for(const auto& descriptor : descriptors) {
        const edm::InputTag& inputTag = descriptor.first;
        ntuple::METTree& METTree = *descriptor.second.get();
        METTree.RunId() = iEvent.id().run();
        METTree.LumiBlock() = iEvent.id().luminosityBlock();
        METTree.EventId() = iEvent.id().event();

        edm::Handle<std::vector<pat::MET> > mets;
        iEvent.getByLabel(inputTag, mets);

        if (!mets.isValid()) {
            edm::LogError("METBlock") << "Error: Failed to get METCollection for label: " << inputTag;
            throw std::runtime_error("Failed to get METCollection.");
        }
        edm::LogInfo("METBlock") << "Total # METs: " << mets->size();
        for (const pat::MET& MET : *mets) {
            METTree.pt() = MET.pt();
            METTree.phi() = MET.phi();
            METTree.sumEt()  = MET.sumEt();
            METTree.pt_uncorrected() = MET.uncorrectedPt(pat::MET::uncorrALL);
            METTree.phi_uncorrected() = MET.uncorrectedPhi(pat::MET::uncorrALL);
            METTree.sumEt_uncorrected() = MET.sumEt() - MET.corSumEt(pat::MET::uncorrALL);
            METTree.significanceMatrix() = ntuple::SignificanceMatrixToVector(MET.getSignificanceMatrix());
            METTree.Fill();
        }
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
