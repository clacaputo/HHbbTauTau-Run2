/*!
 * \file GenMETBlock.cc
 * \author Original author: Subir Sarkar
 * \author Contributing author: Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Contributing author: Maria Teresa Grippo (University of Siena, INFN Pisa)
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

#include "HHbbTauTau/TreeProduction/interface/GenMET.h"

class GenMETBlock : public edm::EDAnalyzer {
public:
    explicit GenMETBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("genMETSrc")),
        genMETTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { genMETTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
  int _verbosity;
  edm::InputTag _inputTag;
  ntuple::GenMETTree genMETTree;
};

void GenMETBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    genMETTree.RunId() = iEvent.id().run();
    genMETTree.LumiBlock() = iEvent.id().luminosityBlock();
    genMETTree.EventId() = iEvent.id().event();

    edm::Handle<reco::GenMETCollection> mets;
    iEvent.getByLabel(_inputTag, mets);
    if (!mets.isValid()) {
        edm::LogError("GenMETBlock") << "Error >>  Failed to get GenMETCollection for label: "
                                     << _inputTag;
        throw std::runtime_error("Failed to get GenMETCollection.");
    }
      edm::LogInfo("GenMETBlock") << "Total # GenMETs: " << mets->size();
      for (const reco::GenMET& genMET : *mets) {
        genMETTree.met()    = genMET.pt();
        genMETTree.metphi() = genMET.phi();
        genMETTree.sumet()  = genMET.sumEt();

        genMETTree.Fill();
    } 
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMETBlock);
