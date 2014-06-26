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

#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/GenEvent.h"

class GenEventBlock : public edm::EDAnalyzer {
public:
    explicit GenEventBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getUntrackedParameter<int>("verbosity", 0)),
        _genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
        _storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
        _pdfWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFWeightsInputTag")) {}

private:
    virtual void endJob() { genEventTree.Write(); }
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);

private:
    int _verbosity;
    const edm::InputTag _genEvtInfoInputTag;
    const bool          _storePDFWeights;
    const edm::InputTag _pdfWeightsInputTag;

    ntuple::GenEventTree genEventTree;
};

void GenEventBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    genEventTree.RunId() = iEvent.id().run();
    genEventTree.LumiBlock() = iEvent.id().luminosityBlock();
    genEventTree.EventId() = iEvent.id().event();

    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(_genEvtInfoInputTag, genEvtInfoProduct);
    if (genEvtInfoProduct.isValid()) {
      edm::LogInfo("GenEventBlock") << "Success >> obtained GenEventInfoProduct for label:" 
                                    << _genEvtInfoInputTag;
      genEventTree.processID() = genEvtInfoProduct->signalProcessID();
      genEventTree.ptHat()     = genEvtInfoProduct->hasBinningValues()
                         ? genEvtInfoProduct->binningValues()[0] : 0.;
    } 
    else {
      edm::LogError("GenEventBlock") << "Error >> Failed to get GenEventInfoProduct for label: " 
                                     << _genEvtInfoInputTag;
    }
    // PDF Weights Part
    if (_storePDFWeights) {
      edm::Handle<std::vector<double> > pdfWeightsHandle;
      iEvent.getByLabel(_pdfWeightsInputTag, pdfWeightsHandle);
      if (pdfWeightsHandle.isValid()) {
	edm::LogInfo("GenEventBlock") << "Success >> obtained PDF handle for label: " 
                                      << _pdfWeightsInputTag;
	for(double weight : *pdfWeightsHandle)
		genEventTree.pdfWeights().push_back(weight);
      } 
      else {
	edm::LogError("GenEventBlock") << "Error >>  Failed to get PDF handle for label: " 
                                       << _pdfWeightsInputTag;
      }
    }
    genEventTree.Fill();
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
