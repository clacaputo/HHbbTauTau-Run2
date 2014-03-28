/*!
 * \file TreeExtractor.h
 * \brief Definition of TreeExtractor class that extracts ntuple information.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-28 created
 */

#pragma once

#include "EventDescriptor.h"

namespace analysis {

class TreeExtractor{
public:
    TreeExtractor(const std::string& inputFileName, bool extractMCtruth)
        :inputFile(new TFile(inputFileName.c_str(), "READ")), eventTree(*inputFile),
          electronTree(*inputFile),muonTree(*inputFile),tauTree(*inputFile),jetTree(*inputFile),
          vertexTree(*inputFile),genParticleTree(nullptr),current_entry(-1){
        if(inputFile->IsZombie())
            throw std::runtime_error("Input file not found.");

        if (extractMCtruth)
            genParticleTree = std::shared_ptr<ntuple::GenParticleTree>(new ntuple::GenParticleTree(*inputFile));
    }

    bool ExtractNext(EventDescriptor& descriptor)
    {
        descriptor.Clear();
        if (current_entry + 1 >= eventTree.GetEntries()) return false;
        ++current_entry;
        eventTree.GetEntry(current_entry);
        descriptor.eventInfo = eventTree.data;
        currentEventId = eventTree.data.EventId;

        if (genParticleTree)
            ReadObjects(*genParticleTree, descriptor.genParticles);

        ReadObjects(electronTree, descriptor.electrons);
        ReadObjects(muonTree, descriptor.muons);
        ReadObjects(tauTree, descriptor.taus);
        ReadObjects(jetTree, descriptor.jets);
        ReadObjects(vertexTree, descriptor.vertices);

        return true;
    }

private:
    template<typename Tree, typename Container>
    void ReadObjects(Tree& tree, Container& container)
    {
        if (tree.GetEntries() <= 0)
            return;

        for(Long64_t n = tree.GetReadEntry(); n < tree.GetEntries();) {
            if (currentEventId != tree.EventId()) break;
            container.push_back(tree.data);
            if(tree.GetEntry(++n) < 0)
                throw std::runtime_error("An I/O error while reading tree.");
        }
    }

private:
    std::shared_ptr<TFile> inputFile;
    ntuple::EventTree eventTree;
    ntuple::ElectronTree electronTree;
    ntuple::MuonTree muonTree;
    ntuple::TauTree tauTree;
    ntuple::JetTree jetTree;
    ntuple::VertexTree vertexTree;
    std::shared_ptr<ntuple::GenParticleTree> genParticleTree;
    Long64_t current_entry;
    UInt_t currentEventId;
};

} // analysis
