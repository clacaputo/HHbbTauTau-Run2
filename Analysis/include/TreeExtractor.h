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

        if (extractMCtruth){
            genParticleTree = new ntuple::GenParticleTree(*inputFile);
        }

        if (electronTree.GetEntries() > 0)
            electronTree.GetEntry(0);
        if (muonTree.GetEntries() > 0)
            muonTree.GetEntry(0);
        if (tauTree.GetEntries() > 0)
            tauTree.GetEntry(0);
        if (jetTree.GetEntries() > 0)
            jetTree.GetEntry(0);
        if (vertexTree.GetEntries() > 0)
            vertexTree.GetEntry(0);
        if (genParticleTree && genParticleTree->GetEntries() > 0)
            genParticleTree->GetEntry(0);
    }

    ~TreeExtractor() {
        delete genParticleTree;
        delete inputFile;
    }

    bool ExtractNext(EventDescriptor& descriptor)
    {
        descriptor.Clear();
        if (current_entry + 1 >= eventTree.GetEntries()) return false;
        ++current_entry;
        eventTree.GetEntry(current_entry);
        descriptor.eventInfo = ntuple::Event(eventTree);
        const ULong64_t currentEventId = eventTree.EventId();

        if (genParticleTree && genParticleTree->GetEntries() > 0){
            for(Long64_t n = genParticleTree->RootTree().GetReadEntry(); n < genParticleTree->GetEntries(); ++n)
            {
                genParticleTree->GetEntry(n);
                const ULong64_t eventId = genParticleTree->EventId();
                if (currentEventId != eventId) break;
                descriptor.genParticles.push_back(ntuple::GenParticle(*genParticleTree));
            }
        }

        if (electronTree.GetEntries() > 0){
            for(Long64_t n = electronTree.RootTree().GetReadEntry(); n < electronTree.GetEntries(); ++n)
            {
                electronTree.GetEntry(n);
                const ULong64_t eventId = electronTree.EventId();
                if (currentEventId != eventId) break;
                descriptor.electrons.push_back(ntuple::Electron(electronTree));
            }
        }

        if (muonTree.GetEntries() > 0){
            for(Long64_t n = muonTree.RootTree().GetReadEntry(); n < muonTree.GetEntries(); ++n)
            {
                muonTree.GetEntry(n);
                const ULong64_t eventId = muonTree.EventId();
                if (currentEventId != eventId) break;
                descriptor.muons.push_back(ntuple::Muon(muonTree));
            }
        }

        if (tauTree.GetEntries() > 0){
            for(Long64_t n = tauTree.RootTree().GetReadEntry(); n < tauTree.GetEntries(); ++n)
            {
                tauTree.GetEntry(n);
                const ULong64_t eventId = tauTree.EventId();
                if (currentEventId != eventId) break;
                descriptor.taus.push_back(ntuple::Tau(tauTree));
            }
        }

        if (jetTree.GetEntries() > 0){
            for(Long64_t n = jetTree.RootTree().GetReadEntry(); n < jetTree.GetEntries(); ++n)
            {
                jetTree.GetEntry(n);
                const ULong64_t eventId = jetTree.EventId();
                if (currentEventId != eventId) break;
                descriptor.jets.push_back(ntuple::Jet(jetTree));
            }
        }

        if (vertexTree.GetEntries() > 0){
            for(Long64_t n = vertexTree.RootTree().GetReadEntry(); n < vertexTree.GetEntries(); ++n)
            {
                vertexTree.GetEntry(n);
                const ULong64_t eventId = vertexTree.EventId();
                if (currentEventId != eventId) break;
                descriptor.vertices.push_back(ntuple::Vertex(vertexTree));
            }
        }

        return true;
    }
private:
    TFile* inputFile;
    ntuple::EventTree eventTree;
    ntuple::ElectronTree electronTree;
    ntuple::MuonTree muonTree;
    ntuple::TauTree tauTree;
    ntuple::JetTree jetTree;
    ntuple::VertexTree vertexTree;
    ntuple::GenParticleTree* genParticleTree;
    Long64_t current_entry;
};

}
