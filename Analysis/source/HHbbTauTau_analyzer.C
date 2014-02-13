/*
 * \file HHbbTauTau_analyzer.C
 * \brief Analysis HHbbTauTau.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-12 created
 */

#include <vector>
#include <boost/shared_ptr.hpp>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TObject.h>

namespace ntuple {
    using namespace std;

    #define EventTree_cxx
    #include "../include/EventTree.h"

    void EventTree::Loop(){} //just declared
}

#include "../include/AnalyzerData.h"
#include "../include/Particles.h"
#include "../include/Htautau_Summer13.h"

class SignalAnalyzerData : public root_ext::AnalyzerData {
public:
    SignalAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}
    TH1D_ENTRY(Pt_muon, 100, 0, 100)
    TH1D_ENTRY(Pt_tau, 100, 0, 100)
    TH1D_ENTRY(N_realtau, 100, 0, 100)

};

using ntuple::EventTree;

class HHbbTauTau_analyzer
{
public:

    typedef std::vector<Int_t> IndexVector;

    HHbbTauTau_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0)
        : anaData(outputFileName), maxNumberOfEvents(_maxNumberOfEvents)
    {
        TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
        TTree* inputTree = dynamic_cast<TTree*> (inputFile->Get("treeCreator/vhtree"));
        eventTree = boost::shared_ptr<EventTree>(new EventTree(inputTree));
        anaData.getOutputFile().cd();
        std::cout << "starting analyzer" << std::endl;
    }

    void Run()
    {
        if (eventTree->fChain == nullptr) return;

        const Long64_t nentries = maxNumberOfEvents ? std::min(eventTree->fChain->GetEntriesFast(),maxNumberOfEvents)
                                                    : eventTree->fChain->GetEntriesFast();
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
           if (eventTree->LoadTree(jentry) < 0)
               throw std::runtime_error("cannot read entry");
           eventTree->fChain->GetEntry(jentry);
           ProcessEvent();
        }
    }

private:

    void ProcessEvent()
    {
        const IndexVector muons = CollectMuon();
        for (unsigned n = 0; n < muons.size(); ++n){
            std::cout << "Muon pt = " << eventTree->Muon_pt[muons.at(n)] << std::endl;
        }
    }

    IndexVector CollectMuon()
    {
        using namespace cuts::muonID;
        IndexVector selected;
        for (Int_t n = 0; n < eventTree->nMuon; ++n){
            if (eventTree->Muon_pt[n] <= pt) continue;
            if (std::abs(eventTree->Muon_eta[n]) >= eta) continue;
            if (isTrackerMuon && !eventTree->Muon_isTrackerMuon[n]) continue;
            if (isGlobalMuonPromptTight && !eventTree->Muon_isGlobalMuonPromptTight[n]) continue;
            if (isPFMuon && !eventTree->Muon_isPFMuon[n]) continue;
            if (eventTree->Muon_nChambers[n] <= nChambers) continue;
            if (eventTree->Muon_nMatchedStations[n] <= nMatched_Stations) continue;
            if (eventTree->Muon_trackerLayersWithMeasurement[n] <= trackerLayersWithMeasurement) continue;
            if (eventTree->Muon_pixHits[n] <= pixHits) continue;
            if (eventTree->Muon_globalChi2[n] >= globalChiSquare) continue;
            if (std::abs(eventTree->Muon_dB[n]) >= dB) continue;
            if (eventTree->Muon_pfRelIso[n] >= pFRelIso) continue;
            selected.push_back(n);
        }

        auto MuonComparitor = [&](unsigned a, unsigned b) -> bool {
            return eventTree->Muon_pt[a] >  eventTree->Muon_pt[b];
        };

        std::sort(selected.begin(), selected.end(), MuonComparitor);

        return selected;
    }

private:
    boost::shared_ptr<EventTree> eventTree;
    SignalAnalyzerData anaData;
    Long64_t maxNumberOfEvents;
};
