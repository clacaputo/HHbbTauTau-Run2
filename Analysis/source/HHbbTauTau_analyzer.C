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
#include "../include/Cuts.h"

class SignalAnalyzerData : public root_ext::AnalyzerData {
public:
    SignalAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}
    TH1D_ENTRY(Pt_muon, 100, 0, 100)
    TH1D_ENTRY(Pt_tau, 100, 0, 100)
    TH1D_ENTRY(N_realtau, 100, 0, 100)
    TH1D_ENTRY_FIX(MuonSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(MuonCounter, 1, 15, -0.5)
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
        eventTree = new EventTree(inputTree);
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
        const IndexVector muons = CollectMuons();
        for (unsigned n = 0; n < muons.size(); ++n){
//            std::cout << "Muon pt = " << eventTree->Muon_pt[muons.at(n)] << std::endl;
        }
//        std::cout << std::endl;
    }

    IndexVector CollectMuons()
    {
        using namespace cuts::muonID;
        IndexVector selected;
        anaData.MuonCounter().Reset();
        for (Int_t n = 0; n < eventTree->nMuon; ++n){
            try {
                    int param_id = -1;
                    const auto apply_cut = [&](bool expected, const std::string& label) {
                        cuts::apply_cut(expected, anaData.MuonCounter(), ++param_id, anaData.MuonSelection(), label);
                    };

                    apply_cut(true, "total");
                    apply_cut(eventTree->Muon_pt[n] > pt, "pt");
                    apply_cut(std::abs(eventTree->Muon_eta[n]) < eta, "eta");
                    apply_cut(!isTrackerMuon || eventTree->Muon_isTrackerMuon[n], "tracker");
                    apply_cut(!isGlobalMuonPromptTight || eventTree->Muon_isGlobalMuonPromptTight[n], "tight");
                    apply_cut(!isPFMuon || eventTree->Muon_isPFMuon[n], "PF");
                    apply_cut(eventTree->Muon_nChambers[n] > nChambers, "chamers");
                    apply_cut(eventTree->Muon_nMatchedStations[n] > nMatched_Stations, "stations");
                    apply_cut(eventTree->Muon_trackerLayersWithMeasurement[n] > trackerLayersWithMeasurement, "layers");
                    apply_cut(eventTree->Muon_pixHits[n] > pixHits, "pix_hits");
                    apply_cut(eventTree->Muon_globalChi2[n] < globalChiSquare, "chi2");
                    apply_cut(std::abs(eventTree->Muon_dB[n]) < dB, "dB");
                    apply_cut(eventTree->Muon_pfRelIso[n] < pFRelIso, "pFRelIso");
                    selected.push_back(n);
            } catch(cuts::cut_failed&) {}
        }

        cuts::fill_selection_histogram(anaData.MuonSelection(), anaData.MuonCounter());

        const auto muonPtComparitor = [&](unsigned a, unsigned b) -> bool
            { return eventTree->Muon_pt[a] >  eventTree->Muon_pt[b]; };

        std::sort(selected.begin(), selected.end(), muonPtComparitor);

        return selected;
    }

private:
    EventTree* eventTree;
    SignalAnalyzerData anaData;
    Long64_t maxNumberOfEvents;
};
