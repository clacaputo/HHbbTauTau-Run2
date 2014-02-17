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
#include "../include/CutTools.h"
#include "../include/GenParticle.h"

class SignalAnalyzerData : public root_ext::AnalyzerData {
public:
    SignalAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}
    TH1D_ENTRY(Pt_muon, 100, 0, 100)
    TH1D_ENTRY(Pt_tau, 100, 0, 100)
    TH1D_ENTRY(N_realtau, 100, 0, 100)
    TH1D_ENTRY_FIX(MuonSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(Counter, 1, 100, -0.5)
};

using ntuple::EventTree;

class HHbbTauTau_analyzer
{
public:
    HHbbTauTau_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0)
        : anaData(outputFileName), maxNumberOfEvents(_maxNumberOfEvents)
    {
        TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
        TTree* inputTree = dynamic_cast<TTree*> (inputFile->Get("treeCreator/vhtree"));
        event = new EventTree(inputTree);
        anaData.getOutputFile().cd();
        std::cout << "starting analyzer" << std::endl;
    }

    void Run()
    {
        if (event->fChain == nullptr) return;

        const Long64_t nentries = maxNumberOfEvents ? std::min(event->fChain->GetEntriesFast(),maxNumberOfEvents)
                                                    : event->fChain->GetEntriesFast();
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
           if (event->LoadTree(jentry) < 0)
               throw std::runtime_error("cannot read entry");
           event->fChain->GetEntry(jentry);
           ProcessEvent();
        }
    }

private:

    void ProcessEvent()
    {
        const analysis::GenEvent genEvent(*event);
        std::cout << "N gen particles = " << genEvent.particles.size() << std::endl;
        const IndexVector muons = CollectMuons();
        for (unsigned n = 0; n < muons.size(); ++n){
//            std::cout << "Muon pt = " << eventTree->Muon_pt[muons.at(n)] << std::endl;
        }
//        std::cout << std::endl;
    }

    IndexVector CollectMuons()
    {
        const auto muonSelector = [&](unsigned id) -> bool
            { return IsMuon(id); };

        const auto muonPtComparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Muon_pt[a] >  event->Muon_pt[b]; };

        return cuts::collect_objects(anaData.Counter(), anaData.MuonSelection(), event->nMuon,
                                     muonSelector, muonPtComparitor);
    }

    bool IsMuon(Int_t id)
    {
        using namespace cuts::muonID;
        int param_id = -1;
        const auto apply_cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.MuonSelection(), label); };

        try {
            apply_cut(true, ">0 mu cand");
            apply_cut(event->Muon_pt[id] > pt, "pt");
            apply_cut(std::abs(event->Muon_eta[id]) < eta, "eta");
            apply_cut(!isTrackerMuon || event->Muon_isTrackerMuon[id], "tracker");
            apply_cut(!isGlobalMuonPromptTight || event->Muon_isGlobalMuonPromptTight[id], "tight");
            apply_cut(!isPFMuon || event->Muon_isPFMuon[id], "PF");
            apply_cut(event->Muon_nChambers[id] > nChambers, "chamers");
            apply_cut(event->Muon_nMatchedStations[id] > nMatched_Stations, "stations");
            apply_cut(event->Muon_trackerLayersWithMeasurement[id] > trackerLayersWithMeasurement, "layers");
            apply_cut(event->Muon_pixHits[id] > pixHits, "pix_hits");
            apply_cut(event->Muon_globalChi2[id] < globalChiSquare, "chi2");
            apply_cut(std::abs(event->Muon_dB[id]) < dB, "dB");
            apply_cut(event->Muon_pfRelIso[id] < pFRelIso, "pFRelIso");
            return true;
        } catch(cuts::cut_failed&) {}
        return false;
    }

private:
    EventTree* event;
    SignalAnalyzerData anaData;
    Long64_t maxNumberOfEvents;
};
