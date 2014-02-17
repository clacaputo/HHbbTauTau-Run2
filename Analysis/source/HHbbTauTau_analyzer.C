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
    TH1D_ENTRY_FIX(EventSelection, 1, 1, -0.5)
    TH1D_ENTRY_FIX(MuonSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(TauSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(ElectronSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(BJetSelection, 1, 15, -0.5)
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
        //        const analysis::GenEvent genEvent(*event);
        //        std::cout << "N gen particles = " << genEvent.particles.size() << std::endl;

        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.EventSelection(), ++param_id, anaData.EventSelection(), label); };

        try {
            cut(true, "total");
            const IndexVector muons = CollectMuons();
            cut(muons.size(), "muon");
            const IndexVector taus = CollectTaus();
            cut(taus.size(), "tau");
            const IndexVector electrons = CollectElectrons();
            cut(!electrons.size(), "no_electron");
            const IndexVector b_jets = CollectBJets();
            cut(b_jets.size() == 2, "2_b_jets");
        } catch(cuts::cut_failed&) {}
    }

    IndexVector CollectMuons()
    {
        const auto selector = [&](unsigned id) -> bool
            { return IsMuon(id); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Muon_pt[a] >  event->Muon_pt[b]; };

        return cuts::collect_objects(anaData.Counter(), anaData.MuonSelection(), event->nMuon, selector, comparitor);
    }

    bool IsMuon(Int_t id)
    {
        using namespace cuts::muonID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.MuonSelection(), label); };

        try {
            cut(true, ">0 mu cand");
            cut(event->Muon_pt[id] > pt, "pt");
            cut(std::abs(event->Muon_eta[id]) < eta, "eta");
            cut(!isTrackerMuon || event->Muon_isTrackerMuon[id], "tracker");
            cut(!isGlobalMuonPromptTight || event->Muon_isGlobalMuonPromptTight[id], "tight");
            cut(!isPFMuon || event->Muon_isPFMuon[id], "PF");
            cut(event->Muon_nChambers[id] > nChambers, "chamers");
            cut(event->Muon_nMatchedStations[id] > nMatched_Stations, "stations");
            cut(event->Muon_trackerLayersWithMeasurement[id] > trackerLayersWithMeasurement, "layers");
            cut(event->Muon_pixHits[id] > pixHits, "pix_hits");
            cut(event->Muon_globalChi2[id] < globalChiSquare, "chi2");
            cut(std::abs(event->Muon_dB[id]) < dB, "dB");
            cut(event->Muon_pfRelIso[id] < pFRelIso, "pFRelIso");
            return true;
        } catch(cuts::cut_failed&) {}
        return false;
    }

    IndexVector CollectTaus()
    {
        const auto selector = [&](unsigned id) -> bool
            { return IsTau(id); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Tau_pt[a] >  event->Tau_pt[b]; };

        return cuts::collect_objects(anaData.Counter(), anaData.TauSelection(), event->nTau, selector, comparitor);
    }

    bool IsTau(Int_t id)
    {
        using namespace cuts::tauID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.TauSelection(), label); };

        try {
            cut(true, ">0 tau cand");
            cut(event->Tau_pt[id] > pt, "pt");
            cut(std::abs(event->Tau_eta[id]) < eta, "eta");
            cut(event->Tau_decayModeFinding[id] > decayModeFinding, "decay_mode");
//            cut(event->Tau_byLooseIsolationDeltaBetaCorr[id] > byLooseIsolationDeltaBetaCorr, "loose_beta");
            cut(event->Tau_againstMuonTight[id] > againstMuonTight, "vs_mu_tight");
            cut(event->Tau_againstElectronLoose[id] > againstElectronLoose, "vs_e_loose");
            return true;
        } catch(cuts::cut_failed&) {}
        return false;
    }

    IndexVector CollectElectrons()
    {
        const auto selector = [&](unsigned id) -> bool
            { return IsElectron(id); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Electron_pt[a] >  event->Electron_pt[b]; };

        return cuts::collect_objects(anaData.Counter(), anaData.ElectronSelection(),
                                     event->nElectron, selector, comparitor);
    }

    bool IsElectron(Int_t id)
    {
        using namespace cuts::electronID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.ElectronSelection(), label); };

        try {
            cut(true, ">0 ele cand");
            cut(event->Electron_pt[id] > pt, "pt");
            const double eta = std::abs(event->Electron_eta[id]);
            cut(eta < eta_high && (eta < eta_CrackVeto_low || eta > eta_CrackVeto_high), "eta");
            // cut dz_pv
            cut(event->Electron_missingHits[id] < missingHits, "mis_hits");
            cut(event->Electron_hasMatchedConv[id] < hasMatchedConv, "has_conv");
            cut(event->Electron_dB[id] < dB, "dB");
            const size_t pt_index = event->Electron_pt[id] < ref_pt ? 0 : 1;
            const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
            cut(event->Electron_mvaPOGNonTrig[id] > MVApogNonTrig[pt_index][eta_index], "mva");
            return true;
        } catch(cuts::cut_failed&) {}
        return false;
    }

    IndexVector CollectBJets()
    {
        const auto selector = [&](unsigned id) -> bool
            { return IsBJet(id); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Jet_pt[a] >  event->Jet_pt[b]; };

        return cuts::collect_objects(anaData.Counter(), anaData.BJetSelection(), event->nJet, selector, comparitor);
    }

    bool IsBJet(Int_t id)
    {
        using namespace cuts::btag;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.BJetSelection(), label); };

        try {
            cut(true, ">0 b-jet cand");
            cut(event->Jet_pt[id] > pt, "pt");
            cut(std::abs(event->Jet_eta[id]) < eta, "eta");
            cut(event->Jet_combinedSecondaryVertexBTag[id] > CSV, "CSV");
            return true;
        } catch(cuts::cut_failed&) {}
        return false;
    }

private:
    EventTree* event;
    SignalAnalyzerData anaData;
    Long64_t maxNumberOfEvents;
};
