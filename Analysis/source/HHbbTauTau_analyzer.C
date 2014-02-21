/*
 * \file HHbbTauTau_analyzer.C
 * \brief Analysis HHbbTauTau.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-12 created
 */

#include <vector>
#include <type_traits>
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

    TH1D_ENTRY_FIX(EventSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(MuonSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(TauSelection, 1, 15, -0.5)
    TH1D_ENTRY_FIX(ElectronSelection, 1, 15, -0.5)
    TH1D_MAP_ENTRY_FIX(BJetSelection, 1, 15, -0.5)

    TH1D_ENTRY_FIX(EventSelection_relative, 1, 15, -0.5)
    TH1D_ENTRY_FIX(MuonSelection_relative, 1, 15, -0.5)
    TH1D_ENTRY_FIX(TauSelection_relative, 1, 15, -0.5)
    TH1D_ENTRY_FIX(ElectronSelection_relative, 1, 15, -0.5)
    TH1D_MAP_ENTRY_FIX(BJetSelection_relative, 1, 15, -0.5)

    TH1D_ENTRY_FIX(Counter, 1, 100, -0.5)
};

using ntuple::EventTree;

namespace cuts {
namespace HHbbTauTau {

const unsigned N_bjet_Loose = 2;
const unsigned N_bjet_Medium = 1;
}
}

#define X(name) \
    cuts::fill_histogram( event->name[id], \
    _anaData.getSmart< std::remove_reference< decltype(event->name[id]) >::type >(#name) )


#define XX(name, suffix) \
    cuts::fill_histogram( event->name[id], \
    _anaData.getSmart< std::remove_reference< decltype(event->name[id]) >::type >(#name, suffix) )


class HHbbTauTau_analyzer
{
public:
    HHbbTauTau_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0)
        : anaData(outputFileName), anaDataBeforeCut(anaData.getOutputFile(), "before_cut"),
          anaDataAfterCut(anaData.getOutputFile(), "after_cut"),
          maxNumberOfEvents(_maxNumberOfEvents)
    {
        TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
        if(inputFile->IsZombie())
            throw std::runtime_error("Input file not found.");

        TTree* inputTree = dynamic_cast<TTree*> (inputFile->Get("treeCreator/vhtree"));
        event = new EventTree(inputTree);
        anaData.getOutputFile().cd();
        std::cout << "Starting analyzer...\n";
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
            try {
                ProcessEvent();
            } catch(cuts::cut_failed&) {}
        }

        FillRelativeSelectionHistograms();
    }

private:

    void ProcessEvent()
    {
//        const analysis::GenEvent genEvent(*event);
//        std::cout << "N gen particles = " << genEvent.particles.size() << std::endl;

        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { cuts::apply_cut(expected, anaData.EventSelection(), ++param_id, anaData.EventSelection(), label); };

        const auto try_cut = [&](bool expected, const std::string& label)
            { try { cut(expected, label); } catch(cuts::cut_failed&) {} };

        cut(true, "total");
        const IndexVector muons = CollectMuons();
        cut(muons.size(), "muon");
        const IndexVector taus = CollectTaus();
        cut(taus.size(), "tau");
        const IndexPairVector mu_tau_pairs = FindCompatibleLeptonCombinations(muons, taus);
        cut(mu_tau_pairs.size(), "mu_tau");
        const IndexVector electrons = CollectElectrons();
        cut(!electrons.size(), "no_electron");
        const IndexVector b_jets_loose = CollectBJets(cuts::btag::CSVL, "loose");
        const IndexVector b_jets_medium = CollectBJets(cuts::btag::CSVM, "medium");
        try_cut(b_jets_loose.size() == 2, "2b_loose");
        try_cut(b_jets_loose.size() == 2 && b_jets_medium.size() >= 1, "1b_loose+1b_medium");
        try_cut(b_jets_medium.size() == 2, "2b_medium");
        try_cut(b_jets_loose.size() >= 2, ">=2b_loose");
        try_cut(b_jets_loose.size() >= 2 && b_jets_medium.size() >= 1, ">=1b_loose+>=1b_medium");
        try_cut(b_jets_medium.size() >= 2, ">=2b_medium");
    }

    IndexVector CollectMuons()
    {
        const auto selector = [&](unsigned id) { SelectMuon(id, true, anaDataBeforeCut); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Muon_pt[a] >  event->Muon_pt[b]; };

        const auto selected = cuts::collect_objects(anaData.Counter(), anaData.MuonSelection(), event->nMuon, selector,
                                                    comparitor);
        for(Int_t id : selected) SelectMuon(id, false, anaDataAfterCut);
        return selected;
    }

    void SelectMuon(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::muonID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { if(apply_cut) cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.MuonSelection(), label); };

        cut(true, ">0 mu cand");
        cut(X(Muon_pt) > pt, "pt");
        cut(std::abs( X(Muon_eta) ) < eta, "eta");
        cut(!isTrackerMuon || X(Muon_isTrackerMuon), "tracker");
        cut(!isGlobalMuonPromptTight || X(Muon_isGlobalMuonPromptTight), "tight");
        cut(!isPFMuon || X(Muon_isPFMuon), "PF");
        cut(X(Muon_nChambers) > nChambers, "chamers");
        cut(X(Muon_nMatchedStations) > nMatched_Stations, "stations");
        cut(X(Muon_trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(Muon_pixHits) > pixHits, "pix_hits");
        cut(X(Muon_globalChi2) < globalChiSquare, "chi2");
        cut(std::abs( X(Muon_dB) ) < dB, "dB");
        cut(X(Muon_pfRelIso) < pFRelIso, "pFRelIso");
    }

    IndexVector CollectTaus()
    {
        const auto selector = [&](unsigned id) { SelectTau(id, true, anaDataBeforeCut); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Tau_pt[a] >  event->Tau_pt[b]; };

        const auto selected = cuts::collect_objects(anaData.Counter(), anaData.TauSelection(), event->nTau, selector,
                                                    comparitor);
        for(Int_t id : selected) SelectTau(id, false, anaDataAfterCut);
        return selected;
    }

    void SelectTau(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::tauID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { if(apply_cut) cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.TauSelection(), label); };

        cut(true, ">0 tau cand");
        cut(X(Tau_pt) > pt, "pt");
        cut(std::abs( X(Tau_eta) ) < eta, "eta");
        cut(X(Tau_decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits) > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(Tau_againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(Tau_againstElectronLoose) > againstElectronLoose, "vs_e_loose");
    }

    IndexVector CollectElectrons()
    {
        const auto selector = [&](unsigned id) { SelectElectron(id, true, anaDataBeforeCut); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Electron_pt[a] >  event->Electron_pt[b]; };

        const auto selected = cuts::collect_objects(anaData.Counter(), anaData.ElectronSelection(),
                                     event->nElectron, selector, comparitor);
        for(Int_t id : selected) SelectElectron(id, false, anaDataAfterCut);
        return selected;

    }

    void SelectElectron(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::electronID;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { if(apply_cut)
                cuts::apply_cut(expected, anaData.Counter(), ++param_id, anaData.ElectronSelection(), label); };

        cut(true, ">0 ele cand");
        cut(X(Electron_pt) > pt, "pt");
        const double eta = std::abs( X(Electron_eta) );
        cut(eta < eta_high && (eta < eta_CrackVeto_low || eta > eta_CrackVeto_high), "eta");
        // cut dz_pv
        cut(X(Electron_missingHits) < missingHits, "mis_hits");
        cut(X(Electron_hasMatchedConv) < hasMatchedConv, "has_conv");
        cut(X(Electron_dB) < dB, "dB");
        const size_t pt_index = event->Electron_pt[id] < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(Electron_mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
    }

    IndexVector CollectBJets(double csv, const std::string& selection_label)
    {
        const auto selector = [&](unsigned id) { SelectBJet(id, csv, selection_label, true, anaDataBeforeCut); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return event->Jet_pt[a] >  event->Jet_pt[b]; };

        const auto selected = cuts::collect_objects(anaData.Counter(), anaData.BJetSelection(selection_label),
                                     event->nJet, selector, comparitor);
        for(Int_t id : selected) SelectBJet(id, csv, selection_label, false, anaDataAfterCut);
        return selected;
    }

    void SelectBJet(Int_t id, double csv, const std::string& selection_label, bool apply_cut,
                    root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::btag;
        int param_id = -1;
        const auto cut = [&](bool expected, const std::string& label)
            { if(apply_cut) cuts::apply_cut(expected, anaData.Counter(), ++param_id,
                              anaData.BJetSelection(selection_label), label); };

        cut(true, ">0 b-jet cand");
        cut(XX(Jet_pt, selection_label) > pt, "pt");
        cut(std::abs( XX(Jet_eta, selection_label) ) < eta, "eta");
        cut(XX(Jet_combinedSecondaryVertexBTag, selection_label) > csv, "CSV");
    }

    IndexPairVector FindCompatibleLeptonCombinations(const IndexVector& muons, const IndexVector& taus)
    {
        IndexPairVector result;
        for(auto mu_id : muons) {
            TVector3 mu_direction;
            mu_direction.SetPtEtaPhi(event->Muon_pt[mu_id], event->Muon_eta[mu_id], event->Muon_phi[mu_id]);
            for(auto tau_id : taus) {
                TVector3 tau_direction;
                tau_direction.SetPtEtaPhi(event->Tau_pt[tau_id], event->Tau_eta[tau_id], event->Tau_phi[tau_id]);
                if(tau_direction.DeltaR(mu_direction) > cuts::DeltaR_signalLeptons)
                    result.push_back(IndexPair(mu_id, tau_id));
            }
        }
        return result;
    }

    void FillRelativeSelectionHistograms()
    {
        cuts::fill_relative_selection_histogram(anaData.EventSelection(), anaData.EventSelection_relative(), 5);
        cuts::fill_relative_selection_histogram(anaData.MuonSelection(), anaData.MuonSelection_relative());
        cuts::fill_relative_selection_histogram(anaData.TauSelection(), anaData.TauSelection_relative());
        cuts::fill_relative_selection_histogram(anaData.ElectronSelection(), anaData.ElectronSelection_relative());
        cuts::fill_relative_selection_histogram(anaData.BJetSelection("loose"), anaData.BJetSelection_relative("loose"));
        cuts::fill_relative_selection_histogram(anaData.BJetSelection("medium"), anaData.BJetSelection_relative("medium"));
    }

private:
    EventTree* event;
    SignalAnalyzerData anaData;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut;
    Long64_t maxNumberOfEvents;
};
