/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#pragma once

#define EventTree_cxx
#include "../include/EventTree.h"
void EventTree::Loop(){} //just declared

#include "../include/AnalyzerData.h"
#include "../include/Particles.h"
#include "../include/Htautau_Summer13.h"
#include "../include/CutTools.h"
#include "../include/GenParticle.h"
#include "../include/MCfinalState.h"
#include "../include/Candidate.h"

#define SELECTION_ENTRY(name, n_bins, ...) \
    template<typename ...Args> \
    root_ext::SmartHistogram< TH1D >& name(const Args& ...args) { \
        auto event = &name##_event(args...); \
        auto eff_rel = &name##_effRel(args...); \
        auto eff_abs = &name##_effAbs(args...); \
        selectionDescriptors[event] = SelectionDescriptor(event, eff_rel, eff_abs, ##__VA_ARGS__); \
        return *event; } \
    TH1D_ENTRY_FIX(name##_event, 1, n_bins, -0.5) \
    TH1D_ENTRY_FIX(name##_effRel, 1, n_bins, -0.5) \
    TH1D_ENTRY_FIX(name##_effAbs, 1, n_bins, -0.5) \
    /**/

#define X(name, ...) \
    cuts::fill_histogram( event->name[id], _anaData.Get(&event->name[id], #name, ##__VA_ARGS__) )

namespace analysis {

class SignalAnalyzerData : public root_ext::AnalyzerData {
public:
    struct SelectionDescriptor {
        root_ext::SmartHistogram<TH1D> *event, *eff_relative, *eff_absolute;
        Int_t fix_bin;
        SelectionDescriptor(root_ext::SmartHistogram<TH1D> *_event = nullptr,
                            root_ext::SmartHistogram<TH1D> *eff_rel = nullptr,
                            root_ext::SmartHistogram<TH1D> *eff_abs = nullptr,
                            Int_t fix = std::numeric_limits<Int_t>::max())
            : event(_event), eff_relative(eff_rel), eff_absolute(eff_abs), fix_bin(fix) {}
    };

public:
    SignalAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}
    virtual ~SignalAnalyzerData()
    {
        Erase(Counter().Name());
        for(const auto& desc : selectionDescriptors) {
            cuts::fill_relative_selection_histogram(*desc.second.event, *desc.second.eff_relative, desc.second.fix_bin);
            cuts::fill_absolute_selection_histogram(*desc.second.event, *desc.second.eff_absolute);
        }
    }

    SELECTION_ENTRY(MuonSelection, 15)
    SELECTION_ENTRY(TauSelection, 15)
    SELECTION_ENTRY(ElectronSelection, 15)
    SELECTION_ENTRY(BJetSelection, 15)

    TH1D_ENTRY_FIX(Counter, 1, 100, -0.5)
    ENTRY_1D(float, Resonance_Mass)
    ENTRY_1D(float, Resonance_Pt)
    ENTRY_1D(float, Resonance_Eta)
    ENTRY_1D(float, Resonance_Phi)
    ENTRY_1D(float, Bjets_Pt_MC)
    ENTRY_1D(float, Higgs_leptonic_MC_Pt)
    ENTRY_1D(float, Higgs_BB_MC_Pt)
    ENTRY_2D(float, DR_bjets_vs_HiggsPt_MC)
    ENTRY_2D(float, DR_Higgs_vs_ResonancePt_MC)

protected:
    std::map<root_ext::SmartHistogram<TH1D>*, SelectionDescriptor> selectionDescriptors;
};

class BaseAnalyzer {
public:
    BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0, bool _useMCtruth = false)
        : outputFile(new TFile(outputFileName.c_str(),"RECREATE")), anaDataBeforeCut(*outputFile, "before_cut"),
          anaDataAfterCut(*outputFile, "after_cut"),
          maxNumberOfEvents(_maxNumberOfEvents), useMCtruth(_useMCtruth)
    {
        TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
        if(inputFile->IsZombie())
            throw std::runtime_error("Input file not found.");

        TTree* inputTree = dynamic_cast<TTree*> (inputFile->Get("treeCreator/vhtree"));
        event = boost::shared_ptr<EventTree>(new EventTree(inputTree,useMCtruth));
        std::cout << "Starting analyzer...\n";
    }

    virtual ~BaseAnalyzer() {}

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
    }

protected:
    virtual SignalAnalyzerData& GetAnaData() = 0;
    virtual void ProcessEvent() = 0;

    template<typename BaseSelector>
    CandidateVector CollectObjects(TH1D& selection_histogram, Int_t n_objects, const BaseSelector base_selector)
    {
        const auto selector = [&](unsigned id) -> analysis::Candidate
            { return base_selector(id, true, anaDataBeforeCut); };

        const auto selected = cuts::collect_objects(GetAnaData().Counter(), selection_histogram, n_objects, selector);
        for(const Candidate& candidate : selected)
            base_selector(candidate.index, false, anaDataAfterCut);
        return selected;
    }

    CandidateVector CollectMuons()
    {
        const auto base_selector = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData) -> Candidate
            { return SelectMuon(id, enabled, _anaData); };
        return CollectObjects(GetAnaData().MuonSelection(), event->nMuon, base_selector);
    }

    CandidateVector CollectTaus()
    {
        const auto base_selector = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData) -> Candidate
            { return SelectTau(id, enabled, _anaData); };
        return CollectObjects(GetAnaData().TauSelection(), event->nTau, base_selector);
    }

    CandidateVector CollectElectrons()
    {
        const auto base_selector = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData) -> Candidate
            { return SelectElectron(id, enabled, _anaData); };
        return CollectObjects(GetAnaData().ElectronSelection(), event->nElectron, base_selector);
    }

    CandidateVector CollectBJets(double csv, const std::string& selection_label)
    {
        const auto base_selector = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData) -> Candidate
            { return SelectBJet(id, csv, selection_label, enabled, _anaData); };
        return CollectObjects(GetAnaData().BJetSelection(selection_label), event->nJet, base_selector);
    }

    virtual Candidate SelectMuon(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData) = 0;
    virtual Candidate SelectTau(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData) = 0;
    virtual Candidate SelectElectron(Int_t id, bool enabled, root_ext::AnalyzerData& _anaData) = 0;

    Candidate SelectBJet(Int_t id, double csv, const std::string& selection_label, bool enabled,
                         root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::btag;
        cuts::Cutter cut(GetAnaData().Counter(), GetAnaData().BJetSelection(selection_label), enabled);

        cut(true, ">0 b-jet cand");
        cut(X(Jet_pt, selection_label) > pt, "pt");
        cut(std::abs( X(Jet_eta, selection_label) ) < eta, "eta");
        cut(X(Jet_combinedSecondaryVertexBTag, selection_label) > csv, "CSV");
        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(event->Jet_pt[id], event->Jet_eta[id], event->Jet_phi[id],
                              event->Jet_energy[id]);
        return analysis::Candidate(analysis::Candidate::Bjet, id, momentum);
    }

    template<typename Histogram>
    static CandidateVector FindCompatibleObjects(const CandidateVector& objects1, const CandidateVector& objects2,
                                                 double minDeltaR, Candidate::Type type, Histogram& mass)
    {
        CandidateVector result;
        for(const Candidate& object1 : objects1) {
            for(const Candidate& object2 : objects2) {
                const Candidate candidate(type, object1, object2);
                if(object2.momentum.DeltaR(object1.momentum) > minDeltaR) {
                    result.push_back(candidate);
                    mass.Fill(candidate.momentum.M());
                }
            }
        }
        return result;
    }

    void FindAnalysisFinalState(finalState::bbTauTau& final_state)
    {
        static const ParticleCodes resonanceCodes = { particles::Radion };
        static const ParticleCodes resonanceDecay = { particles::Higgs, particles::Higgs };
        static const ParticleCodes2D HiggsDecays = { { particles::b, particles::b },
                                                     { particles::tau, particles::tau } };

        genEvent = boost::shared_ptr<GenEvent>(new GenEvent(*event));

        const GenParticleSet resonances = genEvent->GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        GenParticleVector HiggsBosons;
        if(!FindDecayProducts(*final_state.resonance, resonanceDecay,HiggsBosons))
            throw std::runtime_error("Resonance does not decay into 2 Higgs");

        GenParticleVector2D HiggsDecayProducts;
        GenParticleIndexVector HiggsIndexes;
        if(!FindDecayProducts2D(HiggsBosons,HiggsDecays,HiggsDecayProducts,HiggsIndexes))
            throw std::runtime_error("NOT HH -> bb tautau");

        final_state.b_jets = HiggsDecayProducts.at(0);
        final_state.taus = HiggsDecayProducts.at(1);

        final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(1));
        final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(0));

        GetAnaData().Resonance_Mass().Fill(final_state.resonance->momentum.M());
        GetAnaData().Resonance_Pt().Fill(final_state.resonance->momentum.Pt());
        GetAnaData().Resonance_Eta().Fill(final_state.resonance->momentum.Eta());
        GetAnaData().Resonance_Phi().Fill(final_state.resonance->momentum.Phi());

        GetAnaData().Bjets_Pt_MC().Fill(final_state.b_jets.at(0)->momentum.Pt());
        GetAnaData().Bjets_Pt_MC().Fill(final_state.b_jets.at(1)->momentum.Pt());
        GetAnaData().Higgs_leptonic_MC_Pt().Fill(final_state.Higgs_TauTau->momentum.Pt());
        GetAnaData().Higgs_BB_MC_Pt().Fill(final_state.Higgs_BB->momentum.Pt());
        GetAnaData().DR_bjets_vs_HiggsPt_MC().Fill(final_state.Higgs_BB->momentum.Pt(),
                                         final_state.b_jets.at(0)->momentum.DeltaR(final_state.b_jets.at(1)->momentum));
        GetAnaData().DR_Higgs_vs_ResonancePt_MC().Fill(final_state.resonance->momentum.Pt(),
                                         final_state.Higgs_TauTau->momentum.DeltaR(final_state.Higgs_BB->momentum));
    }

protected:
    boost::shared_ptr<EventTree> event;
    boost::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut;
    Long64_t maxNumberOfEvents;
    bool useMCtruth;
    boost::shared_ptr<GenEvent> genEvent;
};

} // analysis
