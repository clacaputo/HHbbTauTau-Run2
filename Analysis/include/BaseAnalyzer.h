/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#pragma once
#include "TreeExtractor.h"
#include "../include/AnalyzerData.h"
#include "../include/Particles.h"
#include "../include/Htautau_Summer13.h"
#include "../include/CutTools.h"
#include "../include/GenParticle.h"
#include "../include/MCfinalState.h"
#include "../include/Candidate.h"
#include <chrono>
#include <iomanip>
#include <functional>

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

#define X(name) \
    cuts::fill_histogram( object.name, _anaData.Get(&object.name, #name, selection_label) )

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

class Timer {
public:
    typedef std::chrono::high_resolution_clock clock;
    Timer(unsigned _report_interval)
        : start(clock::now()), block_start(start), report_interval(_report_interval) {}

    void Report(size_t event_id, bool final_report = false)
    {
        using namespace std::chrono;
        const auto now = clock::now();
        const auto since_last_report = duration_cast<seconds>(now - block_start).count();
        if(!final_report && since_last_report < report_interval) return;

        const auto since_start = duration_cast<seconds>(now - start).count();
        const double speed = ((double) event_id) / since_start;
        if(final_report)
            std::cout << "Total: ";
        std::cout << "time = " << since_start << " seconds, events processed = " << event_id
                  << ", average speed = " << std::setprecision(1) << std::fixed << speed << " events/s\n";
        block_start = now;
    }

private:
    clock::time_point start, block_start;
    unsigned report_interval;
};

class BaseAnalyzer {
public:
    BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                        size_t _maxNumberOfEvents = 0, bool _useMCtruth = false)
        : treeExtractor(inputFileName, _useMCtruth), outputFile(new TFile(outputFileName.c_str(),"RECREATE")),
          anaDataBeforeCut(*outputFile, "before_cut"),anaDataAfterCut(*outputFile, "after_cut"),
          maxNumberOfEvents(_maxNumberOfEvents), useMCtruth(_useMCtruth)
    {
        std::cout << "Starting analyzer...\n";
    }

    virtual ~BaseAnalyzer() {}

    void Run()
    {
        Timer timer(10);
        size_t n = 0;
        for(; !maxNumberOfEvents || n < maxNumberOfEvents; ++n) {
            if(!treeExtractor.ExtractNext(event))
                break;
            timer.Report(n);
            try {
                ProcessEvent();
            } catch(cuts::cut_failed&) {}
        }
        timer.Report(n, true);
    }

protected:
    typedef std::function< Candidate (size_t, bool, root_ext::AnalyzerData&) > BaseSelector;

    virtual SignalAnalyzerData& GetAnaData() = 0;
    virtual void ProcessEvent() = 0;

    CandidateVector CollectObjects(TH1D& selection_histogram, size_t n_objects, const BaseSelector& base_selector)
    {
        const auto selector = [&](size_t id) -> analysis::Candidate
            { return base_selector(id, true, anaDataBeforeCut); };

        const auto selected = cuts::collect_objects(GetAnaData().Counter(), selection_histogram, n_objects, selector);
        for(const Candidate& candidate : selected)
            base_selector(candidate.index, false, anaDataAfterCut);
        return selected;
    }

    template<typename BaseSelectorMethod>
    CandidateVector CollectObjects(TH1D& selection_histogram, size_t n_objects, bool signal,
                                   BaseSelectorMethod signal_selector_method,
                                   BaseSelectorMethod bkg_selector_method)
    {
        const BaseSelector base_selector_signal = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return (this->*signal_selector_method)(id, enabled, _anaData); };
        const BaseSelector base_selector_bkg = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return (this->*bkg_selector_method)(id, enabled, _anaData); };
        const auto base_selector = signal ? base_selector_signal : base_selector_bkg;

        return CollectObjects(selection_histogram, n_objects, base_selector);

    }

    CandidateVector CollectMuons(bool signal = true)
    {
        return CollectObjects(GetAnaData().MuonSelection(), event.muons().size(), signal,
                              &BaseAnalyzer::SelectMuon, &BaseAnalyzer::SelectBackgroundMuon);
    }

    CandidateVector CollectTaus(bool signal = true)
    {
        return CollectObjects(GetAnaData().TauSelection(), event.taus().size(), signal,
                              &BaseAnalyzer::SelectTau, &BaseAnalyzer::SelectBackgroundTau);
    }

    CandidateVector CollectElectrons(bool signal = true)
    {
        return CollectObjects(GetAnaData().ElectronSelection(), event.electrons().size(), signal,
                              &BaseAnalyzer::SelectElectron, &BaseAnalyzer::SelectBackgroundElectron);
    }

    CandidateVector CollectBJets(double csv, const std::string& selection_label, bool signal = true)
    {
        const BaseSelector base_selector_signal = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBJet(id, enabled, _anaData, csv, selection_label); };
        const BaseSelector base_selector_bkg = [&](unsigned id, bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBackgroundBJet(id, enabled, _anaData, csv, selection_label); };
        const auto base_selector = signal ? base_selector_signal : base_selector_bkg;

        return CollectObjects(GetAnaData().BJetSelection(selection_label), event.jets().size(), base_selector);
    }

    virtual Candidate SelectMuon(size_t id, bool enabled, root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Muon selection for signal not implemented");
    }
    virtual Candidate SelectTau(size_t id, bool enabled, root_ext::AnalyzerData& _anaData) {
        throw std::runtime_error("Tau selection for signal not implemented");
    }
    virtual Candidate SelectElectron(size_t id, bool enabled, root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Electron selection for signal not implemented");
    }
    virtual Candidate SelectBackgroundMuon(size_t id, bool enabled, root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Muon selection for background not implemented");
    }
    virtual Candidate SelectBackgroundTau(size_t id, bool enabled, root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Tau selection for background not implemented");
    }
    virtual Candidate SelectBackgroundElectron(size_t id, bool enabled, root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Electron selection for background not implemented");
    }

    Candidate SelectBJet(size_t id, bool enabled, root_ext::AnalyzerData& _anaData,
                         double csv, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::btag::signal;
        cuts::Cutter cut(GetAnaData().Counter(), GetAnaData().BJetSelection(selection_label), enabled);

        const ntuple::Jet& object = event.jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags) > CSV, "CSV");

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
    }

    Candidate SelectBackgroundBJet(size_t id, bool enabled, root_ext::AnalyzerData& _anaData,
                                   double csv, const std::string& _selection_label)
    {
        using namespace cuts::Htautau_Summer13::btag::veto;
        const std::string selection_label = _selection_label + "_bkg";
        cuts::Cutter cut(GetAnaData().Counter(), GetAnaData().BJetSelection(selection_label), enabled);

        const ntuple::Jet& object = event.jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags) > CSV, "CSV");
        cut(X(passLooseID) == passLooseID, "passLooseID");
        //DeltaR with leptons is missing

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
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

        genEvent.Initialize(event.genParticles());

        const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        GenParticlePtrVector HiggsBosons;
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
    EventDescriptor event;
    TreeExtractor treeExtractor;
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut;
    size_t maxNumberOfEvents;
    bool useMCtruth;
    GenEvent genEvent;
};

} // analysis
