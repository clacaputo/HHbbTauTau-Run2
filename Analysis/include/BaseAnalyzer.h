/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#pragma once
#include "TreeExtractor.h"
#include "AnalyzerData.h"
#include "Particles.h"
#include "Htautau_Summer13.h"
#include "CutTools.h"
#include "GenParticle.h"
#include "MCfinalState.h"
#include "Candidate.h"
#include "custom_cuts.h"
#include "SVfit.h"
#include <chrono>
#include <iomanip>
#include <functional>

#include "RunReport.h"

#define SELECTION_ENTRY(name) \
    ENTRY_1D(cuts::ObjectSelector, name) \
    /**/

/*
#define X(name, ...) \
    cuts::fill_histogram( object.name, _anaData.Get(&object.name, #name, selection_label), weight)

#define Y(name, ...) \
    cuts::fill_histogram( name, _anaData.Get(&name, #name, selection_label), weight)
*/

#define X(name, ...) \
    cuts::fill_histogram( object.name, _anaData.Get((TH1D*)nullptr, #name, selection_label, ##__VA_ARGS__), weight)

#define Y(name, ...) \
    cuts::fill_histogram( name, _anaData.Get((TH1D*)nullptr, #name, selection_label, ##__VA_ARGS__), weight)

namespace analysis {

class SignalAnalyzerData : public root_ext::AnalyzerData {

public:
    SignalAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}

    SELECTION_ENTRY(EventSelection)
    SELECTION_ENTRY(VertexSelection)
    SELECTION_ENTRY(MuonSelection)
    SELECTION_ENTRY(TauSelection)
    SELECTION_ENTRY(ElectronSelection)
    SELECTION_ENTRY(BJetSelection)
    SELECTION_ENTRY(MuonSelectionBkg)
    SELECTION_ENTRY(TauSelectionBkg)
    SELECTION_ENTRY(ElectronSelectionBkg)
    SELECTION_ENTRY(BJetSelectionBkg)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)

    ENTRY_1D(float, Resonance_Mass)
    ENTRY_1D(float, Resonance_Pt)
    ENTRY_1D(float, Resonance_Eta)
    ENTRY_1D(float, Resonance_Phi)
    ENTRY_1D(float, Bjets_Pt_MC)
    ENTRY_1D(float, Higgs_leptonic_MC_Pt)
    ENTRY_1D(float, Higgs_BB_MC_Pt)
    ENTRY_2D(float, DR_bjets_vs_HiggsPt_MC)
    ENTRY_2D(float, DR_Higgs_vs_ResonancePt_MC)


};

class Timer {
public:
    typedef std::chrono::high_resolution_clock clock;
    Timer(unsigned _report_interval)
        : start(clock::now()), block_start(start), report_interval(_report_interval)
    {
        std::cout << "Starting analyzer...\n";
    }

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
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0, bool _useMCtruth = false,
                 const std::string& reweightFileName = "none")
        : timer(10), treeExtractor(_prefix == "none" ? "" : _prefix, inputFileName, _useMCtruth),
          outputFile(new TFile(outputFileName.c_str(),"RECREATE")),
          anaDataBeforeCut(*outputFile, "before_cut"), anaDataAfterCut(*outputFile, "after_cut"),
          anaDataFinalSelection(*outputFile, "final_selection"), runReport(outputFileName + ".txt"),
          maxNumberOfEvents(_maxNumberOfEvents), useMCtruth(_useMCtruth), weight(1)
    {
        TH1::SetDefaultSumw2();
        if (reweightFileName != "none"){
            std::shared_ptr<TFile> reweightFile(new TFile(reweightFileName.c_str(),"READ"));
            if(reweightFile->IsZombie()){
                std::ostringstream ss;
                ss << "reweight file " << reweightFileName << " not found." ;
                throw std::runtime_error(ss.str());
            }
            TObject* originalWeights = reweightFile->Get("weights");
            if (!originalWeights)
                throw std::runtime_error("histograms with weights not found");
            outputFile->cd();
            weights = std::shared_ptr<TH1D>(static_cast<TH1D*>(originalWeights->Clone("PUweights")));
        }

    }

    virtual ~BaseAnalyzer() {}

    virtual void Run()
    {
        size_t n = 0;
        for(; !maxNumberOfEvents || n < maxNumberOfEvents; ++n) {
            if(!treeExtractor.ExtractNext(event))
                break;
            timer.Report(n);
            try {
                ProcessEvent();
            } catch(cuts::cut_failed&){}
            GetAnaData().EventSelection().fill_selection(weight);
        }
        timer.Report(n, true);
        runReport.Report();
    }

protected:
    typedef std::function< Candidate (size_t, cuts::ObjectSelector&, bool, root_ext::AnalyzerData&) > BaseSelector;

    virtual SignalAnalyzerData& GetAnaData() = 0;

    virtual void ProcessEvent()
    {
        runReport.AddEvent(event.eventId());
        if (weights){
            const ntuple::Event& eventInfo = event.eventInfo();
            bool foundBX = false;
            for (unsigned n = 0; n < eventInfo.bunchCrossing.size(); ++n){
                if (eventInfo.bunchCrossing.at(n) == 0){
                    SetWeight(eventInfo.nPU.at(n));
                    foundBX = true;
                    break;
                }
            }
            if (!foundBX)
                throw std::runtime_error("in-time BX not found");
        }
    }

    void SetWeight(Int_t nPU)
    {
        static const Int_t maxAvailablePU = 70;
        static const double defaultWeight = 1.0;
        if (!weights) return;
        const Int_t bin = weights->FindBin(nPU);
        if (bin < 1 || bin > maxAvailablePU){
//            std::ostringstream ss;
//            ss << "no weight for this nPU " << nPU;
            //throw std::runtime_error(ss.str());
            weight = defaultWeight;
        }
        else
            weight = weights->GetBinContent(bin);
    }

    bool HaveTriggerFired(const std::vector<std::string>& interestingPaths) const
    {
        for (const ntuple::Trigger& trigger : event.triggers()){
            for (unsigned n = 0; n < trigger.hltpaths.size(); ++n){
                for (unsigned k = 0; k < interestingPaths.size(); ++k){
                    const std::string& triggerPath = trigger.hltpaths.at(n);
                    const std::string& interestingPath = interestingPaths.at(k);
                    std::size_t found = triggerPath.find(interestingPath);
                    if (found != std::string::npos && trigger.hltresults.at(n) == 1 &&
                            trigger.hltprescales.at(n) == 1) return true;
                }
            }
        }
        return false;
    }

    template<typename ObjectType, typename BaseSelectorType>
    std::vector<ObjectType> CollectObjects(cuts::ObjectSelector& objectSelector, size_t n_objects,
                                           const BaseSelectorType& base_selector, const std::string& hist_name)
    {
        const auto selector = [&](size_t id) -> ObjectType
            { return base_selector(id, objectSelector, true, anaDataBeforeCut); };

        const auto selected = objectSelector.collect_objects<ObjectType>(weight, n_objects,selector);
        for(const ObjectType& candidate : selected)
            base_selector(candidate.index, objectSelector, false, anaDataAfterCut);
        GetAnaData().N_objects(hist_name).Fill(selected.size(),weight);
        GetAnaData().N_objects(hist_name + "_ntuple").Fill(n_objects,weight);
        return selected;
    }

    template<typename BaseSelectorMethod>
    CandidateVector CollectObjects(cuts::ObjectSelector& objectSelector, size_t n_objects, bool signal,
                                   BaseSelectorMethod signal_selector_method,
                                   BaseSelectorMethod bkg_selector_method, const std::string& hist_name)
    {
        const BaseSelector base_selector_signal = [&](unsigned id, cuts::ObjectSelector& _objectSelector, bool enabled,
                root_ext::AnalyzerData& _anaData)
                -> Candidate { return (this->*signal_selector_method)(id, _objectSelector, enabled, _anaData); };
        const BaseSelector base_selector_bkg = [&](unsigned id, cuts::ObjectSelector& _objectSelector, bool enabled,
                root_ext::AnalyzerData& _anaData)
                -> Candidate { return (this->*bkg_selector_method)(id, _objectSelector, enabled, _anaData); };
        const auto base_selector = signal ? base_selector_signal : base_selector_bkg;
        const std::string real_hist_name = signal ? hist_name : hist_name + "_bkg";
        return CollectObjects<Candidate>(objectSelector, n_objects, base_selector, real_hist_name);

    }

    CandidateVector CollectMuons(bool signal = true)
    {
        auto& objectSelector = signal ? GetAnaData().MuonSelection() : GetAnaData().MuonSelectionBkg();
        return CollectObjects(objectSelector, event.muons().size(), signal,
                              &BaseAnalyzer::SelectMuon, &BaseAnalyzer::SelectBackgroundMuon, "muons");
    }

    CandidateVector CollectTaus(bool signal = true)
    {
        auto& objectSelector = signal ? GetAnaData().TauSelection() : GetAnaData().TauSelectionBkg();
        return CollectObjects(objectSelector, event.taus().size(), signal,
                              &BaseAnalyzer::SelectTau, &BaseAnalyzer::SelectBackgroundTau, "taus");
    }

    CandidateVector CollectElectrons(bool signal = true)
    {
        auto& objectSelector = signal ? GetAnaData().ElectronSelection() : GetAnaData().ElectronSelectionBkg();
        return CollectObjects(objectSelector, event.electrons().size(), signal,
                              &BaseAnalyzer::SelectElectron, &BaseAnalyzer::SelectBackgroundElectron, "electrons");
    }

    CandidateVector CollectBJets(double csv, const std::string& selection_label, bool signal = true)
    {
        const BaseSelector base_selector_signal = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBJet(id, _objectSelector, enabled, _anaData, csv, selection_label); };
        const BaseSelector base_selector_bkg = [&](unsigned id, cuts::ObjectSelector& _objectSelector,
                bool enabled, root_ext::AnalyzerData& _anaData)
                -> Candidate { return SelectBackgroundBJet(id, _objectSelector, enabled, _anaData); };
        const auto base_selector = signal ? base_selector_signal : base_selector_bkg;
        auto& objectSelector = signal ? GetAnaData().BJetSelection(selection_label)
                                          : GetAnaData().BJetSelectionBkg(selection_label);
        return CollectObjects<Candidate>(objectSelector, event.jets().size(), base_selector, "bjets");
    }

    VertexVector CollectVertices()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector& _objectSelector, bool enabled,
                root_ext::AnalyzerData& _anaData)
                -> Vertex { return SelectVertex(id, _objectSelector, enabled, _anaData); };

        return CollectObjects<Vertex>(GetAnaData().VertexSelection(), event.vertices().size(), base_selector, "vertices");
    }

    virtual Candidate SelectMuon(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                 root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Muon selection for signal not implemented");
    }
    virtual Candidate SelectTau(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                root_ext::AnalyzerData& _anaData) {
        throw std::runtime_error("Tau selection for signal not implemented");
    }
    virtual Candidate SelectElectron(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                     root_ext::AnalyzerData& _anaData){
        throw std::runtime_error("Electron selection for signal not implemented");
    }


    Vertex SelectVertex(size_t id, cuts::ObjectSelector& objectSelector, bool enabled, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::vertex;
        const std::string selection_label = "vertex";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Vertex& object = event.vertices().at(id);

        cut(true, ">0 vertex");
        cut(X(ndf, 350, 0.0, 350.0) >= ndf, "ndf");
        cut(std::abs( X(z, 600, -30.0, 30.0) ) < z, "z");
        const double r_vertex = std::sqrt(object.x*object.x+object.y*object.y);
        cut(std::abs( Y(r_vertex, 500, 0.0, 0.5) ) < r, "r");

        return analysis::Vertex(id, object);
    }

    Candidate SelectBJet(size_t id, cuts::ObjectSelector& objectSelector, bool enabled, root_ext::AnalyzerData& _anaData,
                         double csv, const std::string& _selection_label)
    {
        using namespace cuts::Htautau_Summer13::btag::signal;
        cuts::Cutter cut(objectSelector, enabled);
        const std::string selection_label = "bjet_" + _selection_label;
        const ntuple::Jet& object = event.jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags, 130, -11.0, 2.0) > csv, "CSV");

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
    }

    Candidate SelectBackgroundBJet(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                   root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::btag::veto;
        const std::string selection_label = "bjet_bkg";
        cuts::Cutter cut(objectSelector, enabled);

        const ntuple::Jet& object = event.jets().at(id);
        cut(true, ">0 b-jet cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(combinedSecondaryVertexBJetTags, 130, -11.0, 2.0) > CSV, "CSV");
        cut(X(passLooseID, 2, -0.5, 1.5) == passLooseID, "passLooseID");

        return analysis::Candidate(analysis::Candidate::Bjet, id, object);
    }

    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                         root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::electronID::veto;
        const std::string selection_label = "electron_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Electron& object = event.electrons().at(id);

        cut(true, ">0 ele cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        const double eta = std::abs( X(eta, 120, -6.0, 6.0) );
        cut(eta < eta_high && (eta < cuts::Htautau_Summer13::electronID::eta_CrackVeto_low ||
                               eta > cuts::Htautau_Summer13::electronID::eta_CrackVeto_high), "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig, 300, -1.5, 1.5) > MVApogNonTrig[pt_index][eta_index], "mva");

        return analysis::Candidate(analysis::Candidate::Electron, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                     root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::muonID::veto;
        const std::string selection_label = "muon_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Muon& object = event.muons().at(id);

        cut(true, ">0 mu cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(isTightMuon, 2, -0.5, 1.5) == isTightMuon, "tight");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");
        cut(X(pfRelIso, 1000, 0.0, 100.0) < pFRelIso, "pFRelIso");

        return analysis::Candidate(analysis::Candidate::Mu, id, object,object.charge);
    }

    virtual analysis::Candidate SelectBackgroundTau(size_t id, cuts::ObjectSelector& objectSelector, bool enabled,
                                                    root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::Htautau_Summer13::tauID::veto;
        const std::string selection_label = "tau_bkg";
        cuts::Cutter cut(objectSelector, enabled);
        const ntuple::Tau& object = event.taus().at(id);

        cut(true, ">0 tau cand");
        cut(X(pt, 1000, 0.0, 1000.0) > pt, "pt");
        cut(std::abs( X(eta, 120, -6.0, 6.0) ) < eta, "eta");
        cut(X(decayModeFinding, 2, -0.5, 1.5) > decayModeFinding, "decay_mode");
        cut(X(byLooseCombinedIsolationDeltaBetaCorr3Hits, 2, -0.5, 1.5) >
            LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ, 6000, 0.0, 60.0)  < dz, "dz");

        return analysis::Candidate(analysis::Candidate::Tau, id, object,object.charge);
    }


    CandidateVector FindCompatibleObjects(const CandidateVector& objects1, const CandidateVector& objects2,
                                          double minDeltaR, Candidate::Type type, const std::string& hist_name,
                                          int expectedCharge = Candidate::UnknownCharge())
    {
        CandidateVector result;
        for(const Candidate& object1 : objects1) {
            for(const Candidate& object2 : objects2) {
                if(object2.momentum.DeltaR(object1.momentum) > minDeltaR) {
                    const Candidate candidate(type, object1, object2);
                    if (candidate.charge != expectedCharge) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),weight);
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),weight);
        return result;
    }


    CandidateVector FindCompatibleObjects(const CandidateVector& objects, double minDeltaR, Candidate::Type type,
                                          const std::string& hist_name, int expectedCharge = Candidate::UnknownCharge())
    {
        CandidateVector result;
        for (unsigned n = 0; n < objects.size(); ++n){
            for (unsigned k = n+1; k < objects.size(); ++k){
                if(objects.at(n).momentum.DeltaR(objects.at(k).momentum) > minDeltaR) {
                    const Candidate candidate(type, objects.at(n), objects.at(k));
                    if (candidate.charge != expectedCharge) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),weight);
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),weight);
        return result;
    }

    CandidateVector FilterBackground(const CandidateVector& candidates, const CandidateVector& backgroundCandidates,
                                            double minDeltaR, const std::string& hist_name)
    {
        CandidateVector result;
        for (const Candidate& candidate : candidates){
            if (FilterBackground(candidate, backgroundCandidates,minDeltaR)){
                result.push_back(candidate);
                GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),weight);
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),weight);
        return result;
    }

    static bool FilterBackground(const Candidate& candidate, const CandidateVector& backgroundCandidates,
                                            double minDeltaR)
    {
        if(candidate.finalStateDaughters.size()) {
            for (const Candidate& bkg_candidate : backgroundCandidates) {
                bool hasMatchedDaughter = false;
                for(const Candidate* daughter : candidate.finalStateDaughters) {
                    if (bkg_candidate.momentum.DeltaR(daughter->momentum) <= minDeltaR) {
                        hasMatchedDaughter = true;
                        break;
                    }
                }
                if(!hasMatchedDaughter)
                    return false;
            }
            return true;
        }
        for (const Candidate& bkg_candidate : backgroundCandidates){
            if (bkg_candidate.momentum.DeltaR(candidate.momentum) > minDeltaR)
                return false;
        }
        return true;
    }

    CandidateVector ApplyVetos(const CandidateVector& Resonances, cuts::Cutter& cut)
    {
        const auto electrons_bkg = CollectElectrons(false);
        const auto resonances_noEle = FilterBackground(Resonances,electrons_bkg,
                              cuts::Htautau_Summer13::electronID::veto::deltaR_signalObjects,"resonances_noEle");
        cut(resonances_noEle.size(), "no_electrons");


        const auto muons_bkg = CollectMuons(false);
        const auto resonances_noMu = FilterBackground(resonances_noEle,muons_bkg,
                                 cuts::Htautau_Summer13::muonID::veto::deltaR_signalObjects,"resonances_noMu");
        cut(resonances_noMu.size(), "no_muons");


        const auto bjets_bkg = CollectBJets(cuts::Htautau_Summer13::btag::CSVL, "loose",false);
        const auto resonances_noBjets = FilterBackground(resonances_noMu,bjets_bkg,
                              cuts::Htautau_Summer13::btag::veto::deltaR_signalObjects,"resonances_noBjets");
        cut(resonances_noBjets.size(), "no_bjets");


        const auto taus_bkg = CollectTaus(false);
        const auto resonances_noTau = FilterBackground(resonances_noBjets,taus_bkg,
                                 cuts::Htautau_Summer13::tauID::veto::deltaR_signalObjects,"resonances_noTau");
        cut(resonances_noTau.size(), "no_taus");

        return resonances_noTau;
    }

    void HistogramAfterFinalSelection(const CandidateVector& finalSelectionResonances)
    {
        for (const Candidate& resonance : finalSelectionResonances){
            for (const Candidate* finalStateDaughter : resonance.finalStateDaughters){
                if (finalStateDaughter->type == Candidate::Mu)
                    SelectMuon(finalStateDaughter->index,GetAnaData().MuonSelection(), false, anaDataFinalSelection);
                if (finalStateDaughter->type == Candidate::Tau)
                    SelectTau(finalStateDaughter->index,GetAnaData().TauSelection(), false, anaDataFinalSelection);
                if (finalStateDaughter->type == Candidate::Bjet)
                    SelectBJet(finalStateDaughter->index,GetAnaData().BJetSelection("loose"), false, anaDataFinalSelection,
                               cuts::Htautau_Summer13::btag::CSVL, "loose");
                if (finalStateDaughter->type == Candidate::Electron)
                    SelectElectron(finalStateDaughter->index,GetAnaData().ElectronSelection(), false, anaDataFinalSelection);
            }
            GetAnaData().Mass("Final_H_tautau").Fill(resonance.daughters.at(0)->momentum.M(),weight);
            GetAnaData().Mass("Final_H_bb").Fill(resonance.daughters.at(1)->momentum.M(),weight);
        }
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
    Timer timer;
    EventDescriptor event;
    TreeExtractor treeExtractor;
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
    RunReport runReport;
    size_t maxNumberOfEvents;
    bool useMCtruth;
    GenEvent genEvent;
    Vertex primaryVertex;
    double weight;
    std::shared_ptr<TH1D> weights;
};

} // analysis
