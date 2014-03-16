/*!
 * \file HHbbTauTau_analyzer.C
 * \brief Analysis HHbbTauTau.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-12 created
 */

#include <vector>
#include <set>
#include <type_traits>
#include <boost/shared_ptr.hpp>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TObject.h>

#define EventTree_cxx

namespace ntuple {
    using namespace std;

    #include "../include/EventTree.h"

    void EventTree::Loop(){} //just declared
}

namespace ntupleMC {
    using namespace std;


    #include "../include/EventTreeMC.h"

    void EventTree::Loop(){} //just declared
}

#include "../include/AnalyzerData.h"
#include "../include/Particles.h"
#include "../include/Htautau_Summer13.h"
#include "../include/CutTools.h"
#include "../include/GenParticle.h"


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
    TH1D_ENTRY_FIX(name##_effAbs, 1, n_bins, -0.5)\
    /**/

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
    SignalAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}
    ~SignalAnalyzerData()
    {
        Erase(Counter().Name());
        for(const auto& desc : selectionDescriptors){
            cuts::fill_relative_selection_histogram(*desc.second.event, *desc.second.eff_relative, desc.second.fix_bin);
            cuts::fill_absolute_selection_histogram(*desc.second.event, *desc.second.eff_absolute);
        }
    }

    SELECTION_ENTRY(EventSelection, 15, 5)
    SELECTION_ENTRY(MuonSelection, 15)
    SELECTION_ENTRY(TauSelection, 15)
    SELECTION_ENTRY(ElectronSelection, 15)
    SELECTION_ENTRY(BJetSelection, 15)

    TH1D_ENTRY_FIX(Counter, 1, 100, -0.5)
    ENTRY_1D(double, Radion_Mass)
    ENTRY_1D(double, Radion_Pt)
    ENTRY_1D(double, Radion_Eta)
    ENTRY_1D(double, Radion_Phi)
    ENTRY_1D(double, Mu_tau_mass)
    ENTRY_1D(double, BB_mass)
    ENTRY_1D(double,Tau_Pt_MC)
    ENTRY_1D(double,Mu_Pt_MC)
    ENTRY_1D(double,Bjets_Pt_MC)
    ENTRY_1D(double,Higgs_MuTau_MC_Pt)
    ENTRY_1D(double,Higgs_BB_MC_Pt)
    ENTRY_2D(double, DR_bjets_vs_HiggsPt_MC)
    ENTRY_2D(double, DR_Higgs_vs_RadionPt_MC)

private:
    std::map<root_ext::SmartHistogram<TH1D>*, SelectionDescriptor> selectionDescriptors;
};

using ntupleMC::EventTree;

namespace cuts {
namespace HHbbTauTau {

const unsigned N_bjet_Loose = 2;
const unsigned N_bjet_Medium = 1;
}
}

#define X(name, ...) \
    cuts::fill_histogram( event->name[id], _anaData.Get(&event->name[id], #name, ##__VA_ARGS__) )

class HHbbTauTau_analyzer
{
public:
    HHbbTauTau_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0, bool _useMCtruth = false)
        : anaData(outputFileName), anaDataBeforeCut(anaData.getOutputFile(), "before_cut"),
          anaDataAfterCut(anaData.getOutputFile(), "after_cut"),
          maxNumberOfEvents(_maxNumberOfEvents), useMCtruth(_useMCtruth)
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
    }

private:

    void ProcessEvent()
    {
        if (useMCtruth){
            static const analysis::ParticleCodes resonanceCodes = { particles::Radion };
            static const analysis::ParticleCodes resonanceDecay = { particles::Higgs, particles::Higgs };
            static const analysis::ParticleCodes2D HiggsDecays = { {particles::b, particles::b},
                                                                   {particles::tau, particles::tau}};
            static const analysis::ParticleCodes TauMuonicDecay = {particles::mu, particles::nu_mu, particles::nu_tau};
            static const analysis::ParticleCodes TauElectronDecay = {particles::e, particles::nu_e, particles::nu_tau};

            const analysis::GenEvent genEvent(*event);
            //std::cout << "N gen particles = " << genEvent.genParticles.size() << std::endl;

            const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
            if (resonances.size() != 1)
                throw std::runtime_error("not one resonance per event");

            const analysis::GenParticle& resonance = **resonances.begin();
            if (resonance.MissingDaughter)
                throw std::runtime_error("resonance has no daughters");
//            std::cout << "N resonance daughters = " << resonance.daughters.size() << std::endl;
//            if (resonance.daughters.size() != 2)
//                throw std::runtime_error("the candidate resonance has not 2 daughters");

//            genEvent.Print();
            analysis::GenParticleVector HiggsBosons;
            if(!analysis::FindDecayProducts(resonance, resonanceDecay,HiggsBosons))
                throw std::runtime_error("Resonance does not decay into 2 Higgs");

            analysis::GenParticleVector2D HiggsDecayProducts;
            analysis::GenParticleIndexVector HiggsIndexes;
            if(!analysis::FindDecayProducts2D(HiggsBosons,HiggsDecays,HiggsDecayProducts,HiggsIndexes))
                throw std::runtime_error("NOT HH -> bb tautau");

            const analysis::GenParticleVector& b_jets_MC = HiggsDecayProducts.at(0);
            const analysis::GenParticleVector& taus_MC = HiggsDecayProducts.at(1);

            const analysis::GenParticle* muon_MC = nullptr;
            const analysis::GenParticle* tau_jet_MC = nullptr;
            for (const analysis::GenParticle* tau_MC : taus_MC){
                analysis::GenParticleVector TauProducts;
                if (analysis::FindDecayProducts(*tau_MC,TauMuonicDecay,TauProducts)){
                    muon_MC = TauProducts.at(0);
                    continue;
                }
                if (!analysis::FindDecayProducts(*tau_MC,TauElectronDecay,TauProducts)){
                    tau_jet_MC = tau_MC;
                }
            }

            if (!muon_MC || !tau_jet_MC) {
                //std::cout << "not our final state mu-tau" << std::endl;
                return;
            }

            const analysis::GenParticle& Higgs_MuTau_MC = *HiggsBosons.at(HiggsIndexes.at(1));
            const analysis::GenParticle& Higgs_BB_MC = *HiggsBosons.at(HiggsIndexes.at(0));

            anaData.Radion_Mass().Fill(resonance.momentum.M());
            anaData.Radion_Pt().Fill(resonance.momentum.Pt());
            anaData.Radion_Eta().Fill(resonance.momentum.Eta());
            anaData.Radion_Phi().Fill(resonance.momentum.Phi());

            anaData.Tau_Pt_MC().Fill(tau_jet_MC->momentum.Pt());
            anaData.Mu_Pt_MC().Fill(muon_MC->momentum.Pt());
            anaData.Bjets_Pt_MC().Fill(b_jets_MC.at(0)->momentum.Pt());
            anaData.Bjets_Pt_MC().Fill(b_jets_MC.at(1)->momentum.Pt());
            anaData.Higgs_MuTau_MC_Pt().Fill(Higgs_MuTau_MC.momentum.Pt());
            anaData.Higgs_BB_MC_Pt().Fill(Higgs_BB_MC.momentum.Pt());
            anaData.DR_bjets_vs_HiggsPt_MC().Fill(Higgs_BB_MC.momentum.Pt(),
                                                  b_jets_MC.at(0)->momentum.DeltaR(b_jets_MC.at(1)->momentum));
            anaData.DR_Higgs_vs_RadionPt_MC().Fill(resonance.momentum.Pt(),
                                                   Higgs_MuTau_MC.momentum.DeltaR(Higgs_BB_MC.momentum));

        }

        cuts::Cutter cut(anaData.EventSelection(), anaData.EventSelection());

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
        cut.test(b_jets_loose.size() == 2, "2b_loose");
        cut.test(b_jets_loose.size() == 2 && b_jets_medium.size() >= 1, "1b_loose+1b_medium");
        if (cut.test(b_jets_medium.size() == 2, "2b_medium")){
            TLorentzVector bjet1_momentum;
            bjet1_momentum.SetPtEtaPhiE(event->Jet_pt[b_jets_medium.at(0)], event->Jet_eta[b_jets_medium.at(0)],
                    event->Jet_phi[b_jets_medium.at(0)],event->Jet_energy[b_jets_medium.at(0)]);
            TLorentzVector bjet2_momentum;
            bjet2_momentum.SetPtEtaPhiE(event->Jet_pt[b_jets_medium.at(1)], event->Jet_eta[b_jets_medium.at(1)],
                    event->Jet_phi[b_jets_medium.at(1)],event->Jet_energy[b_jets_medium.at(1)]);
            TLorentzVector bb_momentum = bjet1_momentum + bjet2_momentum;
            anaData.BB_mass().Fill(bb_momentum.M());
        }
        cut.test(b_jets_loose.size() >= 2, ">=2b_loose");
        cut.test(b_jets_loose.size() >= 2 && b_jets_medium.size() >= 1, ">=1b_loose+>=1b_medium");
        cut.test(b_jets_medium.size() >= 2, ">=2b_medium");

    }

    template<typename BaseSelector, typename ValueType>
    IndexVector CollectObjects(TH1D& selection_histogram, Int_t n_objects,
                               const BaseSelector base_selector, const ValueType* values_to_compare)
    {
        const auto selector = [&](unsigned id) { base_selector(id, true, anaDataBeforeCut); };

        const auto comparitor = [&](unsigned a, unsigned b) -> bool
            { return values_to_compare[a] > values_to_compare[b]; };

        const auto selected = cuts::collect_objects(anaData.Counter(), selection_histogram, n_objects, selector,
                                                    comparitor);
        for(Int_t id : selected) base_selector(id, false, anaDataAfterCut);
        return selected;
    }

    IndexVector CollectMuons()
    {
        const auto base_selector = [&](unsigned id, bool apply_cut, root_ext::AnalyzerData& _anaData)
            { SelectMuon(id, apply_cut, _anaData); };
        return CollectObjects(anaData.MuonSelection(), event->nMuon, base_selector, event->Muon_pt);
    }

    IndexVector CollectTaus()
    {
        const auto base_selector = [&](unsigned id, bool apply_cut, root_ext::AnalyzerData& _anaData)
            { SelectTau(id, apply_cut, _anaData); };
        return CollectObjects(anaData.TauSelection(), event->nTau, base_selector, event->Tau_pt);
    }

    IndexVector CollectElectrons()
    {
        const auto base_selector = [&](unsigned id, bool apply_cut, root_ext::AnalyzerData& _anaData)
            { SelectElectron(id, apply_cut, _anaData); };
        return CollectObjects(anaData.ElectronSelection(), event->nElectron, base_selector, event->Electron_pt);
    }

    IndexVector CollectBJets(double csv, const std::string& selection_label)
    {
        const auto base_selector = [&](unsigned id, bool apply_cut, root_ext::AnalyzerData& _anaData)
            { SelectBJet(id, csv, selection_label, apply_cut, _anaData); };
        return CollectObjects(anaData.BJetSelection(selection_label), event->nJet, base_selector, event->Jet_pt);
    }

    void SelectMuon(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::muonID;
        cuts::Cutter cut(anaData.Counter(), anaData.MuonSelection(), apply_cut);

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

    void SelectTau(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::tauID;
        cuts::Cutter cut(anaData.Counter(), anaData.TauSelection(), apply_cut);

        cut(true, ">0 tau cand");
        cut(X(Tau_pt) > pt, "pt");
        cut(std::abs( X(Tau_eta) ) < eta, "eta");
        cut(X(Tau_decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits) > LooseCombinedIsolationDeltaBetaCorr3Hits, "looseIso3Hits");
        cut(X(Tau_againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(Tau_againstElectronLoose) > againstElectronLoose, "vs_e_loose");
    }

    void SelectElectron(Int_t id, bool apply_cut, root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::electronID;
        cuts::Cutter cut(anaData.Counter(), anaData.ElectronSelection(), apply_cut);

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

    void SelectBJet(Int_t id, double csv, const std::string& selection_label, bool apply_cut,
                    root_ext::AnalyzerData& _anaData)
    {
        using namespace cuts::btag;
        cuts::Cutter cut(anaData.Counter(), anaData.BJetSelection(selection_label), apply_cut);

        cut(true, ">0 b-jet cand");
        cut(X(Jet_pt, selection_label) > pt, "pt");
        cut(std::abs( X(Jet_eta, selection_label) ) < eta, "eta");
        cut(X(Jet_combinedSecondaryVertexBTag, selection_label) > csv, "CSV");
    }

    IndexPairVector FindCompatibleLeptonCombinations(const IndexVector& muons, const IndexVector& taus)
    {
        IndexPairVector result;
        for(auto mu_id : muons) {
            TLorentzVector mu_momentum;
            mu_momentum.SetPtEtaPhiE(event->Muon_pt[mu_id], event->Muon_eta[mu_id], event->Muon_phi[mu_id],
                                     event->Muon_energy[mu_id]);
            for(auto tau_id : taus) {
                TLorentzVector tau_momentum;
                tau_momentum.SetPtEtaPhiE(event->Tau_pt[tau_id], event->Tau_eta[tau_id], event->Tau_phi[tau_id],
                                         event->Tau_energy[tau_id]);
                TLorentzVector mu_tau_momentum = mu_momentum + tau_momentum;
                if(tau_momentum.DeltaR(mu_momentum) > cuts::DeltaR_signalLeptons){
                    result.push_back(IndexPair(mu_id, tau_id));
                    anaData.Mu_tau_mass().Fill(mu_tau_momentum.M());
                }
            }
        }
        return result;
    }

private:
    EventTree* event;
    SignalAnalyzerData anaData;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut;
    Long64_t maxNumberOfEvents;
    bool useMCtruth;
};
