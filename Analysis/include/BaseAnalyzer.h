/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#pragma once
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
    SignalAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}
    ~SignalAnalyzerData()
    {
        Erase(Counter().Name());
        for(const auto& desc : selectionDescriptors){
            cuts::fill_relative_selection_histogram(*desc.second.event, *desc.second.eff_relative, desc.second.fix_bin);
            cuts::fill_absolute_selection_histogram(*desc.second.event, *desc.second.eff_absolute);
        }
    }

    TH1D_ENTRY_FIX(Counter, 1, 100, -0.5)
    ENTRY_1D(double, Resonance_Mass)
    ENTRY_1D(double, Resonance_Pt)
    ENTRY_1D(double, Resonance_Eta)
    ENTRY_1D(double, Resonance_Phi)
    ENTRY_1D(double,Bjets_Pt_MC)
    ENTRY_1D(double,Higgs_leptonic_MC_Pt)
    ENTRY_1D(double,Higgs_BB_MC_Pt)
    ENTRY_2D(double, DR_bjets_vs_HiggsPt_MC)
    ENTRY_2D(double, DR_Higgs_vs_ResonancePt_MC)

protected:
    std::map<root_ext::SmartHistogram<TH1D>*, SelectionDescriptor> selectionDescriptors;
};

#define X(name, ...) \
    cuts::fill_histogram( event->name[id], _anaData.Get(&event->name[id], #name, ##__VA_ARGS__) )

class BaseAnalyzer
{
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
        event = boost::shared_ptr<ntuple::EventTree>(new ntuple::EventTree(inputTree,useMCtruth));
        std::cout << "Starting analyzer...\n";
    }

    virtual ~BaseAnalyzer(){}

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
    analysis::CandidateVector CollectObjects(TH1D& selection_histogram, Int_t n_objects,
                               const BaseSelector base_selector)
    {
        const auto selector = [&](unsigned id) -> analysis::Candidate
            { return base_selector(id, true, anaDataBeforeCut); };


        const auto selected = cuts::collect_objects(GetAnaData().Counter(), selection_histogram, n_objects, selector);
        for(const analysis::Candidate& candidate : selected) base_selector(candidate.index, false, anaDataAfterCut);
        return selected;
    }



    void FindAnalysisFinalState(analysis::finalState::bbTauTau& finalState){
        static const analysis::ParticleCodes resonanceCodes = { particles::Radion };
        static const analysis::ParticleCodes resonanceDecay = { particles::Higgs, particles::Higgs };
        static const analysis::ParticleCodes2D HiggsDecays = { {particles::b, particles::b},
                                                               {particles::tau, particles::tau}};

        //const analysis::GenEvent genEvent(*event);
        genEvent = boost::shared_ptr<analysis::GenEvent>(new analysis::GenEvent(*event));
        //std::cout << "N gen particles = " << genEvent->genParticles.size() << std::endl;

        const analysis::GenParticleSet resonances = genEvent->GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        finalState.resonance = *resonances.begin();

//            genEvent->Print();
        analysis::GenParticleVector HiggsBosons;
        if(!analysis::FindDecayProducts(*finalState.resonance, resonanceDecay,HiggsBosons))
            throw std::runtime_error("Resonance does not decay into 2 Higgs");

        analysis::GenParticleVector2D HiggsDecayProducts;
        analysis::GenParticleIndexVector HiggsIndexes;
        if(!analysis::FindDecayProducts2D(HiggsBosons,HiggsDecays,HiggsDecayProducts,HiggsIndexes))
            throw std::runtime_error("NOT HH -> bb tautau");

        finalState.b_jets = HiggsDecayProducts.at(0);
        finalState.taus = HiggsDecayProducts.at(1);

        finalState.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(1));
        finalState.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(0));

        GetAnaData().Resonance_Mass().Fill(finalState.resonance->momentum.M());
        GetAnaData().Resonance_Pt().Fill(finalState.resonance->momentum.Pt());
        GetAnaData().Resonance_Eta().Fill(finalState.resonance->momentum.Eta());
        GetAnaData().Resonance_Phi().Fill(finalState.resonance->momentum.Phi());


        GetAnaData().Bjets_Pt_MC().Fill(finalState.b_jets.at(0)->momentum.Pt());
        GetAnaData().Bjets_Pt_MC().Fill(finalState.b_jets.at(1)->momentum.Pt());
        GetAnaData().Higgs_leptonic_MC_Pt().Fill(finalState.Higgs_TauTau->momentum.Pt());
        GetAnaData().Higgs_BB_MC_Pt().Fill(finalState.Higgs_BB->momentum.Pt());
        GetAnaData().DR_bjets_vs_HiggsPt_MC().Fill(finalState.Higgs_BB->momentum.Pt(),
                                              finalState.b_jets.at(0)->momentum.DeltaR(finalState.b_jets.at(1)->momentum));
        GetAnaData().DR_Higgs_vs_ResonancePt_MC().Fill(finalState.resonance->momentum.Pt(),
                                               finalState.Higgs_TauTau->momentum.DeltaR(finalState.Higgs_BB->momentum));

    }

protected:
    boost::shared_ptr<ntuple::EventTree> event;
    boost::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut;
    Long64_t maxNumberOfEvents;
    bool useMCtruth;
    boost::shared_ptr<analysis::GenEvent> genEvent;
};



