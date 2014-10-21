/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau and H->tautau analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-20 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#ifndef __APPLE__
#define override
#endif

#include <iomanip>
#include <functional>
#include <string>
#include <iostream>

#include "AnalysisBase/include/TreeExtractor.h"
#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/CutTools.h"
#include "AnalysisBase/include/GenParticle.h"
#include "AnalysisBase/include/MCfinalState.h"
#include "AnalysisBase/include/Candidate.h"
#include "AnalysisBase/include/RunReport.h"
#include "AnalysisBase/include/Tools.h"
#include "AnalysisBase/include/AnalysisTools.h"
#include "AnalysisBase/include/AnalysisTypes.h"

#include "Htautau_Summer13.h"
#include "Htautau_TriggerEfficiency.h"
#include "BTagWeight.h"
#include "Config.h"
#include "custom_cuts.h"
#include "MvaMet.h"
#include "RecoilCorrection.h"
#include "SVfit.h"
#include "KinFit.h"

#define SELECTION_ENTRY(name) \
    ENTRY_1D(cuts::ObjectSelector, name) \
    /**/

/*
#define X(name, ...) \
    cuts::fill_histogram( object.name, _anaData.Get(&object.name, #name, selection_label), eventWeight)

#define Y(name, ...) \
    cuts::fill_histogram( name, _anaData.Get(&name, #name, selection_label), eventWeight)
*/

#define X(name, ...) \
    cuts::fill_histogram( object.name, _anaData.Get((TH1D*)nullptr, #name, selection_label, ##__VA_ARGS__), eventWeight)

#define Y(name, ...) \
    cuts::fill_histogram( name, _anaData.Get((TH1D*)nullptr, #name, selection_label, ##__VA_ARGS__), eventWeight)

namespace analysis {

class BaseAnalyzerData : public root_ext::AnalyzerData {

public:
    BaseAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(DeltaRmin1_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin1_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRbjets_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmin_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_visible, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_original, 600, 0, 6)
    TH1D_ENTRY(deltaPtMax, 250, 0, 5)
    TH1D_ENTRY(deltaPtMax_vis, 250, 0, 5)
    TH1D_ENTRY(MinPtBjetsMC, 20, 0, 200)
    TH1D_ENTRY(MassBB_MC, 30, 0, 300)
    TH1D_ENTRY(MassBB_MCvis, 30, 0, 300)
    TH1D_ENTRY(goodElectronsFromZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(goodMuonsFromZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZmm, 5, -0.5, 4.5)

};

class BaseAnalyzer {
public:
    BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName, const std::string& configFileName,
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : config(configFileName),
          outputFile(new TFile(outputFileName.c_str(), "RECREATE")),
          anaDataBeforeCut(*outputFile, "before_cut"), anaDataAfterCut(*outputFile, "after_cut"),
          anaDataFinalSelection(*outputFile, "final_selection"),
          maxNumberOfEvents(_maxNumberOfEvents),
          eventWeight(1), PUweight(1), triggerWeight(1), IDweight(1), IsoWeight(1),
          mvaMetProducer(config.MvaMet_dZcut(), config.MvaMet_inputFileNameU(), config.MvaMet_inputFileNameDPhi(),
                         config.MvaMet_inputFileNameCovU1(), config.MvaMet_inputFileNameCovU2())
    {
        if ( _prefix != "external" ){
            timer = std::shared_ptr<tools::Timer>(new tools::Timer(config.ReportInterval()));
            treeExtractor = std::shared_ptr<TreeExtractor>(
                        new TreeExtractor(_prefix == "none" ? "" : _prefix, inputFileName, config.extractMCtruth(),
                                          config.MaxTreeVersion()));

        }
        TH1::SetDefaultSumw2();
        if(config.ApplyPUreweight()){
            //std::cout << "I'm here" << std::endl;
            pu_weights = LoadPUWeights(config.PUreweight_fileName(), outputFile);
        }
        if(config.ApplyRecoilCorrection())
            recoilCorrectionProducer = std::shared_ptr<RecoilCorrectionProducer>(
                        new RecoilCorrectionProducer(config.RecoilCorrection_fileCorrectTo(),
                                                     config.RecoilCorrection_fileZmmData(),
                                                     config.RecoilCorrection_fileZmmMC()));
    }

    virtual ~BaseAnalyzer() {}

    virtual void Run()
    {
        size_t n = 0;
        auto _event = std::shared_ptr<EventDescriptor>(new EventDescriptor());
        if (!treeExtractor || !timer)
            throw exception("treeExtractor not initialized");
        for(; ( !maxNumberOfEvents || n < maxNumberOfEvents ) && treeExtractor->ExtractNext(*_event); ++n) {
            timer->Report(n);
//            std::cout << "event = " << _event->eventId().eventId << std::endl;
            if(config.RunSingleEvent() && _event->eventId().eventId != config.SingleEventId()) continue;
            try {
                ProcessEvent(_event);
            } catch(cuts::cut_failed&){}
            GetAnaData().Selection("event").fill_selection(eventWeight);
            if(config.RunSingleEvent()) break;
        }
        timer->Report(n, true);
    }


    virtual void ProcessEvent(std::shared_ptr<const EventDescriptor> _event)
    {
        event = _event;
        GetAnaData().getOutputFile().cd();
        eventWeight = 1;
        if (pu_weights){
            const ntuple::Event& eventInfo = event->eventInfo();
            const size_t bxIndex = tools::find_index(eventInfo.bunchCrossing, 0);
            if(bxIndex >= eventInfo.bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            //SetPUWeight(eventInfo.nPU.at(bxIndex));
            SetPUWeight(eventInfo.trueNInt.at(bxIndex));
        }
        eventWeight *= PUweight;
        if (config.isDYEmbeddedSample()) eventWeight *= event->genEvent().embeddedWeight;
    }

    double GetEventWeight() const
    {
        return eventWeight;
    }

protected:
    virtual BaseAnalyzerData& GetAnaData() = 0;

    void SetPUWeight(float nPU)
    {
        if (!pu_weights) return;
        const Int_t bin = pu_weights->FindBin(nPU);
        const Int_t maxBin = pu_weights->FindBin(config.PUreweight_maxAvailablePU());
        const bool goodBin = bin >= 1 && bin <= maxBin;
        PUweight = goodBin ? pu_weights->GetBinContent(bin) : config.PUreweight_defaultWeight();
    }

    void CalculateFullEventWeight(const Candidate& candidate)
    {
        CalculateTriggerWeights(candidate);
        CalculateIsoWeights(candidate);
        CalculateIdWeights(candidate);
        CalculateDMWeights(candidate);
        CalculateFakeWeights(candidate);
        for (double weight : triggerWeights){
            eventWeight *= weight;
        }
        for (double weight : IsoWeights){
            eventWeight *= weight;
        }
        for (double weight : IDweights){
            eventWeight *= weight;
        }
        for (double weight : DMweights){
            eventWeight *= weight;
        }
//        for (double weight : fakeWeights){
//            eventWeight *= weight;
//        }
    }

    virtual void CalculateTriggerWeights(const Candidate& candidate) = 0;
    virtual void CalculateIsoWeights(const Candidate& candidate) = 0;
    virtual void CalculateIdWeights(const Candidate& candidate) = 0;
    virtual void CalculateDMWeights(const Candidate& candidate) = 0;
    virtual void CalculateFakeWeights(const Candidate& candidate) = 0;

    bool HaveTriggerFired(const std::vector<std::string>& interestinghltPaths) const
    {
        for (const ntuple::Trigger& trigger : event->triggers()){
            for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
                if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1)
                    return true;
            }
        }
        return false;
    }

    std::vector<std::string> CollectPathsForTriggerFired(const std::vector<std::string>& interestinghltPaths)
    {
        std::vector<std::string> firedPaths;
        for (const ntuple::Trigger& trigger : event->triggers()){
            for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
                if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1){
                    std::string firedPath = trigger.hltpaths.at(n);
                    firedPaths.push_back(firedPath);
                }
            }
        }
        return firedPaths;
    }

    template<typename ObjectType, typename BaseSelectorType>
    std::vector<ObjectType> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                           size_t n_objects)
    {
        cuts::ObjectSelector& objectSelector = GetAnaData().Selection(selection_label);
        const auto selector = [&](size_t id) -> ObjectType
            { return base_selector(id, &objectSelector, anaDataBeforeCut, selection_label); };

        const auto selected = objectSelector.collect_objects<ObjectType>(eventWeight, n_objects,selector);
        for(const ObjectType& candidate : selected)
            base_selector(candidate.index, nullptr, anaDataAfterCut, selection_label);
        GetAnaData().N_objects(selection_label).Fill(selected.size(),eventWeight);
        GetAnaData().N_objects(selection_label + "_ntuple").Fill(n_objects,eventWeight);
        return selected;
    }

    template<typename BaseSelectorMethod>
    CandidateVector CollectCandidateObjects(const std::string& selection_label, BaseSelectorMethod selector_method,
                                            size_t n_objects)
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> Candidate
            { return (this->*selector_method)(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<Candidate>(selection_label, base_selector, n_objects);
    }

    CandidateVector CollectMuons()
    {
        return CollectCandidateObjects("muons", &BaseAnalyzer::SelectMuon, event->muons().size());
    }

    CandidateVector CollectSignalMuons()
    {
        return CollectCandidateObjects("muons", &BaseAnalyzer::SelectSignalMuon, event->muons().size());
    }

    virtual Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                 const std::string& selection_label)
    {
        throw std::runtime_error("Muon selection for signal not implemented");
    }

    virtual Candidate SelectSignalMuon(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                 const std::string& selection_label)
    {
        throw std::runtime_error("Signal muon selection for signal not implemented");
    }

    CandidateVector CollectBackgroundMuons()
    {
        return CollectCandidateObjects("muons_bkg", &BaseAnalyzer::SelectBackgroundMuon, event->muons().size());
    }

    virtual Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        throw std::runtime_error("Muon selection for background not implemented");
    }

    CandidateVector CollectTaus()
    {
        return CollectCandidateObjects("taus", &BaseAnalyzer::SelectTau, event->taus().size());
    }

    CandidateVector CollectSignalTaus()
    {
        return CollectCandidateObjects("taus", &BaseAnalyzer::SelectSignalTau, event->taus().size());
    }

    virtual Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                const std::string& selection_label)
    {
        throw std::runtime_error("Tau selection for signal not implemented");
    }

    virtual Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                const std::string& selection_label)
    {
        throw std::runtime_error("Signal tau selection for signal not implemented");
    }

    CandidateVector CollectElectrons()
    {
        return CollectCandidateObjects("electrons", &BaseAnalyzer::SelectElectron, event->electrons().size());
    }

    CandidateVector CollectSignalElectrons()
    {
        return CollectCandidateObjects("electrons", &BaseAnalyzer::SelectSignalElectron, event->electrons().size());
    }

    virtual Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                     const std::string& selection_label)
    {
        throw std::runtime_error("Electron selection for signal not implemented");
    }

    virtual Candidate SelectSignalElectron(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                           const std::string& selection_label)
    {
        throw std::runtime_error("Electron selection for signal not implemented");
    }

    CandidateVector CollectBackgroundElectrons()
    {
        return CollectCandidateObjects("electrons_bkg", &BaseAnalyzer::SelectBackgroundElectron,
                                       event->electrons().size());
    }

    virtual Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                               root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        throw std::runtime_error("Electron selection for background not implemented");
    }

    VertexVector CollectVertices()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> Vertex
            { return SelectVertex(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<Vertex>("vertices", base_selector, event->vertices().size());
    }

    Vertex SelectVertex(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                        const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::vertex;
        cuts::Cutter cut(objectSelector);
        const ntuple::Vertex& object = event->vertices().at(id);

        cut(true, ">0 vertex");
        cut(X(ndf) > ndf, "ndf");
        cut(std::abs( X(z) ) < z, "z");
        const double r_vertex = std::sqrt(object.x*object.x+object.y*object.y);
        cut(std::abs( Y(r_vertex) ) < r, "r");

        return analysis::Vertex(id, object);
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
                    if (expectedCharge != Candidate::UnknownCharge() && candidate.charge != expectedCharge) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),eventWeight);
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),eventWeight);
        return result;
    }


    CandidateVector FindCompatibleObjects(const CandidateVector& objects, double minDeltaR, Candidate::Type type,
                                          const std::string& hist_name, int expectedCharge = Candidate::UnknownCharge())
    {
        CandidateVector result;
        for (unsigned n = 0; n < objects.size(); ++n){
            for (unsigned k = n+1; k < objects.size(); ++k){
//                std::cout << "first tau momentum " << objects.at(n).momentum << std::endl;
//                std::cout << "second tau momentum " << objects.at(k).momentum << std::endl;
//                std::cout << "DeltaR " << objects.at(n).momentum.DeltaR(objects.at(k).momentum) << std::endl;
//                std::cout << "first tau charge " << objects.at(n).charge << std::endl;
//                std::cout << "second tau charge " << objects.at(k).charge << std::endl;
                if(objects.at(n).momentum.DeltaR(objects.at(k).momentum) > minDeltaR) {
                    const Candidate candidate(type, objects.at(n), objects.at(k));
                    if (expectedCharge != Candidate::UnknownCharge() && candidate.charge != expectedCharge ) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),eventWeight);
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),eventWeight);
        return result;
    }

    CandidateVector FilterCompatibleObjects(const CandidateVector& objectsToFilter,
                                            const Candidate& referenceObject,
                                          double minDeltaR)
    {
        CandidateVector result;
        for(const Candidate& filterObject : objectsToFilter) {
            bool allDaughterPassed = true;
            for (const Candidate& daughter : referenceObject.finalStateDaughters){
                if(filterObject.momentum.DeltaR(daughter.momentum) <= minDeltaR) {
                    allDaughterPassed = false;
                    break;
                }
            }
            if (allDaughterPassed) result.push_back(filterObject);
        }
        return result;
    }

    ntuple::TauVector ApplyTauCorrections(const VisibleGenObjectVector& hadronic_taus, bool useLegacyCorrections)
    {
        using namespace cuts::Htautau_Summer13::tauCorrections;

        ntuple::TauVector correctedTaus;

        for(const ntuple::Tau& tau : event->taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);

            const bool hasMCmatch = FindMatchedObjects(momentum, hadronic_taus, deltaR).size() != 0;
            const double scaleFactor = MomentumScaleFactor(hasMCmatch, momentum.Pt(),
                                   ntuple::tau_id::ConvertToHadronicDecayMode(tau.decayMode), useLegacyCorrections);
            const TLorentzVector correctedMomentum = momentum * scaleFactor;
            ntuple::Tau correctedTau(tau);
            correctedTau.pt = correctedMomentum.Pt();
            correctedTau.eta = correctedMomentum.Eta();
            correctedTau.phi = correctedMomentum.Phi();
            correctedTau.mass = correctedMomentum.M();
            correctedTaus.push_back(correctedTau);
        }
        return correctedTaus;
    }

    ntuple::MET ApplyTauCorrectionsToMVAMET(const ntuple::MET& metMVA, const ntuple::TauVector& correctedTaus)
    {
        TLorentzVector sumCorrectedTaus, sumTaus;
        for(const ntuple::Tau& tau : event->taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
            sumTaus += momentum;
        }

        for (const ntuple::Tau& correctedTau : correctedTaus) {
            TLorentzVector correctedMomentum;
            correctedMomentum.SetPtEtaPhiM(correctedTau.pt, correctedTau.eta, correctedTau.phi,correctedTau.mass);
            sumCorrectedTaus += correctedMomentum;
        }

        TLorentzVector met, metCorrected;
        met.SetPtEtaPhiM(metMVA.pt, 0, metMVA.phi, 0.);
        metCorrected = met + sumTaus - sumCorrectedTaus;
        ntuple::MET correctedMET = metMVA;
        correctedMET.pt = metCorrected.Pt();
        correctedMET.phi = metCorrected.Phi();
        return correctedMET;
    }

    analysis::CandidateVector ApplyTriggerMatch(const analysis::CandidateVector& higgses,
                                                        const std::vector<std::string>& hltPaths,
                                                        bool useStandardTriggerMatch)
    {
        analysis::CandidateVector triggeredHiggses;
        for (const auto& higgs : higgses){
            if(!useStandardTriggerMatch && analysis::HaveTriggerMatched(event->triggerObjects(), hltPaths, higgs,
                                                                        cuts::Htautau_Summer13::DeltaR_triggerMatch))
                triggeredHiggses.push_back(higgs);
            if (useStandardTriggerMatch && analysis::HaveTriggerMatched(*event,hltPaths,higgs))
                triggeredHiggses.push_back(higgs);
        }
        return triggeredHiggses;
    }

    ntuple::MET ApplyRecoilCorrections(const Candidate& higgs, const GenParticle* resonance,
                                       const size_t njets, const ntuple::MET& correctedMET)
    {
        if (config.ApplyRecoilCorrection()){
            if(!resonance && config.ApplyRecoilCorrectionForW())
                resonance = FindWboson();
            if(resonance)
                return recoilCorrectionProducer->ApplyCorrection(correctedMET, higgs.momentum, resonance->momentum, njets);
        }
        return correctedMET;
    }

    const GenParticle* FindWboson()
    {
        static const particles::ParticleCodes Wcode = { particles::W_plus };
        static const particles::ParticleCodes WDecay_tau = { particles::tau, particles::nu_tau };
        static const particles::ParticleCodes WDecay_electron = { particles::e, particles::nu_e };
        static const particles::ParticleCodes WDecay_muon = { particles::mu, particles::nu_mu };

        const analysis::GenParticleSet Wparticles_all = genEvent.GetParticles(Wcode);

        analysis::GenParticleSet Wparticles;
        for(const GenParticle* w : Wparticles_all) {
            if(w->mothers.size() == 1) {
                const GenParticle* mother = w->mothers.at(0);
                if(mother->pdg.Code == particles::W_plus && mother->status == particles::HardInteractionProduct)
                    Wparticles.insert(w);
            }
         }

        if (Wparticles.size() == 0) return nullptr;

        if (Wparticles.size() > 1)
            throw exception("more than 1 W in the event");

        const GenParticle* Wboson = *Wparticles.begin();
        while(Wboson->daughters.size() == 1 && Wboson->daughters.front()->pdg.Code == particles::W_plus)
            Wboson = Wboson->daughters.front();

        analysis::GenParticlePtrVector WProducts;
        if(analysis::FindDecayProducts(*Wboson, WDecay_tau, WProducts,true))
            return Wboson;
        if (!FindDecayProducts(*Wboson, WDecay_electron, WProducts,true)
                && !FindDecayProducts(*Wboson, WDecay_muon, WProducts,true))
            throw exception("not leptonic W decay");


        return nullptr;
    }

    const analysis::Candidate& SelectSemiLeptonicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const double first_Pt1 = first.daughters.at(0).momentum.Pt();
            const double first_Pt2 = first.daughters.at(1).momentum.Pt();
            const double first_sumPt = first_Pt1 + first_Pt2;
            const double second_Pt1 = second.daughters.at(0).momentum.Pt();
            const double second_Pt2 = second.daughters.at(1).momentum.Pt();
            const double second_sumPt = second_Pt1 + second_Pt2;

            return first_sumPt < second_sumPt;
        };
        return *std::max_element(higgses.begin(), higgses.end(), higgsSelector) ;
    }


protected:
    Config config;
    std::shared_ptr<tools::Timer> timer;
    std::shared_ptr<const EventDescriptor> event;
    std::shared_ptr<TreeExtractor> treeExtractor;
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
    size_t maxNumberOfEvents;
    GenEvent genEvent;
    Vertex primaryVertex;
    double eventWeight, PUweight, triggerWeight, IDweight, IsoWeight;
    std::vector<double> triggerWeights, IDweights, IsoWeights, DMweights, fakeWeights;
    std::shared_ptr<TH1D> pu_weights;
    MvaMetProducer mvaMetProducer;
    std::shared_ptr<RecoilCorrectionProducer> recoilCorrectionProducer;
};

} // analysis

namespace make_tools {
template<typename T>
struct Factory {
    static T* Make(int argc, char *argv[])
    {
        std::cout << "Command line: ";
        for(int n = 0; n < argc; ++n)
            std::cout << argv[n] << " ";
        std::cout << std::endl;
        if(argc < 4 || argc > 6)
            throw std::runtime_error("Invalid number of command line arguments.");

        int n = 0;
        const std::string inputFileName = argv[++n];
        const std::string outputFileName = argv[++n];
        const std::string configFileName = argv[++n];
        if(argc <= n)
            return new T(inputFileName, outputFileName, configFileName);

        const std::string prefix = argv[++n];
        if(argc <= n)
            return new T(inputFileName, outputFileName, configFileName, prefix);

        char c;
        size_t maxNumberOfEvents;
        std::istringstream ss_nEvents(argv[++n]);
        ss_nEvents >> c >> maxNumberOfEvents;
        if(c != '@')
            throw std::runtime_error("Bad command line format.");

        return new T(inputFileName, outputFileName, configFileName, prefix, maxNumberOfEvents);
    }
};
} // make_tools
