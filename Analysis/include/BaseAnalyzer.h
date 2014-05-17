/*!
 * \file BaseAnalyzer.h
 * \brief Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau and H->tautau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-20 created
 */

#pragma once
#include "TreeExtractor.h"
#include "AnalyzerData.h"
#include "Htautau_Summer13.h"
#include "CutTools.h"
#include "GenParticle.h"
#include "MCfinalState.h"
#include "Candidate.h"
#include "custom_cuts.h"
#include "PostRecoilCorrection.h"
#include "SVfit.h"
#include <iomanip>
#include <functional>
#include "RunReport.h"
#include "Tools.h"
#include "AnalysisTools.h"

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

class BaseAnalyzerData : public root_ext::AnalyzerData {

public:
    BaseAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
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
        if (reweightFileName != "none")
            weights = LoadPUWeights(reweightFileName, outputFile);
    }

    virtual ~BaseAnalyzer() {}

    virtual void Run()
    {
        size_t n = 0;
        for(; ( !maxNumberOfEvents || n < maxNumberOfEvents ) && treeExtractor.ExtractNext(event); ++n) {
            timer.Report(n);
            //if (event.eventId().eventId != 101874) continue;
            try {
                ProcessEvent();
            } catch(cuts::cut_failed&){}
            GetAnaData().Selection("event").fill_selection(weight);
            //break;
        }
        timer.Report(n, true);
        runReport.Report();
    }

protected:
    virtual BaseAnalyzerData& GetAnaData() = 0;

    virtual void ProcessEvent()
    {
        runReport.AddEvent(event.eventId());
        if (weights){
            const ntuple::Event& eventInfo = event.eventInfo();
            const size_t bxIndex = tools::find_index(eventInfo.bunchCrossing, 0);
            if(bxIndex >= eventInfo.bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            //SetWeight(eventInfo.nPU.at(bxIndex));
            SetWeight(eventInfo.trueNInt.at(bxIndex));
        }
    }

    void SetWeight(float nPU)
    {
        static const float maxAvailablePU = 60;
        static const double defaultWeight = 0.0;
        if (!weights) return;
        const Int_t bin = weights->FindBin(nPU);
		const Int_t maxBin = weights->FindBin(maxAvailablePU);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        weight = goodBin ? weights->GetBinContent(bin) : defaultWeight;
    }

    bool HaveTriggerFired(const std::vector<std::string>& interestinghltPaths) const
    {
        for (const ntuple::Trigger& trigger : event.triggers()){
            for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
                if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1)
                    return true;
            }
        }
        return false;
    }

    template<typename ObjectType, typename BaseSelectorType>
    std::vector<ObjectType> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                           size_t n_objects)
    {
        cuts::ObjectSelector& objectSelector = GetAnaData().Selection(selection_label);
        const auto selector = [&](size_t id) -> ObjectType
            { return base_selector(id, &objectSelector, anaDataBeforeCut, selection_label); };

        const auto selected = objectSelector.collect_objects<ObjectType>(weight, n_objects,selector);
        for(const ObjectType& candidate : selected)
            base_selector(candidate.index, nullptr, anaDataAfterCut, selection_label);
        GetAnaData().N_objects(selection_label).Fill(selected.size(),weight);
        GetAnaData().N_objects(selection_label + "_ntuple").Fill(n_objects,weight);
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
        return CollectCandidateObjects("muons", &BaseAnalyzer::SelectMuon, event.muons().size());
    }

    virtual Candidate SelectMuon(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                 const std::string& selection_label)
    {
        throw std::runtime_error("Muon selection for signal not implemented");
    }

    CandidateVector CollectBackgroundMuons()
    {
        return CollectCandidateObjects("muons_bkg", &BaseAnalyzer::SelectBackgroundMuon, event.muons().size());
    }

    virtual Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                           root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        throw std::runtime_error("Muon selection for background not implemented");
    }

    CandidateVector CollectTaus()
    {
        return CollectCandidateObjects("taus", &BaseAnalyzer::SelectTau, event.taus().size());
    }

    virtual Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                const std::string& selection_label)
    {
        throw std::runtime_error("Tau selection for signal not implemented");
    }


    CandidateVector CollectElectrons()
    {
        return CollectCandidateObjects("electrons", &BaseAnalyzer::SelectElectron, event.electrons().size());
    }

    virtual Candidate SelectElectron(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                                     const std::string& selection_label)
    {
        throw std::runtime_error("Electron selection for signal not implemented");
    }

    CandidateVector CollectBackgroundElectrons()
    {
        return CollectCandidateObjects("electrons_bkg", &BaseAnalyzer::SelectBackgroundElectron,
                                       event.electrons().size());
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
        return CollectObjects<Vertex>("vertices", base_selector, event.vertices().size());
    }

    Vertex SelectVertex(size_t id, cuts::ObjectSelector* objectSelector, root_ext::AnalyzerData& _anaData,
                        const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::vertex;
        cuts::Cutter cut(objectSelector);
        const ntuple::Vertex& object = event.vertices().at(id);

        cut(true, ">0 vertex");
        cut(X(ndf, 350, 0.0, 350.0) > ndf, "ndf");
        cut(std::abs( X(z, 600, -30.0, 30.0) ) < z, "z");
        const double r_vertex = std::sqrt(object.x*object.x+object.y*object.y);
        cut(std::abs( Y(r_vertex, 500, 0.0, 0.5) ) < r, "r");

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
//                std::cout << "first tau momentum " << objects.at(n).momentum << std::endl;
//                std::cout << "second tau momentum " << objects.at(k).momentum << std::endl;
//                std::cout << "DeltaR " << objects.at(n).momentum.DeltaR(objects.at(k).momentum) << std::endl;
//                std::cout << "first tau charge " << objects.at(n).charge << std::endl;
//                std::cout << "second tau charge " << objects.at(k).charge << std::endl;
                if(objects.at(n).momentum.DeltaR(objects.at(k).momentum) > minDeltaR) {
                    const Candidate candidate(type, objects.at(n), objects.at(k));
                    if (expectedCharge != Candidate::UnknownCharge() && candidate.charge != expectedCharge ) continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate.momentum.M(),weight);
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(),weight);
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


protected:
    tools::Timer timer;
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
