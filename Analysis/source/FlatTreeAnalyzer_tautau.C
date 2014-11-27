/*!
 * \file FlatTreeAnalyzer_tautau.C
 * \brief Analyze flat-tree for tau-tau channel for HHbbtautau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-03 created
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

#include "Analysis/include/BaseFlatTreeAnalyzer.h"

class FlatAnalyzerData_tautau : public analysis::FlatAnalyzerData {
public:
    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(iso_tau1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1)
    TH1D_ENTRY_EX(iso_tau2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1)
    TH1D_ENTRY_EX(H_tt_charge, 8, -4, 4, "H_{#tau#tau} Charge", "Events", false, 1.1)

    virtual void Fill(const analysis::FlatEventInfo& eventInfo, double weight, bool fill_all,
                      analysis::EventEnergyScale eventEnergyScale) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight, fill_all, eventEnergyScale);
        if(!fill_all) return;
        if (eventEnergyScale != analysis::EventEnergyScale::Central) return;
        const ntuple::Flat& event = *eventInfo.event;
        mt_1().Fill(event.mt_1, weight);
        iso_tau1().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1,weight);
        iso_tau2().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2,weight);
        int diTau_charge = event.q_1*event.q_2;
        H_tt_charge().Fill(diTau_charge,weight);
    }
};

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& _inputPath,
                            const std::string& outputFileName, const std::string& signal_list)
          : BaseFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()), _inputPath,
                                 outputFileName)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::TauTau; }

    virtual std::shared_ptr<analysis::FlatAnalyzerData> MakeAnaData() override
    {
        return std::shared_ptr<FlatAnalyzerData_tautau>(new FlatAnalyzerData_tautau());
    }

    virtual const analysis::EventRegionSet& EssentialEventRegions() override
    {
        static const analysis::EventRegionSet regions = { analysis::EventRegion::OS_Isolated,
                                                          analysis::EventRegion::OS_AntiIsolated };
        return regions;
    }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory /*eventCategory*/) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;

        if(!event.againstElectronLooseMVA_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                    && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                         event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

        if(iso) return os ? EventRegion::OS_Isolated : EventRegion::SS_Isolated;
        return os ? EventRegion::OS_AntiIsolated : EventRegion::SS_AntiIsolated;
    }

    virtual analysis::PhysicalValue CalculateQCDScaleFactor(const std::string& hist_name,
                                                            analysis::EventCategory eventCategory) override
    {
        analysis::EventCategory refEventCategory = eventCategory;
        if(CategoriesToRelax().count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        const auto iso_antiIso_sf = CalculateQCDScaleFactorEx(hist_name, refEventCategory, refEventCategory,
                                                              analysis::EventRegion::SS_Isolated,
                                                              analysis::EventRegion::SS_AntiIsolated);
        std::cout << "iso_antiIso_sf: " << iso_antiIso_sf << "\n";
        if(refEventCategory == eventCategory)
            return iso_antiIso_sf;

        const auto medium_loose_sf = CalculateQCDScaleFactorEx(hist_name, eventCategory, refEventCategory,
                                                               analysis::EventRegion::SS_AntiIsolated,
                                                               analysis::EventRegion::SS_AntiIsolated);
        std::cout << "medium_loose_sf: " << medium_loose_sf << "\n";
        return iso_antiIso_sf * medium_loose_sf;
    }

    virtual void EstimateQCD(const std::string& hist_name, analysis::EventCategory eventCategory,
                             const analysis::PhysicalValue& scale_factor) override
    {
        analysis::EventCategory refEventCategory = eventCategory;
        if(CategoriesToRelax().count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        return EstimateQCDEx(hist_name, eventCategory, refEventCategory, analysis::EventRegion::OS_AntiIsolated,
                             scale_factor);
    }

    virtual PhysicalValueMap CalculateWjetsScaleFactors(analysis::EventCategory eventCategory,
                                                        const std::string& hist_name) override
    {
        static const analysis::PhysicalValue one(1, 0.001);

        PhysicalValueMap valueMap;

        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(analysis::DataCategoryType::WJets_MC);

        analysis::EventCategory refEventCategory = eventCategory;
        if(CategoriesToRelax().count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        for (analysis::EventRegion eventRegion : analysis::QcdRegions) {
            analysis::PhysicalValue sf = one;
            if(eventCategory != refEventCategory) {
                try {
                    const auto n_medium = CalculateFullIntegral(eventCategory, eventRegion, hist_name, wjets_mc_categories);
                    const auto n_ref = CalculateFullIntegral(refEventCategory, eventRegion, hist_name, wjets_mc_categories);
                    sf = n_medium / n_ref;
                } catch(analysis::exception& ex) {
                    std::cerr << ex.message() << " Using default value " << one << "." << std::endl;
                }
            }
            valueMap[eventRegion] = sf;
        }

        return valueMap;
    }
};
