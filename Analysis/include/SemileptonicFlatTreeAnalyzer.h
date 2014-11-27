/*!
 * \file SemileptonicFlatTreeAnalyzer.h
 * \brief Definition of the base class for semi-leptonic flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-11-20 created
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

#include "BaseFlatTreeAnalyzer.h"

namespace analysis {

class SemileptonicFlatTreeAnalyzer : public BaseFlatTreeAnalyzer {
public:
    SemileptonicFlatTreeAnalyzer(const DataCategoryCollection& _dataCategoryCollection, const std::string& _inputPath,
                                 const std::string& _outputFileName)
         : BaseFlatTreeAnalyzer(_dataCategoryCollection, _inputPath, _outputFileName)
    {
    }

protected:
    bool IsHighMtRegion(const ntuple::Flat& event, analysis::EventCategory eventCategory)
    {
        using namespace cuts;
        if (eventCategory == analysis::EventCategory::TwoJets_TwoBtag)
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion_low &&
                    event.mt_1 < WjetsBackgroundEstimation::HighMtRegion_high;
        else
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion;
    }

    virtual PhysicalValue CalculateQCDScaleFactor(const std::string& hist_name, EventCategory eventCategory) override
    {
        static const PhysicalValue sf(1.06, 0.001);
        return sf;
    }

    virtual analysis::PhysicalValue CalculateQCDScaleFactor(const std::string& hist_name,
                                                            analysis::EventCategory eventCategory) override
    {
        static const PhysicalValue sf(1.06, 0.001);
        static const analysis::EventCategorySet categories= {analysis::EventCategory::TwoJets_TwoBtag};
        analysis::EventCategory refEventCategory = eventCategory;
        if(categories.count(eventCategory))
            //refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);
            refEventCategory = analysis::EventCategory::TwoJets_Inclusive;

        const analysis::PhysicalValue yield_SSIso =
                CalculateYieldsForQCD(hist_name,refEventCategory,analysis::EventRegion::SS_Isolated);

        std::cout << "yield_ssIso: " << yield_SSIso << "\n";
        if(refEventCategory == eventCategory)
            return sf*yield_SSIso;

        const analysis::PhysicalValue yield_SSAntiIso_Medium =
                CalculateYieldsForQCD(hist_name,eventCategory,analysis::EventRegion::SS_AntiIsolated);

        const analysis::PhysicalValue yield_SSAntiIso_Loose =
                CalculateYieldsForQCD(hist_name,refEventCategory,analysis::EventRegion::SS_AntiIsolated);

        const auto medium_loose_sf = yield_SSAntiIso_Medium/yield_SSAntiIso_Loose;
        std::cout << "medium_loose_sf: " << medium_loose_sf << "\n";

        return sf*yield_SSIso * medium_loose_sf;
    }


    virtual void EstimateQCD(const std::string& hist_name, analysis::EventCategory eventCategory,
                             const analysis::PhysicalValue& scale_factor) override
    {
        static const analysis::EventCategorySet categories= {analysis::EventCategory::TwoJets_OneBtag,
                                                             analysis::EventCategory::TwoJets_TwoBtag};
        analysis::EventCategory refEventCategory = eventCategory;
        if(categories.count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        auto hist_shape_data = GetHistogram(refEventCategory, data.name, analysis::EventRegion::SS_AntiIsolated, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, analysis::EventRegion::OS_Isolated, *hist_shape_data);
        //here is missing subtraction for inclusive part
        analysis::RenormalizeHistogram(histogram,scale_factor,true);

    }

    virtual PhysicalValueMap CalculateWjetsYields(EventCategory eventCategory,
                                                        const std::string& hist_name) override
    {
        PhysicalValueMap valueMap;
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        static const EventCategorySet categoriesToRelax= {EventCategory::TwoJets_TwoBtag};
        analysis::EventCategory refEventCategory = eventCategory;
        if(categoriesToRelax.count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        for (const auto& eventRegion : HighMt_LowMt_RegionMap){

            const PhysicalValue bkg_yield =
                    CalculateBackgroundIntegral(hist_name,eventCategory,eventRegion.first,wjets.name);

            auto hist_data = GetHistogram(eventCategory, data.name, eventRegion.first, hist_name);
            if(!hist_data)
                throw exception("Unable to find data histograms for Wjet scale factors estimation");
            const PhysicalValue n_HighMt = Integral(*hist_data, true) - bkg_yield;
            if(n_HighMt.value < 0)
                throw exception("Negative number of estimated events in Wjets SF estimation for ")
                        << eventCategory << " " << eventRegion.first << ".";

            PhysicalValue n_HighMt_mc = CalculateFullIntegral(refEventCategory,eventRegion.first,hist_name,wjets_mc_categories);
            PhysicalValue n_LowMt_mc = CalculateFullIntegral(refEventCategory,eventRegion.second,hist_name,wjets_mc_categories);

            PhysicalValue ratio_LowToHighMt = n_LowMt_mc / n_HighMt_mc;

            valueMap[eventRegion.second] = n_HighMt * ratio_LowToHighMt;
        }

        return valueMap;
    }

};

} // namespace analysis
