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
    TH1D_ENTRY(mt_1, 20, 0, 200)

    virtual void Fill(const analysis::FlatEventInfo& eventInfo, double weight, bool fill_all, bool doESvariation = true) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight, fill_all, doESvariation);
        if(!fill_all) return;

        const ntuple::Flat& event = *eventInfo.event;
        mt_1().Fill(event.mt_1, weight);
    }
};

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                            const std::string& outputFileName, const std::string& signal_list)
          : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, ChannelId(), signal_list)
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
                                                analysis::EventRegion::OS_NotIsolated };
        return regions;
    }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;

        if(event.againstElectronLooseMVA_2 <= againstElectronLooseMVA3
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                    && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                         event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

        if(iso) return os ? EventRegion::OS_Isolated : EventRegion::SS_Isolated;
        return os ? EventRegion::OS_NotIsolated : EventRegion::SS_NotIsolated;
    }

    virtual analysis::PhysicalValue CalculateQCDScaleFactor(analysis::EventCategory eventCategory,
                                                            const std::string& hist_name) override
    {
        return CalculateQCDScaleFactor(eventCategory, hist_name, analysis::EventRegion::SS_Isolated,
                                       analysis::EventRegion::SS_NotIsolated);
    }

    virtual analysis::PhysicalValue CalculateQCDScaleFactor(analysis::EventCategory eventCategory, const std::string& hist_name,
                                                  analysis::EventRegion num_eventRegion, analysis::EventRegion den_eventRegion) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;
        std::cout<<"Ciao2"<<std::endl;
        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        const analysis::EventCategory loose_category = analysis::Medium_Loose_CategoryMap.at(eventCategory);

        auto hist_data_SS_AIso_Medium = GetHistogram(eventCategory, data.name, den_eventRegion, hist_name);
        auto hist_data_SS_Iso_Loose   = GetHistogram(loose_category, data.name, num_eventRegion, hist_name);
        auto hist_data_SS_AIso_Loose  = GetHistogram(loose_category, data.name, den_eventRegion, hist_name);

        if(!hist_data_SS_AIso_Medium || !hist_data_SS_Iso_Loose || !hist_data_SS_AIso_Loose)
            throw analysis::exception("Unable to find histograms for QCD scale factor estimation");

        TH1D& hist_SS_AIso_Medium = CloneHistogram(eventCategory, qcd.name, den_eventRegion, *hist_data_SS_AIso_Medium);
        TH1D& hist_SS_Iso_Loose = CloneHistogram(loose_category, qcd.name, num_eventRegion, *hist_data_SS_Iso_Loose);
        TH1D& hist_SS_AIso_Loose = CloneHistogram(loose_category, qcd.name, den_eventRegion, *hist_data_SS_AIso_Loose);

        SubtractBackgroundHistograms(hist_SS_AIso_Medium, eventCategory, den_eventRegion, qcd.name, true);
        SubtractBackgroundHistograms(hist_SS_Iso_Loose, loose_category, num_eventRegion, qcd.name, true);
        SubtractBackgroundHistograms(hist_SS_AIso_Loose, loose_category, den_eventRegion, qcd.name, true);

        const analysis::PhysicalValue n_SS_AIso_Medium = analysis::Integral(hist_SS_AIso_Medium, false);
        const analysis::PhysicalValue n_SS_Iso_Loose   = analysis::Integral(hist_SS_Iso_Loose, false);
        const analysis::PhysicalValue n_SS_AIso_Loose  = analysis::Integral(hist_SS_AIso_Loose, false);
        if(n_SS_AIso_Medium.value < 0 || n_SS_Iso_Loose.value < 0 || n_SS_AIso_Loose.value < 0)
            throw analysis::exception("Negative number of estimated events in QCD SF estimation for ") << eventCategory;

        const analysis::PhysicalValue bTagEff = n_SS_AIso_Medium/n_SS_AIso_Loose;
        const analysis::PhysicalValue isoEff_looseRegion = n_SS_Iso_Loose/n_SS_AIso_Loose;
        std::cout<<" BTag Eff = " << bTagEff << "isoEff = " << isoEff_looseRegion <<std::endl;

        return bTagEff*isoEff_looseRegion;
    }

    virtual void EstimateQCD(analysis::EventCategory eventCategory, const std::string& hist_name,
                             const analysis::PhysicalValue& scale_factor) override
    {
        return EstimateQCD(eventCategory,hist_name,scale_factor,
                                                 analysis::EventRegion::OS_NotIsolated);
    }

    virtual void EstimateQCD(analysis::EventCategory eventCategory, const std::string& hist_name,
                             const analysis::PhysicalValue& scale_factor, analysis::EventRegion shapeRegion) override
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        const analysis::EventCategory loose_category = analysis::Medium_Loose_CategoryMap.at(eventCategory);

        auto hist_shape_data = GetHistogram(loose_category, data.name, shapeRegion, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_Isolated, *hist_shape_data);
        SubtractBackgroundHistograms(histogram, loose_category, shapeRegion, qcd.name);
        histogram.Scale(scale_factor.value);

    }

    virtual PhysicalValueMap CalculateWjetsScaleFactors(analysis::EventCategory eventCategory,
                                                                   const std::string& hist_name) override
    {

        PhysicalValueMap valueMap;

        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(analysis::DataCategoryType::WJets_MC);
        std::cout<<"Event Category ------>>>>>>>  "<<eventCategory<<std::endl;
        if(analysis::TwoJetsEventCategories_LooseBjets.count(eventCategory)) {
             static const analysis::PhysicalValue v(1, 0.001);
             for (analysis::EventRegion eventRegion : analysis::QcdRegions){
                 const analysis::EventRegion loose_eventRegion = analysis::MediumBtag_LooseBtag_RegionMap.at(eventRegion);
                  valueMap[loose_eventRegion] = v;
                }
        } else {

             const analysis::EventCategory loose_category = analysis::Medium_Loose_CategoryMap.at(eventCategory);
             for (analysis::EventRegion eventRegion : analysis::QcdRegions){
                 //analysis::PhysicalValue region_SF(0,0);
                 TH1D* hist_loose=nullptr;
                 TH1D* hist_medium=nullptr;

                 const analysis::EventRegion loose_eventRegion = analysis::MediumBtag_LooseBtag_RegionMap.at(eventRegion);
                 for (const analysis::DataCategory* wjets_category : wjets_mc_categories){

                     auto loose_hist_mc = GetHistogram(loose_category, wjets_category->name, eventRegion, hist_name);
                     auto medium_hist_mc = GetHistogram(eventCategory, wjets_category->name, eventRegion, hist_name);
                     std::cout<<"loose_hist_mc ----> "<<loose_hist_mc<<" medium_hist_mc ----> "<<medium_hist_mc<<std::endl;
                     if (!loose_hist_mc && !medium_hist_mc) continue;
                     if (!hist_loose && !hist_medium) {

                        hist_loose  = &CloneHistogram(loose_category, wjets_category->name, loose_eventRegion, *loose_hist_mc);
                        hist_medium = &CloneHistogram(eventCategory, wjets_category->name, loose_eventRegion, *medium_hist_mc);
                     } else{
                         hist_loose->Add(loose_hist_mc);
                         hist_medium->Add(medium_hist_mc);
                     }
                     std::cout<<"Event Region:  "<<eventRegion<<"WJets Category:   "<<wjets_category->name<<"  Integral Loose:  "<<analysis::Integral(*hist_loose,false)
                             <<" Integral Medium: "<<analysis::Integral(*hist_medium,false)<<std::endl;
                     //region_SF = region_SF + analysis::Integral(hist_loose,false)/analysis::Integral(hist_medium,false);
                 }
                valueMap[eventRegion] = analysis::Integral(*hist_medium,false)/analysis::Integral(*hist_loose,false);
             }
         }

        return valueMap;
    }

};
