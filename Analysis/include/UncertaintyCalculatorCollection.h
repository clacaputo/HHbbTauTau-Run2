/*!
 * \file UncertaintyCalculatorCollection.h
 * \brief Collection of methods to calculate uncertainties.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2015-02-04 created
 *
 * Copyright 2015 Konstantin Androsov <konstantin.androsov@gmail.com>,
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

#include "UncertaintyConfiguration.h"
#include "FlatAnalyzerDataCollection.h"

namespace analysis {
namespace limits {

namespace uncertainty_names {

const std::string Btag_efficiency = "eff_b";
const std::string Btag_fake = "fake_b";
const std::string TTbar_normalization = "ttbarNorm";
const std::string ZTT_extrapolation = "extrap_ztt";
const std::string ZLL_FakeTau = "ZLL_FakeTau";
const std::string JetFakeTau = "JetFakeTau";
const std::string LeptonFakeTau = "LeptonFakeTau";
const std::string QCD = "QCDSyst";

} // namespace uncertainty_names

class UncertaintyCalculatorCollection {
private:
    typedef std::function< UncertaintyInterval (EventCategory, const std::string&) > UncertaintyCalculator;
    typedef std::map<std::string, UncertaintyCalculator> UncertaintyCalculatorMap;

    template<typename MethodPtr>
    void Bind(const std::string& name, MethodPtr method_ptr)
    {
        using namespace std::placeholders;
        calculator_map[name] = std::bind(method_ptr, this, _1, _2);
    }

    static double DefaultPrecision() { return 0.001; }

public:
    UncertaintyCalculatorCollection(const UncertaintyDescriptorCollection& _uncertainties,
                                    const DataCategoryCollection& _dataCategories,
                                    const FlatAnalyzerDataCollectionReader& _reader,
                                    EventSubCategory _eventSubCategory, const std::string& _referenceHistName)
        : uncertainties(&_uncertainties), dataCategories(&_dataCategories), reader(&_reader),
          eventSubCategory(_eventSubCategory), referenceHistName(_referenceHistName)
    {
        using namespace uncertainty_names;

        Bind(Btag_efficiency, &UncertaintyCalculatorCollection::CalculateBtagEfficiencyUnc);
        Bind(Btag_fake, &UncertaintyCalculatorCollection::CalculateBtagFakeUnc);
        Bind(TTbar_normalization, &UncertaintyCalculatorCollection::CalculateTTUnc);
        Bind(ZTT_extrapolation, &UncertaintyCalculatorCollection::CalculateZTTextrapUnc);
        Bind(ZLL_FakeTau, &UncertaintyCalculatorCollection::CalculateZFakeTauUnc);
        Bind(QCD, &UncertaintyCalculatorCollection::CalculateQCDUnc);
    }

    UncertaintyInterval Calculate(const std::string& unc_name, EventCategory event_category,
                                  const std::string& sample_name)
    {
        if(!calculator_map.count(unc_name))
            throw exception("Calculator for uncertainty '") << unc_name << "' not found.";
        return calculator_map.at(unc_name)(event_category, sample_name);
    }

    static PhysicalValue CombineUpDownUncertainties(const UncertaintyInterval& unc)
    {
        return ( unc.up + (PhysicalValue::Two - unc.down) ) / PhysicalValue::Two;
    }

private:
    PhysicalValue GetYield(EventCategory eventCategory, const std::string& dataCategoryName,
                           EventRegion eventRegion = EventRegion::OS_Isolated,
                           EventEnergyScale eventEnergyScale = EventEnergyScale::Central) const
    {
        const FlatAnalyzerDataId id(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
        auto hist = reader->GetHistogram<TH1D>(id, referenceHistName);
        if(hist)
            return Integral(*hist, true);
        std::cout << "Histogram '" << referenceHistName << "' in " << id << " not found. Considering zero yield.\n";
        return PhysicalValue::Zero;
    }

    const std::string& GetDataCategoryName(const std::string& datacard) const
    {
        return dataCategories->FindCategoryForDatacard(datacard).name;
    }

    void AddUncertainty(PhysicalValue& physicalValue, const std::string& unc_name) const
    {
        physicalValue.AddSystematicUncertainty(unc_name, uncertainties->Get(unc_name).value - 1);
    }

    UncertaintyInterval CalculateEnergyScaleRelatedUncertainty(EventCategory eventCategory,
                                                               const std::string& dataCategoryName,
                                                               EventEnergyScale up_energy_scale,
                                                               EventEnergyScale down_energy_scale) const
    {
        const PhysicalValue n_central = GetYield(eventCategory, dataCategoryName,
                                                 EventRegion::OS_Isolated, EventEnergyScale::Central);
        const PhysicalValue n_up = GetYield(eventCategory, dataCategoryName,
                                            EventRegion::OS_Isolated, up_energy_scale);
        const PhysicalValue n_down = GetYield(eventCategory, dataCategoryName,
                                              EventRegion::OS_Isolated, down_energy_scale);

        return UncertaintyInterval(n_down / n_central, n_up / n_central);
    }


    UncertaintyInterval CalculateBtagEfficiencyUnc(EventCategory eventCategory,
                                                   const std::string& sample_name) const
    {
        const auto& dataCategoryName = GetDataCategoryName(sample_name);
        return CalculateEnergyScaleRelatedUncertainty(eventCategory, dataCategoryName,
                                                      EventEnergyScale::BtagEfficiencyUp,
                                                      EventEnergyScale::BtagEfficiencyDown);
    }

    UncertaintyInterval CalculateBtagFakeUnc(EventCategory eventCategory,
                                             const std::string& sample_name) const
    {
        const auto& dataCategoryName = GetDataCategoryName(sample_name);
        return CalculateEnergyScaleRelatedUncertainty(eventCategory, dataCategoryName,
                                                      EventEnergyScale::BtagFakeUp,
                                                      EventEnergyScale::BtagFakeDown);
    }

    const PhysicalValue& GetZTT_SF(EventCategory eventCategory)
    {
        if(!ztt_sf_map.count(eventCategory)) {
            const std::string& DY_emb_name = dataCategories->GetUniqueCategory(DataCategoryType::Embedded).name;
            const std::string& TT_emb_name = dataCategories->GetUniqueCategory(DataCategoryType::TT_Embedded).name;

            PhysicalValue DY_cat = GetYield(eventCategory, DY_emb_name);
            PhysicalValue TT_cat = GetYield(eventCategory, TT_emb_name);
            PhysicalValue DY_incl = GetYield(EventCategory::Inclusive, DY_emb_name);
            PhysicalValue TT_incl = GetYield(EventCategory::Inclusive, TT_emb_name);

            AddUncertainty(TT_cat, uncertainty_names::TTbar_normalization);
            AddUncertainty(TT_incl, uncertainty_names::TTbar_normalization);

            const PhysicalValue sf = (DY_cat - TT_cat) / (DY_incl - TT_incl);
            std::cout << "ZTT SF = " << sf << std::endl;
            ztt_sf_map[eventCategory] = sf;
        }
        return ztt_sf_map.at(eventCategory);
    }

    UncertaintyInterval CalculateTTUnc(EventCategory eventCategory, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::TTbar_normalization << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(eventCategory);
        const double unc = sf.GetRelativeSystematicUncertainty(uncertainty_names::TTbar_normalization);
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZTTextrapUnc(EventCategory eventCategory, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZTT_extrapolation << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(eventCategory);
        const double unc = sf.GetRelativeStatisticalError();
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZFakeTauUnc(EventCategory eventCategory, const std::string& sample_name) const
    {
        if(sample_name != "ZLL")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZLL_FakeTau << " uncertainty calculator.";

        const std::string ZJ_name = dataCategories->GetUniqueCategory(DataCategoryType::ZJ_MC).name;
        const std::string ZL_name = dataCategories->GetUniqueCategory(DataCategoryType::ZL_MC).name;

        PhysicalValue n_ZJ = GetYield(eventCategory, ZJ_name);
        PhysicalValue n_ZL = GetYield(eventCategory, ZL_name);

        AddUncertainty(n_ZJ, uncertainty_names::JetFakeTau);
        AddUncertainty(n_ZL, uncertainty_names::LeptonFakeTau);

        const PhysicalValue n_ZLL = n_ZJ + n_ZL;
        std::cout << "ZLL yield = " << n_ZLL << std::endl;

        const double unc = n_ZLL.GetRelativeFullError();
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateQCDUnc(EventCategory eventCategory, const std::string& sample_name) const
    {
        if(sample_name != "QCD")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::QCD << " uncertainty calculator.";

        const std::map<EventCategory, double> QCD_values = {
            { EventCategory::TwoJets_ZeroBtag, 0.1 },
            { EventCategory::TwoJets_OneBtag, 0.2 },
            { EventCategory::TwoJets_TwoBtag, 0.4 }
        };

        const double unc = QCD_values.at(eventCategory);
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

private:
    const UncertaintyDescriptorCollection* uncertainties;
    const DataCategoryCollection* dataCategories;
    const FlatAnalyzerDataCollectionReader* reader;
    EventSubCategory eventSubCategory;
    std::string referenceHistName;
    UncertaintyCalculatorMap calculator_map;
    std::map<EventCategory, PhysicalValue> ztt_sf_map;
};

} // namespace limits
} // namespace analysis
