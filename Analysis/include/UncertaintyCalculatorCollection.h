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
    typedef std::function< UncertaintyInterval (const std::string&, const std::string&) > UncertaintyCalculator;
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
                                    std::shared_ptr<TFile> _shape_file)
        : uncertainties(&_uncertainties), shape_file(_shape_file)
    {
        using namespace uncertainty_names;

        Bind(Btag_efficiency, &UncertaintyCalculatorCollection::CalculateBtagEfficiencyUnc);
        Bind(Btag_fake, &UncertaintyCalculatorCollection::CalculateBtagFakeUnc);
        Bind(TTbar_normalization, &UncertaintyCalculatorCollection::CalculateTTUnc);
        Bind(ZTT_extrapolation, &UncertaintyCalculatorCollection::CalculateZTTextrapUnc);
        Bind(ZLL_FakeTau, &UncertaintyCalculatorCollection::CalculateZFakeTauUnc);
        Bind(QCD, &UncertaintyCalculatorCollection::CalculateQCDUnc);
    }

    UncertaintyInterval Calculate(const std::string& unc_name, const std::string& full_category_name,
                                  const std::string& sample_name)
    {
        if(!calculator_map.count(unc_name))
            throw exception("Calculator for uncertainty '") << unc_name << "' not found.";
        return calculator_map.at(unc_name)(full_category_name, sample_name);
    }

    static PhysicalValue CombineUpDownUncertainties(const UncertaintyInterval& unc)
    {
        return ( unc.up + (PhysicalValue::Two - unc.down) ) / PhysicalValue::Two;
    }

private:
    std::shared_ptr<TH1D> LoadHistogram(const std::string& hist_name) const
    {
        std::shared_ptr<TH1D> hist(static_cast<TH1D*>(shape_file->Get(hist_name.c_str())));
        if(!hist)
            throw exception("Histogram '") << hist_name << "' not found.";
        return hist;
    }

    void AddUncertainty(PhysicalValue& physicalValue, const std::string& unc_name) const
    {
        physicalValue.AddSystematicUncertainty(unc_name, uncertainties->Get(unc_name).value - 1);
    }

    UncertaintyInterval CalculateBtagEfficiencyUnc(const std::string& full_category_name,
                                                   const std::string& sample_name) const
    {
        const std::string hist_name = full_category_name + "/" + sample_name;
        const std::string hist_name_scale_prefix = hist_name + "_CMS_scale_btagEff_8TeV";
        const std::string hist_up_name = hist_name_scale_prefix + "Up";
        const std::string hist_down_name = hist_name_scale_prefix + "Down";

        auto hist = LoadHistogram(hist_name);
        auto hist_up = LoadHistogram(hist_up_name);
        auto hist_down = LoadHistogram(hist_down_name);

        const PhysicalValue n_central = Integral(*hist, true);
        const PhysicalValue n_up = Integral(*hist_up, true);
        const PhysicalValue n_down = Integral(*hist_down, true);

        return UncertaintyInterval(n_down/n_central, n_up/n_central);
    }

    UncertaintyInterval CalculateBtagFakeUnc(const std::string& full_category_name,
                                             const std::string& sample_name) const
    {
        const std::string hist_name = full_category_name + "/" + sample_name;
        const std::string hist_name_scale_prefix = hist_name + "_CMS_scale_btagFake_8TeV";
        const std::string hist_up_name = hist_name_scale_prefix + "Up";
        const std::string hist_down_name = hist_name_scale_prefix + "Down";

        auto hist = LoadHistogram(hist_name);
        auto hist_up = LoadHistogram(hist_up_name);
        auto hist_down = LoadHistogram(hist_down_name);

        const PhysicalValue n_central = Integral(*hist, true);
        const PhysicalValue n_up = Integral(*hist_up, true);
        const PhysicalValue n_down = Integral(*hist_down, true);

        return UncertaintyInterval(n_down/n_central, n_up/n_central);
    }

    const PhysicalValue& GetZTT_SF(const std::string& full_category_name)
    {
        if(!ztt_sf_map.count(full_category_name)) {
            const std::string hist_name_DYemb_cat = full_category_name + "/DY_emb";
            const std::string hist_name_TTemb_cat = full_category_name + "/TT_emb";
            const std::string hist_name_DYemb_incl = "tauTau_inclusive/DY_emb";
            const std::string hist_name_TTemb_incl = "tauTau_inclusive/TT_emb";

            auto hist_DYemb_cat = LoadHistogram(hist_name_DYemb_cat);
            auto hist_TTemb_cat = LoadHistogram(hist_name_TTemb_cat);
            auto hist_DYemb_incl = LoadHistogram(hist_name_DYemb_incl);
            auto hist_TTemb_incl = LoadHistogram(hist_name_TTemb_incl);

            PhysicalValue DY_cat = Integral(*hist_DYemb_cat, true);
            PhysicalValue TT_cat = Integral(*hist_TTemb_cat, true);
            PhysicalValue DY_incl = Integral(*hist_DYemb_incl, true);
            PhysicalValue TT_incl = Integral(*hist_TTemb_incl, true);

            AddUncertainty(TT_cat, uncertainty_names::TTbar_normalization);
            AddUncertainty(TT_incl, uncertainty_names::TTbar_normalization);

            const PhysicalValue sf = (DY_cat - TT_cat) / (DY_incl - TT_incl);
            std::cout << "ZTT SF = " << sf << std::endl;
            ztt_sf_map[full_category_name] = sf;
        }
        return ztt_sf_map.at(full_category_name);
    }

    UncertaintyInterval CalculateTTUnc(const std::string& full_category_name, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::TTbar_normalization << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(full_category_name);
        const double unc = sf.GetRelativeSystematicUncertainty(uncertainty_names::TTbar_normalization);
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZTTextrapUnc(const std::string& full_category_name, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZTT_extrapolation << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(full_category_name);
        const double unc = sf.GetRelativeStatisticalError();
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZFakeTauUnc(const std::string& full_category_name,
                                             const std::string& sample_name) const
    {
        if(sample_name != "ZLL")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZLL_FakeTau << " uncertainty calculator.";

       const std::string hist_name_ZJ = full_category_name + "/ZJ";
       const std::string hist_name_ZL = full_category_name + "/ZL";

       auto hist_ZJ = LoadHistogram(hist_name_ZJ);
       auto hist_ZL = LoadHistogram(hist_name_ZL);

       PhysicalValue n_ZJ = Integral(*hist_ZJ, true);
       PhysicalValue n_ZL = Integral(*hist_ZL, true);

       AddUncertainty(n_ZJ, uncertainty_names::JetFakeTau);
       AddUncertainty(n_ZL, uncertainty_names::LeptonFakeTau);

       const PhysicalValue n_ZLL = n_ZJ + n_ZL;
       std::cout << "ZLL yield = " << n_ZLL << std::endl;

       const double unc = n_ZLL.GetRelativeFullError();
       return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateQCDUnc(const std::string& full_category_name, const std::string& sample_name) const
    {
        if(sample_name != "QCD")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::QCD << " uncertainty calculator.";

        const std::map<std::string, double> QCD_values = {
            { "tauTau_2jet0tag", 0.1 }, { "tauTau_2jet1tag", 0.2 }, { "tauTau_2jet2tag", 0.4 }
        };

        const double unc = QCD_values.at(full_category_name);
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

private:
    const UncertaintyDescriptorCollection* uncertainties;
    std::shared_ptr<TFile> shape_file;
    UncertaintyCalculatorMap calculator_map;
    std::map<std::string, PhysicalValue> ztt_sf_map;
};

} // namespace limits
} // namespace analysis
