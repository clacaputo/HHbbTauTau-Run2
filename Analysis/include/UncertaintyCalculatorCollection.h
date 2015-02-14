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

#include <TFile.h>

#include "UncertaintyConfiguration.h"

namespace analysis {
namespace limits {

namespace uncertainties {

namespace names {

const std::string Btag_efficiency = "eff_b";
const std::string Btag_fake = "fake_b";
const std::string TTbar_normalization = "ttbarNorm";
const std::string ZTT_extrapolation = "extrap_ztt";
const std::string ZLL_FakeTau = "ZLL_FakeTau";
const std::string JetFakeTau = "JetFakeTau";
const std::string LeptonFakeTau = "LeptonFakeTau";
const std::string QCD = "QCDSyst";

} // namespace names


const double TTbar_CrossSection = 0.10;
const double JetFakeTau = 0.20;
const double LeptonFakeTau = 0.30;

const std::map<std::string, double> QCD = {
    { "tauTau_2jet0tag", 0.1 }, { "tauTau_2jet1tag", 0.2 }, { "tauTau_2jet2tag", 0.4 }
};

} // namespace uncertainties


class UncertaintyCalculatorCollection {
private:
    typedef std::function< UncertaintyInterval (const std::string&, const std::string&) > UncertaintyCalculator;
    typedef std::map<std::string, UncertaintyCalculator> UncertaintyCalculatorMap;

    template<typename MethodPtr>
    UncertaintyCalculator Bind(MethodPtr method_ptr) const
    {
        using namespace std::placeholders;
        return std::bind(method_ptr, this, _1, _2);
    }

public:
    UncertaintyCalculatorCollection(std::shared_ptr<TFile> _shape_file)
        : shape_file(_shape_file)
    {
        using namespace uncertainties::names;

        calculator_map[Btag_efficiency] = Bind(&UncertaintyCalculatorCollection::CalculateBtagEfficiencyUnc);
        calculator_map[Btag_fake] = Bind(&UncertaintyCalculatorCollection::CalculateBtagFakeUnc);
        calculator_map[TTbar_normalization] = Bind(&UncertaintyCalculatorCollection::CalculateTTUnc);
        calculator_map[ZTT_extrapolation] = Bind(&UncertaintyCalculatorCollection::CalculateZTTextrapUnc);
        calculator_map[ZLL_FakeTau] = Bind(&UncertaintyCalculatorCollection::CalculateZFakeTauUnc);
        calculator_map[QCD] = Bind(&UncertaintyCalculatorCollection::CalculateQCDUnc);
    }

    UncertaintyInterval Calculate(const std::string& unc_name, const std::string& full_category_name,
                                  const std::string& sample_name) const
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

    PhysicalValue CalculateZTT_SF(const std::string& full_category_name) const
    {
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

        TT_cat.AddSystematicUncertainty(uncertainties::names::TTbar_normalization, uncertainties::TTbar_CrossSection);
        TT_incl.AddSystematicUncertainty(uncertainties::names::TTbar_normalization, uncertainties::TTbar_CrossSection);

        return (DY_cat - TT_cat) / (DY_incl - TT_incl);
    }

    UncertaintyInterval CalculateTTUnc(const std::string& full_category_name, const std::string& sample_name) const
    {
        double unc;

        if (sample_name == "TT")
            unc = uncertainties::TTbar_CrossSection;
        else if (sample_name == "ZTT") {
            const PhysicalValue sf = CalculateZTT_SF(full_category_name);
            std::cout << "ZTT SF = " << sf << std::endl;
            unc = sf.GetRelativeSystematicUncertainty(uncertainties::names::TTbar_normalization);
        }
        else
            throw exception("Not found correct sample on which apply ttNorm uncertainty");

        return UncertaintyInterval(PhysicalValue(unc));
    }

    UncertaintyInterval CalculateZTTextrapUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
        const PhysicalValue sf = CalculateZTT_SF(full_category_name);
        std::cout << "ZTT SF = " << sf << std::endl;
        const double unc = sf.GetRelativeStatisticalError();
        return UncertaintyInterval(PhysicalValue(unc));
    }

    UncertaintyInterval CalculateZFakeTauUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
       const std::string hist_name_ZJ = full_category_name + "/ZJ";
       const std::string hist_name_ZL = full_category_name + "/ZL";

       auto hist_ZJ = LoadHistogram(hist_name_ZJ);
       auto hist_ZL = LoadHistogram(hist_name_ZL);

       PhysicalValue n_ZJ = Integral(*hist_ZJ, true);
       PhysicalValue n_ZL = Integral(*hist_ZL, true);

       n_ZJ.AddSystematicUncertainty(uncertainties::names::JetFakeTau, uncertainties::JetFakeTau);
       n_ZL.AddSystematicUncertainty(uncertainties::names::LeptonFakeTau, uncertainties::LeptonFakeTau);

       const PhysicalValue n_ZLL = n_ZJ + n_ZL;
       std::cout << "ZLL yield = " << n_ZLL << std::endl;

       const double unc = n_ZLL.GetRelativeFullError();
       return UncertaintyInterval(PhysicalValue(unc));
    }

    UncertaintyInterval CalculateQCDUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
        const double unc = uncertainties::QCD.at(full_category_name);
        return UncertaintyInterval(PhysicalValue(unc));
    }

private:
    std::shared_ptr<TFile> shape_file;
    UncertaintyCalculatorMap calculator_map;
};

} // namespace limits
} // namespace analysis
