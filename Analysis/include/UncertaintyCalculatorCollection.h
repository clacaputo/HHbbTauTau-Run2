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

class UncertaintyCalculatorCollection {
private:
    typedef std::function< UncertaintyInterval (const std::string&, const std::string&) > UncertaintyCalculator;
    typedef std::map<std::string, UncertaintyCalculator> UncertaintyCalculatorMap;

public:
    UncertaintyCalculatorCollection(std::shared_ptr<TFile> _shape_file)
        : shape_file(_shape_file)
    {
        using namespace std::placeholders;

        calculator_map["eff_b"] = std::bind(&UncertaintyCalculatorCollection::CalculateBtagEfficiencyUnc, this, _1, _2);
        calculator_map["fake_b"] = std::bind(&UncertaintyCalculatorCollection::CalculateBtagFakeUnc, this, _1, _2);
        calculator_map["ttbarNorm"] = std::bind(&UncertaintyCalculatorCollection::CalculateTTUnc, this, _1, _2);
        calculator_map["extrap_ztt"] = std::bind(&UncertaintyCalculatorCollection::CalculateZTTextrapUnc, this, _1, _2);
        calculator_map["ZLL_FakeTau"] = std::bind(&UncertaintyCalculatorCollection::CalculateZFakeTauUnc, this, _1, _2);
        //calculator_map["ZJetFakeTau"] = std::bind(&UncertaintyCalculatorCollection::CalculateZJetFakeTauUnc, this, _1, _2);
        //calculator_map["ZLeptonFakeTau"] = std::bind(&UncertaintyCalculatorCollection::CalculateZLeptonFakeTauUnc, this, _1, _2);
        calculator_map["QCDSyst"] = std::bind(&UncertaintyCalculatorCollection::CalculateQCDUnc, this, _1, _2);

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
        static const PhysicalValue two(2, 0);
        return ( unc.up + (two - unc.down) ) / two;
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

    UncertaintyInterval CalculateTTUnc(const std::string& full_category_name,
                                                   const std::string& sample_name) const
    {
        static const PhysicalValue one(1, 0);

        PhysicalValue total_unc;

        if (sample_name == "TT")
            total_unc.value = 0.10;
        else if (sample_name == "ZTT"){
            const std::string hist_name_DYemb_cat = full_category_name + "/DY_emb";
            const std::string hist_name_TTemb_cat = full_category_name + "/TT_emb";
            const std::string hist_name_DYemb_incl = "tauTau_inclusive/DY_emb";
            const std::string hist_name_TTemb_incl = "tauTau_inclusive/TT_emb";

            auto hist_DYemb_cat = LoadHistogram(hist_name_DYemb_cat);
            auto hist_TTemb_cat = LoadHistogram(hist_name_TTemb_cat);
            auto hist_DYemb_incl = LoadHistogram(hist_name_DYemb_incl);
            auto hist_TTemb_incl = LoadHistogram(hist_name_TTemb_incl);

            const PhysicalValue x_c = Integral(*hist_TTemb_cat, true) / Integral(*hist_DYemb_cat, true);
            const PhysicalValue x_i = Integral(*hist_TTemb_incl, true) / Integral(*hist_DYemb_incl, true);
            const PhysicalValue ratio = (x_i - x_c)/((one - x_i)*(one - x_c));
            total_unc = ratio.Scale(0.1);
        }
        else
            throw exception("Not found correct sample on which apply ttNorm uncertainty");

        return UncertaintyInterval(one - total_unc, one + total_unc);
    }

    UncertaintyInterval CalculateZTTextrapUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
       static const PhysicalValue one(1, 0);

       const std::string hist_name_DYemb_cat = full_category_name + "/DY_emb";
       const std::string hist_name_TTemb_cat = full_category_name + "/TT_emb";
       const std::string hist_name_DYemb_incl = "tauTau_inclusive/DY_emb";
       const std::string hist_name_TTemb_incl = "tauTau_inclusive/TT_emb";

       auto hist_DYemb_cat = LoadHistogram(hist_name_DYemb_cat);
       auto hist_TTemb_cat = LoadHistogram(hist_name_TTemb_cat);
       auto hist_DYemb_incl = LoadHistogram(hist_name_DYemb_incl);
       auto hist_TTemb_incl = LoadHistogram(hist_name_TTemb_incl);

       const PhysicalValue DY_cat = Integral(*hist_DYemb_cat, true);
       const PhysicalValue TT_cat = Integral(*hist_TTemb_cat, true);
       const PhysicalValue DY_incl = Integral(*hist_DYemb_incl, true);
       const PhysicalValue TT_incl = Integral(*hist_TTemb_incl, true);
//       std::cout << "DY_cat: " << DY_cat << std::endl;
//       std::cout << "TT_cat: " << TT_cat << std::endl;
//       std::cout << "DY_incl: " << DY_incl << std::endl;
//       std::cout << "TT_incl: " << TT_incl << std::endl;

       const PhysicalValue scale_factor = (DY_cat - TT_cat)/(DY_incl - TT_incl);
       std::cout << "scale factor: " << scale_factor.value << ", error: " << scale_factor.error <<
                    ", ratio: " << scale_factor.error/scale_factor.value << std::endl;
       const PhysicalValue total_unc(scale_factor.error/scale_factor.value,0);
       std::cout << "total unc: " << total_unc << std::endl;

       const PhysicalValue DY_cat_mod(DY_cat.value, 0);
       const PhysicalValue TT_cat_mod(TT_cat.value, 0);
       const PhysicalValue DY_incl_mod(DY_incl.value, 0);
       const PhysicalValue TT_incl_mod(TT_incl.value, 0);

       const PhysicalValue sf_DY_only = (DY_cat - TT_cat_mod)/(DY_incl - TT_incl_mod);
       const PhysicalValue dy_unc(sf_DY_only.error/sf_DY_only.value,0);
       std::cout << "dy unc: " << dy_unc << ", contrib: " << dy_unc / total_unc << std::endl;

       const PhysicalValue sf_TT_only = (DY_cat_mod - TT_cat)/(DY_incl_mod - TT_incl);
       const PhysicalValue tt_unc(sf_TT_only.error/sf_TT_only.value,0);
       std::cout << "tt unc: " << tt_unc << ", contrib: " << tt_unc / total_unc << std::endl;


       return UncertaintyInterval(one - total_unc, one + total_unc);
    }

    UncertaintyInterval CalculateZFakeTauUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
       static const PhysicalValue one(1, 0);

       const std::string hist_name_ZJ = full_category_name + "/ZJ";
       const std::string hist_name_ZL = full_category_name + "/ZL";
       const std::string hist_name_ZLL = full_category_name + "/" + sample_name;

       auto hist_ZJ = LoadHistogram(hist_name_ZJ);
       auto hist_ZL = LoadHistogram(hist_name_ZL);
       auto hist_ZLL = LoadHistogram(hist_name_ZLL);

       const PhysicalValue n_ZJ = Integral(*hist_ZJ, true);
       const PhysicalValue n_ZL = Integral(*hist_ZL, true);
       const PhysicalValue n_ZLL = Integral(*hist_ZLL, true);

       const double uncertainty = std::sqrt( (sqr(n_ZJ.error) + sqr(0.20 * n_ZJ.value) +
                                              sqr(n_ZL.error) + sqr(0.30 * n_ZL.value)) / sqr(n_ZLL.value));

       const PhysicalValue total_unc(uncertainty,0);
       return UncertaintyInterval(one - total_unc, one + total_unc);
    }

//    UncertaintyInterval CalculateZJetFakeTauUnc(const std::string& full_category_name,
//                                              const std::string& sample_name) const
//    {
//       static const PhysicalValue one(1, 0);

//       const std::string hist_name_ZJ = full_category_name + "/ZJ";

//       auto hist_ZJ = LoadHistogram(hist_name_ZJ);

//       const PhysicalValue n_ZJ = Integral(*hist_ZJ, true);

//       const double uncertainty = std::sqrt( sqr(n_ZJ.error/n_ZJ.value) + sqr(0.20) );
//       const PhysicalValue total_unc(uncertainty,0);

//       return UncertaintyInterval(one - total_unc, one + total_unc);
//    }

//    UncertaintyInterval CalculateZLeptonFakeTauUnc(const std::string& full_category_name,
//                                              const std::string& sample_name) const
//    {
//        static const PhysicalValue one(1, 0);

//        const std::string hist_name_ZL = full_category_name + "/ZL";

//        auto hist_ZL = LoadHistogram(hist_name_ZL);

//        const PhysicalValue n_ZL = Integral(*hist_ZL, true);
//        std::cout << "ZL: " << n_ZL.value << std::endl;
//        static const double tiny_value = 1e-7;
//        if (n_ZL.value < tiny_value)
//            return UncertaintyInterval(one, one);

//        const double uncertainty = std::sqrt( sqr(n_ZL.error/n_ZL.value) + sqr(0.30) );
//        const PhysicalValue total_unc(uncertainty,0);

//        return UncertaintyInterval(one - total_unc, one + total_unc);
//    }

    UncertaintyInterval CalculateQCDUnc(const std::string& full_category_name,
                                              const std::string& sample_name) const
    {
        static const PhysicalValue one(1, 0);

        PhysicalValue total_unc;
        if (full_category_name == "tauTau_2jet0tag")
            total_unc.value = 0.10;
        else if (full_category_name == "tauTau_2jet1tag")
            total_unc.value = 0.20;
        else if (full_category_name == "tauTau_2jet2tag")
            total_unc.value = 0.40;

        return UncertaintyInterval(one - total_unc, one + total_unc);
    }

private:
    std::shared_ptr<TFile> shape_file;
    UncertaintyCalculatorMap calculator_map;
};

} // namespace limits
} // namespace analysis
