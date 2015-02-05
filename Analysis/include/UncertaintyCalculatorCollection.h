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

private:
    std::shared_ptr<TFile> shape_file;
    UncertaintyCalculatorMap calculator_map;
};

} // namespace limits
} // namespace analysis
