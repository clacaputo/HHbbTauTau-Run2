/*!
 * \file Cuts.h
 * \brief Common tools and definitions to apply cuts.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-14 created
 */

#pragma once

#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>

#include <TH1D.h>

namespace cuts {

class cut_failed : public std::exception {
public:
    cut_failed(int parameter_id) noexcept
        : _param_id(parameter_id)
    {
        std::ostringstream ss;
        ss << "Cut requirements are not fulfilled for parameter id = " << _param_id << ".";
        message = ss.str();
    }

    ~cut_failed() noexcept {}

    virtual const char* what() const noexcept { return message.c_str(); }
    int param_id() const noexcept { return _param_id; }

private:
    int _param_id;
    std::string message;
};

void apply_cut(bool expected_condition, TH1D& histogram, int param_id)
{
    if(!expected_condition)
        throw cut_failed(param_id);
    histogram.Fill(param_id);
}

void apply_cut(bool expected_condition, TH1D& counter_histogram, int param_id, TH1D& selection_histogram,
               const std::string& param_label)
{
    apply_cut(expected_condition, counter_histogram, param_id);
    selection_histogram.GetXaxis()->SetBinLabel(param_id + 1, param_label.c_str());
}

void fill_selection_histogram(TH1D& selection_histogram, const TH1D& counter_histogram)
{
    if(selection_histogram.GetNbinsX() != counter_histogram.GetNbinsX())
        throw std::runtime_error("Selection and counter histograms have different number of bins.");
    for(Int_t n = 1; n <= counter_histogram.GetNbinsX(); ++n) {
        if(counter_histogram.GetBinContent(n) > 0.5)
            selection_histogram.AddBinContent(n);
    }
}

} // namespace cuts
