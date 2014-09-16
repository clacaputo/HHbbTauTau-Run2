/*!
 * \file StackedPlotDescriptor.h
 * \brief Code to produce stacked plots using CMS preliminary style.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-16 created
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

#include <iostream>
#include <iomanip>
#include <memory>

#include <TH1.h>
#include <THStack.h>

#include "RootPrintSource.h"
#include "TdrStyle.h"

namespace analysis {

struct HistogramDescriptor {
    std::string name;
    std::string title;
    std::string Xaxis_title;
    std::string Yaxis_title;
    bool useLogY;
};

typedef std::vector<HistogramDescriptor> HistogramDescriptorVector;

inline std::ostream& operator<<(std::ostream& s, const HistogramDescriptor& hist) {
    s << "Name: " << hist.name << ", Title: " << hist.title << ", useLog: " << hist.useLogY  ;
    return s;
}

inline HistogramDescriptorVector ReadHistogramConfigurations(const std::string& config_name)
{
    HistogramDescriptorVector histograms;
    std::ifstream cfg(config_name);
    while (cfg.good()) {
        std::string cfgLine;
        std::getline(cfg, cfgLine);
        if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
        std::istringstream ss(cfgLine);
        HistogramDescriptor hist;
        ss >> hist.name;
        ss >> hist.title;
        ss >> hist.Xaxis_title;
        ss >> hist.Yaxis_title;
        ss >> std::boolalpha >> hist.useLogY;
        histograms.push_back(hist);
    }
    return histograms;
}

class StackedPlotDescriptor {
public:
    typedef std::shared_ptr<TH1D> hist_ptr;
    typedef std::vector<hist_ptr> hist_ptr_vector;

    StackedPlotDescriptor(const analysis::HistogramDescriptor& _hist_descriptor, const std::string& page_title)
        : hist_descriptor(_hist_descriptor),
          data_histogram(nullptr),
          stack(new THStack(hist_descriptor.name.c_str(), hist_descriptor.title.c_str())),
          legend(new TLegend (0.60, 0.65, 0.67, 0.90)),
          text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"))
    {
        page.side.use_log_scaleY = hist_descriptor.useLogY;
        page.side.fit_range_x = false;
        page.title = page_title;
        page.side.axis_titleX = hist_descriptor.Xaxis_title;
        page.side.axis_titleY = hist_descriptor.Yaxis_title;

        page.side.layout.main_pad.right_top.x = 0.95;
        page.side.layout.main_pad.right_top.y = 1;
        page.side.layout.main_pad.left_bottom.x = 0.05;
        page.side.layout.main_pad.left_bottom.y = 0.1;
        page.side.layout.ratio_pad.right_top.x = 0.95;
        page.side.layout.ratio_pad.right_top.y = 0.1;
        page.side.layout.ratio_pad.left_bottom.x = 0.05;
        page.side.layout.ratio_pad.left_bottom.y = 0;

        legend->SetFillColor(0);
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);

        text->SetTextSize(0.03);
        text->SetTextFont(42);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        text->AddText(0.01,0.5, "CMS Preliminary");
        text->AddText(0.25, 0.6, "#sqrt{s} = 8 TeV  L = 19.7 fb^{-1}");
    }

    const std::string& GetTitle() const { return page.title; }

    void AddBackgroundHistogram(TH1D* original_histogram, const std::string& legend_title, Color_t color)
    {
        hist_ptr histogram = PrepareHistogram(original_histogram);
        histogram->SetFillColor(color);
        background_histograms.push_back(histogram);
        stack->Add(histogram.get());
        legend->AddEntry(histogram.get(), legend_title.c_str(), "f");
        if(!sum_backgound_histogram)
            sum_backgound_histogram = hist_ptr( static_cast<TH1D*>(histogram->Clone()) );
        else
            sum_backgound_histogram->Add(histogram.get());
    }

    void AddSignalHistogram(TH1D* original_signal, const std::string& legend_title, Color_t color, double scale_factor)
    {
        hist_ptr histogram = PrepareHistogram(original_signal);
        histogram->SetLineColor(color);
        histogram->Scale(scale_factor);
        signal_histograms.push_back(histogram);
        legend->AddEntry(histogram.get(), legend_title.c_str(), "F");
    }

    void AddDataHistogram(TH1D* original_data, const std::string& legend_title,
                          bool blind, const std::pair<double, double>& blind_region)
    {
        if(data_histogram)
            throw std::runtime_error("Only one data histogram per stack is supported.");

        data_histogram = PrepareHistogram(original_data);
        legend->AddEntry(data_histogram.get(), legend_title.c_str(), "LP");

        if(blind)
            BlindHistogram(data_histogram, blind_region);
    }

    void Draw(TCanvas& canvas)
    {
        cms_tdr::setTDRStyle();

        main_pad = std::shared_ptr<TPad>(root_ext::Adapter::NewPad(page.side.layout.main_pad));
        if(page.side.use_log_scaleX)
            main_pad->SetLogx();
        if(page.side.use_log_scaleY)
            main_pad->SetLogy();
        main_pad->Draw();
        main_pad->cd();

        stack->Draw("HIST");
        const Double_t maxY = std::max(stack->GetMaximum(), data_histogram->GetMaximum());
        stack->SetMaximum(maxY*1.1);
        const Double_t minY = page.side.use_log_scaleY ? 1 : 0;
        stack->SetMinimum(minY);

        stack->GetXaxis()->SetTitle(page.side.axis_titleX.c_str());
        stack->GetYaxis()->SetTitle(page.side.axis_titleY.c_str());

        for(const hist_ptr& signal : signal_histograms)
            signal->Draw("SAME HIST");

        if(data_histogram) {
            data_histogram->SetMarkerStyle(7);
            data_histogram->Draw("samepPE0");
        }

        text->Draw("same");

        legend->Draw("same");

        ratio_pad = std::shared_ptr<TPad>(root_ext::Adapter::NewPad(page.side.layout.ratio_pad));
        ratio_pad->Draw();

        ratio_pad->cd();

        ratio_histogram = hist_ptr(static_cast<TH1D*>(data_histogram->Clone()));
        ratio_histogram->Divide(sum_backgound_histogram.get());

        ratio_histogram->GetYaxis()->SetRangeUser(0.8,1.2);
        ratio_histogram->GetYaxis()->SetNdivisions(3);
        ratio_histogram->GetYaxis()->SetLabelSize(0.2);
        ratio_histogram->GetYaxis()->SetTitleSize(0.25);
        ratio_histogram->GetYaxis()->SetTitleOffset(0.15);
        ratio_histogram->GetYaxis()->SetTitle("Obs/Est");
        ratio_histogram->GetXaxis()->SetNdivisions(-1);
        ratio_histogram->GetXaxis()->SetTitle("");
        ratio_histogram->GetXaxis()->SetLabelSize(0.0001);
        ratio_histogram->SetMarkerStyle(7);
        ratio_histogram->SetMarkerColor(2);

        ratio_histogram->Draw("histp");

        TLine* line = new TLine();
        line->DrawLine(ratio_histogram->GetXaxis()->GetXmin(), 1, ratio_histogram->GetXaxis()->GetXmax(), 1);

        canvas.cd();
        main_pad->Draw();
        canvas.cd();
        ratio_pad->Draw();
    }

private:
    static hist_ptr PrepareHistogram(TH1D* original_histogram)
    {
        hist_ptr histogram( static_cast<TH1D*>(original_histogram->Clone()) );
        histogram->SetLineColor(root_ext::colorMapName.at("black"));
        histogram->SetLineWidth(1.);
        ReweightWithBinWidth(histogram);
        return histogram;
    }

    static void ReweightWithBinWidth(hist_ptr histogram)
    {
        for(Int_t n = 1; n <= histogram->GetNbinsX(); ++n) {
            const double new_value = histogram->GetBinContent(n) / histogram->GetBinWidth(n);
            const double new_bin_error = histogram->GetBinError(n) / histogram->GetBinWidth(n);
            histogram->SetBinContent(n, new_value);
            histogram->SetBinError(n, new_bin_error);
        }
    }

    static void BlindHistogram(hist_ptr histogram, const std::pair<double, double>& blind_region)
    {
        for(Int_t n = 1; n <= histogram->GetNbinsX(); ++n) {
            const double x = histogram->GetBinCenter(n);
            const bool need_blind = x > blind_region.first && x < blind_region.second;
            histogram->SetBinContent(n, need_blind ? 0. : histogram->GetBinContent(n));
            histogram->SetBinError(n, need_blind ? 0. : histogram->GetBinError(n));
        }
    }

private:
    analysis::HistogramDescriptor hist_descriptor;
    root_ext::SingleSidedPage page;
    hist_ptr data_histogram;
    hist_ptr sum_backgound_histogram;
    hist_ptr ratio_histogram;
    hist_ptr_vector background_histograms;
    hist_ptr_vector signal_histograms;
    std::shared_ptr<THStack> stack;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<TPaveText> text;

    std::shared_ptr<TPad> main_pad, ratio_pad;
};

} // namespace analysis
