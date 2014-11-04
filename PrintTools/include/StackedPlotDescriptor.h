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
#include <TLine.h>

#include "RootPrintSource.h"
#include "TdrStyle.h"

namespace analysis {

struct HistogramDescriptor;
typedef std::vector<HistogramDescriptor> HistogramDescriptorVector;

struct HistogramDescriptor {
    std::string name;
    std::string title;
    std::string Xaxis_title;
    std::string Yaxis_title;
    bool useLogY;
    double max_Y;

    static HistogramDescriptorVector ReadFromFile(const std::string& config_name)
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
            ss >> hist.max_Y;
            histograms.push_back(hist);
        }
        return histograms;
    }
};

inline std::ostream& operator<<(std::ostream& s, const HistogramDescriptor& hist) {
    s << "Name: " << hist.name << ", Title: " << hist.title << ", useLog: " << hist.useLogY  ;
    return s;
}

class StackedPlotDescriptor {
public:
    typedef std::shared_ptr<TH1D> hist_ptr;
    typedef std::vector<hist_ptr> hist_ptr_vector;

    StackedPlotDescriptor(const analysis::HistogramDescriptor& _hist_descriptor, const std::string& page_title,
                          bool draw_title)
        : hist_descriptor(_hist_descriptor),
          data_histogram(nullptr),
          stack(new THStack(hist_descriptor.name.c_str(), hist_descriptor.title.c_str())),
          legend(new TLegend (0.6, 0.55, 0.8, 0.90)),
          text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"))
    {
        page.side.use_log_scaleY = hist_descriptor.useLogY;
        page.side.fit_range_x = false;
        page.title = page_title;
        page.side.axis_titleX = hist_descriptor.Xaxis_title;
        page.side.axis_titleY = hist_descriptor.Yaxis_title;
        page.layout.has_title = draw_title;
        if (page.layout.has_title) {
            page.side.layout.main_pad.right_top.x = 0.95;
            page.side.layout.main_pad.right_top.y = 0.95;
            page.side.layout.main_pad.left_bottom.x = 0.05;
            page.side.layout.main_pad.left_bottom.y = 0.25;
            page.side.layout.ratio_pad.right_top.x = 0.95;
            page.side.layout.ratio_pad.right_top.y = 0.3;
            page.side.layout.ratio_pad.left_bottom.x = 0.05;
            page.side.layout.ratio_pad.left_bottom.y = 0.05;
        } else {
            page.side.layout.main_pad.right_top.x = 1;
            page.side.layout.main_pad.right_top.y = 1;
            page.side.layout.main_pad.left_bottom.x = 0.02;
            page.side.layout.main_pad.left_bottom.y = 0.21;
            page.side.layout.ratio_pad.right_top.x = 1;
            page.side.layout.ratio_pad.right_top.y = 0.30;
            page.side.layout.ratio_pad.left_bottom.x = 0.02;
            page.side.layout.ratio_pad.left_bottom.y = 0.02;
        }

        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
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

    void AddSignalHistogram(TH1D* original_signal, const std::string& legend_title, Color_t color, unsigned scale_factor)
    {
        hist_ptr histogram = PrepareHistogram(original_signal);
        histogram->SetLineColor(color);
        histogram->SetLineStyle(kDashed);
        histogram->Scale(scale_factor);
        signal_histograms.push_back(histogram);
        std::ostringstream ss;
        if(scale_factor != 1)
            ss << scale_factor << "x ";
        ss << legend_title;
        legend->AddEntry(histogram.get(), ss.str().c_str(), "F");
    }

    void AddDataHistogram(TH1D* original_data, const std::string& legend_title,
                          bool blind, const std::vector< std::pair<double, double> >& blind_regions)
    {
        if(data_histogram)
            throw std::runtime_error("Only one data histogram per stack is supported.");

        data_histogram = PrepareHistogram(original_data);
        legend->AddEntry(data_histogram.get(), legend_title.c_str(), "LP");

        if(blind)
            BlindHistogram(data_histogram, blind_regions);
    }

    bool NeedDraw() const
    {
        return background_histograms.size() || data_histogram || signal_histograms.size();
    }

    void Draw(TCanvas& canvas)
    {
        cms_tdr::setTDRStyle();

        if(page.layout.has_title) {
            TPaveLabel *title = root_ext::Adapter::NewPaveLabel(page.layout.title_box, page.title);
            title->SetTextFont(page.layout.title_font);
            title->Draw();
        }
        main_pad = std::shared_ptr<TPad>(root_ext::Adapter::NewPad(page.side.layout.main_pad));
        if(page.side.use_log_scaleX)
            main_pad->SetLogx();
        if(page.side.use_log_scaleY)
            main_pad->SetLogy();
        main_pad->Draw();
        main_pad->cd();


        if (background_histograms.size()){
            stack->Draw("HIST");
            if (data_histogram){
                const Double_t maxY = std::max(stack->GetMaximum(), data_histogram->GetMaximum());
                stack->SetMaximum(maxY*hist_descriptor.max_Y);
            }

            const Double_t minY = page.side.use_log_scaleY ? 1 : 0;
            stack->SetMinimum(minY);

            //stack->GetXaxis()->SetTitle(page.side.axis_titleX.c_str());
            stack->GetXaxis()->SetTitleOffset(1.03);
            stack->GetXaxis()->SetTitleSize(0.03);
            stack->GetXaxis()->SetTitle("");
            stack->GetXaxis()->SetLabelSize(0.03);
            stack->GetXaxis()->SetLabelColor(kWhite);

            stack->GetYaxis()->SetTitleSize(0.03);
            stack->GetYaxis()->SetTitleOffset(1.5);
            stack->GetYaxis()->SetLabelSize(0.04);
            stack->GetYaxis()->SetTitle(page.side.axis_titleY.c_str());

        }

        for(const hist_ptr& signal : signal_histograms)
            signal->Draw("SAME HIST");

        if(data_histogram) {
            data_histogram->SetMarkerStyle(7);
            data_histogram->Draw("samepPE0");
        }

        text->Draw("same");

        legend->Draw("same");

        const std::string axis_titleX = page.side.axis_titleX;
        if (data_histogram){
            ratio_pad = std::shared_ptr<TPad>(root_ext::Adapter::NewPad(page.side.layout.ratio_pad));
            ratio_pad->Draw();

            ratio_pad->cd();


            ratio_histogram = hist_ptr(static_cast<TH1D*>(data_histogram->Clone()));
            ratio_histogram->Divide(sum_backgound_histogram.get());

            ratio_histogram->GetYaxis()->SetRangeUser(0.75,1.3);
            ratio_histogram->GetYaxis()->SetNdivisions(505);
            ratio_histogram->GetYaxis()->SetLabelSize(0.09);
            ratio_histogram->GetYaxis()->SetTitleSize(0.12);
            ratio_histogram->GetYaxis()->SetTitleOffset(0.3);
            ratio_histogram->GetYaxis()->SetTitle("Obs/Bkg");
            ratio_histogram->GetXaxis()->SetNdivisions(510);
            ratio_histogram->GetXaxis()->SetTitle(axis_titleX.c_str());
            ratio_histogram->GetXaxis()->SetTitleSize(0.09);
            ratio_histogram->GetXaxis()->SetTitleOffset(0.95);
            //ratio_histogram->GetXaxis()->SetLabelColor(kBlack);
            ratio_histogram->GetXaxis()->SetLabelSize(0.09);
            ratio_histogram->SetMarkerStyle(7);
            ratio_histogram->SetMarkerColor(1);

            ratio_histogram->Draw("E0P");

            TLine* line = new TLine();
            line->SetLineStyle(3);
            line->DrawLine(ratio_histogram->GetXaxis()->GetXmin(), 1, ratio_histogram->GetXaxis()->GetXmax(), 1);
            TLine* line1 = new TLine();
            line1->SetLineStyle(3);
            line1->DrawLine(ratio_histogram->GetXaxis()->GetXmin(), 1.2, ratio_histogram->GetXaxis()->GetXmax(), 1.2);
            TLine* line2 = new TLine();
            line2->SetLineStyle(3);
            line2->DrawLine(ratio_histogram->GetXaxis()->GetXmin(), 0.8, ratio_histogram->GetXaxis()->GetXmax(), 0.8);
            ratio_pad->SetTopMargin(0.04);
            ratio_pad->SetBottomMargin(0.3);
            ratio_pad->Update();
        }

        canvas.cd();
        main_pad->Draw();

        canvas.cd();
        if (data_histogram)
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

    static void BlindHistogram(hist_ptr histogram, const std::vector< std::pair<double, double> >& blind_regions)
    {
        for(Int_t n = 1; n <= histogram->GetNbinsX(); ++n) {
            const double x = histogram->GetBinCenter(n);
            const auto blind_predicate = [&](const std::pair<double, double>& region) -> bool {
                return x > region.first && x < region.second;
            };

            const bool need_blind = std::any_of(blind_regions.begin(), blind_regions.end(), blind_predicate);
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
