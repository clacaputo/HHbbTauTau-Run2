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

#include "AnalysisBase/include/SmartHistogram.h"

namespace analysis {

class StackedPlotDescriptor {
public:
    typedef root_ext::SmartHistogram<TH1D> Histogram;
    typedef std::shared_ptr<Histogram> hist_ptr;
    typedef std::vector<hist_ptr> hist_ptr_vector;

    StackedPlotDescriptor(const std::string& page_title, bool draw_title, const std::string& channelNameLatex,
                          bool _draw_ratio)
        : text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC")), channelText(new TPaveText(0.2, 0.85, 0.25, 0.89, "NDC")),
          channelName(channelNameLatex), draw_ratio(_draw_ratio)
    {
        page.side.fit_range_x = false;
        page.title = page_title;
        page.layout.has_title = draw_title;
        if (draw_ratio){
            if (page.layout.has_title) {
                page.side.layout.main_pad = root_ext::Box(0.02,0.19,0.95,0.94);
                page.side.layout.ratio_pad = root_ext::Box(0.02,0.02,1,0.28);
            } else {
                page.side.layout.main_pad = root_ext::Box(0.02,0.21,1,1);
                page.side.layout.ratio_pad = root_ext::Box(0.02,0.02,1,0.30);
            }
        }
        else {
            if (page.layout.has_title) {
                page.side.layout.main_pad = root_ext::Box(0.02,0.02,0.95,0.94);
            } else {
                page.side.layout.main_pad = root_ext::Box(0.02,0.02,1,1);
            }
        }

        if (draw_ratio){
            legend = std::shared_ptr<TLegend>(new TLegend (0.6, 0.55, 0.85, 0.90));
        }
        else {
            legend = std::shared_ptr<TLegend>(new TLegend (0.55, 0.60, 0.80, 0.95));
        }
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);

        text->SetTextSize(0.03);
        text->SetTextFont(62);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        text->AddText(0.01,0.5, "CMS, 19.7 fb^{-1} at 8 TeV");

        channelText->SetTextSize(0.05);
        channelText->SetTextFont(62);
        channelText->SetFillColor(0);
        channelText->SetBorderSize(0);
        channelText->SetMargin(0.01);
        channelText->SetTextAlign(12); // align left
        channelText->AddText(channelName.c_str());
    }

    const std::string& GetTitle() const { return page.title; }

    void AddBackgroundHistogram(const Histogram& original_histogram, const std::string& legend_title, Color_t color)
    {
        hist_ptr histogram = PrepareHistogram(original_histogram);
        histogram->SetFillColor(color);
        background_histograms.push_back(histogram);
        //stack->Add(histogram.get());
        legend->AddEntry(histogram.get(), legend_title.c_str(), "f");
        if(!sum_backgound_histogram)
            sum_backgound_histogram = hist_ptr(new Histogram(*histogram));
        else
            sum_backgound_histogram->Add(histogram.get());
    }

    void AddSignalHistogram(const Histogram& original_signal, const std::string& legend_title, Color_t color,
                            unsigned scale_factor)
    {
        hist_ptr histogram = PrepareHistogram(original_signal);
        //histogram->SetLineColor(color);
        histogram->SetLineColor(4);
        histogram->SetLineStyle(2);
        histogram->SetLineWidth(3);
        histogram->Scale(scale_factor);
        signal_histograms.push_back(histogram);
        std::ostringstream ss;
        if(scale_factor != 1)
            ss << scale_factor << "x ";
        ss << legend_title;
        legend->AddEntry(histogram.get(), ss.str().c_str(), "F");
    }

    void AddDataHistogram(const Histogram& original_data, const std::string& legend_title,
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


        if (background_histograms.size()) {
            const auto& bkg_hist = background_histograms.front();
            stack = std::shared_ptr<THStack>(new THStack(bkg_hist->GetName(), bkg_hist->GetTitle()));

            for (auto iter = background_histograms.rbegin(); iter != background_histograms.rend(); ++iter){
                stack->Add(iter->get());
            }
            stack->Draw("HIST");

            Double_t maxY = stack->GetMaximum();
            if (data_histogram){
                const Int_t maxBin = data_histogram->GetMaximumBin();
                const Double_t maxData = data_histogram->GetBinContent(maxBin) + data_histogram->GetBinError(maxBin);
                maxY = std::max(maxY, maxData);
            }

            for (const hist_ptr& signal : signal_histograms)
                maxY = std::max(maxY,signal->GetMaximum());

            stack->SetMaximum(maxY * bkg_hist->MaxYDrawScaleFactor());

            const Double_t minY = page.side.use_log_scaleY ? 1 : 0;
            stack->SetMinimum(minY);


            if (draw_ratio){
                stack->GetXaxis()->SetTitle("");
                stack->GetXaxis()->SetLabelColor(kWhite);

            }
            else {
                stack->GetXaxis()->SetTitle(page.side.axis_titleX.c_str());
                stack->GetXaxis()->SetTitleOffset(1.05);
                stack->GetXaxis()->SetTitleSize(0.04);
                stack->GetXaxis()->SetLabelSize(0.04);
            }

            stack->GetYaxis()->SetTitleSize(0.05);
            stack->GetYaxis()->SetTitleOffset(1.6);
            stack->GetYaxis()->SetLabelSize(0.04);
            stack->GetYaxis()->SetTitle(page.side.axis_titleY.c_str());

            sum_backgound_histogram->SetMarkerSize(0);
            sum_backgound_histogram->SetFillColor(13);
            sum_backgound_histogram->SetFillStyle(3013);
            sum_backgound_histogram->SetLineWidth(1);
            sum_backgound_histogram->Draw("e2same");
        }

        for(const hist_ptr& signal : signal_histograms)
            signal->Draw("SAME HIST");

        if(data_histogram) {
            data_histogram->SetMarkerStyle(20);
            data_histogram->SetMarkerSize(1);
            data_histogram->Draw("samepPE0");
        }

        text->Draw("same");
        channelText->Draw("same");

        legend->Draw("same");

        const std::string axis_titleX = page.side.axis_titleX;
        if (data_histogram && draw_ratio){
            ratio_pad = std::shared_ptr<TPad>(root_ext::Adapter::NewPad(page.side.layout.ratio_pad));
            ratio_pad->Draw();

            ratio_pad->cd();


            ratio_histogram = hist_ptr(new Histogram(*data_histogram));
            ratio_histogram->Divide(sum_backgound_histogram.get());

            ratio_histogram->GetYaxis()->SetRangeUser(0.75,1.3);
            ratio_histogram->GetYaxis()->SetNdivisions(505);
            ratio_histogram->GetYaxis()->SetLabelSize(0.11);
            ratio_histogram->GetYaxis()->SetTitleSize(0.14);
            ratio_histogram->GetYaxis()->SetTitleOffset(0.55);
            ratio_histogram->GetYaxis()->SetTitle("Obs/Bkg");
            ratio_histogram->GetXaxis()->SetNdivisions(510);
            ratio_histogram->GetXaxis()->SetTitle(axis_titleX.c_str());
            ratio_histogram->GetXaxis()->SetTitleSize(0.1);
            ratio_histogram->GetXaxis()->SetTitleOffset(0.98);
            //ratio_histogram->GetXaxis()->SetLabelColor(kBlack);
            ratio_histogram->GetXaxis()->SetLabelSize(0.1);
            ratio_histogram->SetMarkerStyle(20);
            ratio_histogram->SetMarkerColor(1);
            ratio_histogram->SetMarkerSize(1);

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
        if (data_histogram && draw_ratio)
            ratio_pad->Draw();
    }

private:
    hist_ptr PrepareHistogram(const Histogram& original_histogram)
    {
        hist_ptr histogram(new Histogram(original_histogram));
        histogram->SetLineColor(root_ext::colorMapName.at("black"));
        histogram->SetLineWidth(1.);
        ReweightWithBinWidth(histogram);
        UpdateDrawInfo(histogram);
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
            histogram->SetBinContent(n, need_blind ? -1 : histogram->GetBinContent(n));
            histogram->SetBinError(n, need_blind ? -1 : histogram->GetBinError(n));
        }
    }

    void UpdateDrawInfo(hist_ptr hist)
    {
        page.side.use_log_scaleY = hist->UseLogY();
        page.side.axis_titleX = hist->GetXTitle();
        page.side.axis_titleY = hist->GetYTitle();
    }

private:
    root_ext::SingleSidedPage page;
    hist_ptr data_histogram;
    hist_ptr sum_backgound_histogram;
    hist_ptr ratio_histogram;
    hist_ptr_vector background_histograms;
    hist_ptr_vector signal_histograms;
    std::shared_ptr<THStack> stack;
    std::shared_ptr<TLegend> legend;
    std::shared_ptr<TPaveText> text;
    std::shared_ptr<TPaveText> channelText;
    std::shared_ptr<TPad> main_pad, ratio_pad;
    std::string channelName;
    bool draw_ratio;
};

} // namespace analysis
