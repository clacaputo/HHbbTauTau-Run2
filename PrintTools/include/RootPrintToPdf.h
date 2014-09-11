/*!
 * \file RootPrintToPdf.h
 * \brief Print ROOT histograms to PDF.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 *
 * Copyright 2013, 2014 Konstantin Androsov <konstantin.androsov@gmail.com>
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

#include <sstream>
#include <stdexcept>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <Rtypes.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <THStack.h>

#include "RootPrintSource.h"
#include "TdrStyle.h"

namespace root_ext {

class PdfPrinter {
public:
    PdfPrinter(const std::string& _output_file_name)
        : output_file_name(_output_file_name)
    {
        canvas = new TCanvas();
        canvas->Print((output_file_name + "[").c_str());
    }

    template<typename Source>
    void Print(const Page& page, const Source& source)
    {
        gROOT->SetStyle(page.layout.global_style.c_str());
        gStyle->SetOptStat(page.layout.stat_options);
        gStyle->SetOptFit(page.layout.fit_options);
        canvas->cd();

        canvas->SetTitle(page.title.c_str());
        if(page.layout.has_title) {
            TPaveLabel *title = Adapter::NewPaveLabel(page.layout.title_box, page.title);
            title->SetTextFont(page.layout.title_font);
            title->Draw();
        }

        Page::RegionCollection page_regions = page.Regions();
        for(Page::RegionCollection::const_iterator iter = page_regions.begin(); iter != page_regions.end(); ++iter)
        {
            canvas->cd();
            DrawHistograms(*(*iter), source);
        }

        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << page.title;
        canvas->Print(output_file_name.c_str(), print_options.str().c_str());
    }


    void PrintStack(const Page& page, THStack& stack, TH1D& data_hist, TH1D&ratioData_Bkg, TLegend& leg, TPaveText& pave_text)
    {
        TH1::SetDefaultSumw2();
        cms_tdr::setTDRStyle();

        canvas->cd();

        canvas->SetTitle(page.title.c_str());
        if(page.layout.has_title) {
            TPaveLabel *title = Adapter::NewPaveLabel(page.layout.title_box, page.title);
            title->SetTextFont(page.layout.title_font);
            title->Draw();
        }

        Page::RegionCollection page_regions = page.Regions();
        for(Page::RegionCollection::const_iterator iter = page_regions.begin(); iter != page_regions.end(); ++iter)
        {
            canvas->cd();
            DrawStack(*(*iter), stack, data_hist, ratioData_Bkg, leg, pave_text);
        }

        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << page.title;
        canvas->Print(output_file_name.c_str(), print_options.str().c_str());
    }

    ~PdfPrinter()
    {
        canvas->Print((output_file_name+"]").c_str());
    }

private:
    template<typename Source>
    void DrawHistograms(const PageSide& page_side, const Source& source)
    {
        typedef root_ext::HistogramPlotter<typename Source::Histogram, typename Source::ValueType> Plotter;
        typedef root_ext::HistogramFitter<typename Source::Histogram, typename Source::ValueType> Fitter;

        TPad* stat_pad = 0;
        if(page_side.layout.has_stat_pad) {
            stat_pad = Adapter::NewPad(page_side.layout.stat_pad);
            stat_pad->Draw();
        }

        TPad *pad = Adapter::NewPad(page_side.layout.main_pad);
        if(page_side.use_log_scaleX)
            pad->SetLogx();
        if(page_side.use_log_scaleY)
            pad->SetLogy();
        pad->Draw();
        pad->cd();

        Plotter plotter(page_side.histogram_title, page_side.axis_titleX, page_side.axis_titleY);
        for(unsigned n = 0; n < source.Size(); ++n)
        {
            const typename Plotter::Entry entry = source.Get(n, page_side.histogram_name);
            plotter.Add(entry);
        }

        Fitter::SetRanges(plotter.Histograms(), page_side.fit_range_x, page_side.fit_range_y, page_side.xRange,
                          page_side.yRange, page_side.use_log_scaleY);
        plotter.Superpose(pad, stat_pad, page_side.layout.has_legend, page_side.layout.legend_pad,
                          page_side.draw_options);
    }

    void DrawStack(const PageSide& page_side, THStack& stack, TH1D& data_hist, TH1D&ratioData_Bkg, TLegend& leg, TPaveText& pave_text)
    {
        TPad* stat_pad = 0;
        if(page_side.layout.has_stat_pad) {
            stat_pad = Adapter::NewPad(page_side.layout.stat_pad);
            stat_pad->Draw();
        }

        TPad *pad = Adapter::NewPad(page_side.layout.main_pad);
        if(page_side.use_log_scaleX)
            pad->SetLogx();
        if(page_side.use_log_scaleY)
            pad->SetLogy();
        pad->Draw();
        pad->cd();

        stack.Draw("HIST");

//        const Int_t firstBin = stack.GetXaxis()->FindBin(page_side.xRange.min);
//        const Int_t lastBin = stack.GetXaxis()->FindBin(page_side.xRange.max);
//        stack.GetXaxis()->SetRange(firstBin,lastBin);
        const Double_t maxY = std::max(stack.GetMaximum(),data_hist.GetMaximum());
        stack.SetMaximum(maxY*1.1);
        const Double_t minY = page_side.use_log_scaleY ? 1 : 0;
        stack.SetMinimum(minY);

        stack.GetXaxis()->SetTitle(page_side.axis_titleX.c_str());
        stack.GetYaxis()->SetTitle(page_side.axis_titleY.c_str());

        data_hist.SetMarkerStyle(7);
        data_hist.Draw("samepPE0");
        pave_text.Draw("same");

//        TPad* padRatio = new TPad("padRatio","",0,0,1,0.1);
//        padRatio->Draw();
//        padRatio->cd();

//        // Draw the ratio of the historgrams

//        ratioData_Bkg.GetYaxis()->SetRangeUser(0.9,1.1);
//        ratioData_Bkg.GetYaxis()->SetNdivisions(3);
//        ratioData_Bkg.GetYaxis()->SetLabelSize(0.1);
//        ratioData_Bkg.GetYaxis()->SetTitleSize(0.1);
//        ratioData_Bkg.GetYaxis()->SetTitleOffset(0.5);
//        ratioData_Bkg.GetYaxis()->SetTitle("Ratio");
//        ratioData_Bkg.GetXaxis()->SetNdivisions(-1);
//        ratioData_Bkg.GetXaxis()->SetTitle("");
//        ratioData_Bkg.GetXaxis()->SetLabelSize(0.0001);
//        ratioData_Bkg.SetMarkerStyle(7);
//        ratioData_Bkg.SetMarkerColor(2);
//        ratioData_Bkg.Draw("histp");
//        TLine* line = new TLine();
//        line->DrawLine(ratioData_Bkg.GetXaxis()->GetXmin(),1,ratioData_Bkg.GetXaxis()->GetXmax(),1);


        if(stat_pad)
            stat_pad->cd();
        leg.Draw("same");
    }

private:
    TCanvas* canvas;
    std::string output_file_name;
};

} // namespace root_ext
