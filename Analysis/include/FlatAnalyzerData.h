/*!
 * \file FlatAnalyzerData.h
 * \brief Definition of histogram containers for flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-12-02 created
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

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau}[GeV]", "dN/dm_{#tau#tau}[1/GeV]", false, 1.5, true)
    TH1D_ENTRY_CUSTOM_EX(m_vis, M_tt_Bins(), "M_{vis}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_bb_slice, M_tt_bbSlice_Bins(), "2DM_{sv}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1, true)

    TH1D_ENTRY_EX(pt_1, 20, 0, 200, "P_{T}(leading#tau_{h})[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 2, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 200, "P_{T}(subleading#tau_{h})[GeV]", "Events", false, 1, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 2, SaveAll)
    TH1D_ENTRY_EX(pt_b1, 20, 0, 200, "P_{T}[GeV](leading_jet)", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(eta_b1, 25, -2.5, 2.5, "#eta(leading_jet)", "Events", false, 2, SaveAll)
    TH1D_ENTRY_EX(csv_b1, 25, 0, 1, "CSV(leading_jet)", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pt_b2, 20, 0, 200, "P_{T}[GeV](subleading_jet)", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(eta_b2, 25, -2.5, 2.5, "#eta(subleading_jet)", "Events", false, 2, SaveAll)
    TH1D_ENTRY_EX(csv_b2, 25, 0, 1, "CSV(subleading_jet)", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pt_H_bb, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pt_H_hh, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(m_bb, 30, 0, 600, "M_{jj}[GeV]", "dN/dm_{bb}[1/GeV]", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_tt, 22, 0., 3.3, "#Delta#Phi_{#tau#tau}[rad]", "Events", false, 1.3, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_bb, 22, 0., 3.3, "#Delta#Phi_{bb}[rad]", "Events", false, 1.8, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_bb_MET, 22, 0., 3.3, "#Delta#Phi_{bb,MET}[rad]", "Events", false, 1.5, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_tt_MET, 22, 0., 3.3, "#Delta#Phi_{#tau#tau,MET}[rad]", "Events", false, 1.5, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_hh, 22, 0., 3.3, "#Delta#Phi_{#tau#taubb}[rad]", "Events", false, 1.5, SaveAll)
    TH1D_ENTRY_EX(DeltaR_tt, 40, 0, 6, "#DeltaR_{#tau#tau}", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(DeltaR_bb, 40, 0, 6, "#DeltaR_{bb}[rad]", "Events", false, 1.7, SaveAll)
    TH1D_ENTRY_EX(DeltaR_hh, 40, 0, 6, "#DeltaR_{#tau#taubb}[rad]", "Events", false, 1.5, SaveAll)
    TH1D_ENTRY_EX(mt_2, 20, 0, 200, "M_{T}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt_MET, 20, 0, 300, "P_{T}[GeV]", "Evnets", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(convergence, 10, -3.5, 6.5, "Fit_convergence", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(chi2, 20, 0, 100, "#chi^{2}", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(fit_probability, 20, 0, 1, "Fit_probability", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pull_balance, 20, -10, 10, "pull_balance", "Events", false, 2, SaveAll)
    TH1D_ENTRY_EX(pull_balance_1, 100, -10, 10, "pull_balance_1", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(pull_balance_2, 100, -10, 10, "pull_balance_1", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(MET, 20, 0, 100, "E_{T}^{miss}[GeV]", "Events", false, 1.1, SaveAll)
    TH2D_ENTRY_EX(csv_b1_vs_ptb1, 20, 0, 200, 25, 0, 1, "P_{T}[GeV](leading_jet)", "CSV(leading_jet)", false, 1, SaveAll)

    static constexpr bool SaveAll = true;

    explicit FlatAnalyzerData(bool _fill_all) : fill_all(_fill_all) {}

    FlatAnalyzerData(TFile& outputFile, const std::string& directoryName, bool _fill_all)
        : AnalyzerData(outputFile, directoryName), fill_all(_fill_all) {}

    typedef root_ext::SmartHistogram<TH1D>& (FlatAnalyzerData::*HistogramAccessor)();

    virtual const std::vector<double>& M_tt_Bins() const
    {
        static const std::vector<double> bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
                                                  160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };
        return bins;
    }

    virtual const std::vector<double>& M_tt_bbSlice_Bins() const
    {
        static const size_t number_of_slices = 5;
        const static std::vector<double> bins = CreateMttSlicedBins(number_of_slices, M_tt_Bins());
        return bins;
    }

    virtual const std::vector<double>& M_ttbb_Bins() const
    {
        static const std::vector<double> bins = { 200, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 500, 550,
                                                  600, 650, 700 };
        return bins;
    }

    virtual void Fill(const FlatEventInfo& eventInfo, double weight)
    {
        const ntuple::Flat& event = *eventInfo.event;

        const double mass_tautau = event.m_sv_MC;
        m_sv().Fill(mass_tautau, weight);
        if (eventInfo.fitResults.has_valid_mass)
            m_ttbb_kinfit().Fill(eventInfo.fitResults.mass, weight);
//        FillSlice(m_bb_slice(), mass_tautau, eventInfo.Hbb.M(), weight);

        if(!fill_all) return;

        pt_1().Fill(event.pt_1, weight);
        eta_1().Fill(event.eta_1, weight);
        pt_2().Fill(event.pt_2, weight);
        eta_2().Fill(event.eta_2, weight);
        DeltaPhi_tt().Fill(std::abs(eventInfo.lepton_momentums.at(0).DeltaPhi(eventInfo.lepton_momentums.at(1))), weight);
        DeltaR_tt().Fill(eventInfo.lepton_momentums.at(0).DeltaR(eventInfo.lepton_momentums.at(1)), weight);
        pt_H_tt().Fill(eventInfo.Htt.Pt(),weight);

        m_vis().Fill(eventInfo.Htt.M(),weight);
        pt_H_tt_MET().Fill(eventInfo.Htt_MET.Pt(), weight);
        DeltaPhi_tt_MET().Fill(std::abs(eventInfo.Htt.DeltaPhi(eventInfo.MET)), weight);
        mt_2().Fill(event.mt_2, weight);
        MET().Fill(eventInfo.MET.Pt(),weight);

        if(!eventInfo.has_bjet_pair) return;
        pt_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(), weight);
        eta_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Eta(), weight);
        csv_b1().Fill(eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.first), weight);
        pt_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Pt(), weight);
        eta_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Eta(), weight);
        csv_b2().Fill(eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.second), weight);
        DeltaPhi_bb().Fill(std::abs(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaPhi(
                                       eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second))), weight);
        DeltaR_bb().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaR(
                                     eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second)), weight);
        pt_H_bb().Fill(eventInfo.Hbb.Pt(),weight);
        m_bb().Fill(eventInfo.Hbb.M(), weight);
        DeltaPhi_bb_MET().Fill(std::abs(eventInfo.Hbb.DeltaPhi(eventInfo.MET)), weight);
        DeltaPhi_hh().Fill(std::abs(eventInfo.Htt.DeltaPhi(eventInfo.Hbb)), weight);
        DeltaR_hh().Fill(eventInfo.Htt.DeltaR(eventInfo.Hbb), weight);
        m_ttbb().Fill(eventInfo.resonance.M(), weight);
        pt_H_hh().Fill(eventInfo.resonance.Pt(), weight);

        convergence().Fill(eventInfo.fitResults.convergence,weight);
        chi2().Fill(eventInfo.fitResults.chi2,weight);
        fit_probability().Fill(eventInfo.fitResults.fit_probability,weight);
        pull_balance().Fill(eventInfo.fitResults.pull_balance,weight);
        pull_balance_1().Fill(eventInfo.fitResults.pull_balance_1,weight);
        pull_balance_2().Fill(eventInfo.fitResults.pull_balance_2,weight);

        csv_b1_vs_ptb1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(),
                              eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.first), weight);
    }

protected:
    static void FillSlice(TH1D& hist, double m_sv, double m_Hbb, double weight)
    {
        static const std::vector<double> slice_regions = { 60, 100, 140, 200, 600 };
        static const double slice_size = 350;

        if(m_sv < 0 || m_sv >= slice_size) return;
        const auto slice_region = std::find_if(slice_regions.begin(), slice_regions.end(),
                                               [&](double x) { return m_Hbb < x; });
        if(slice_region == slice_regions.end()) return;
        const ptrdiff_t slice_id = slice_region - slice_regions.begin();
        const double slice_shift = slice_size * slice_id;
        hist.Fill(m_sv + slice_shift, weight);
    }

    static std::vector<double> CreateMttSlicedBins(size_t number_of_slices, const std::vector<double>& mtt_bins)
    {
        if(!mtt_bins.size())
            throw exception("Invalid mtt bins.");
        if(!number_of_slices)
            throw exception("Invalid number of slices.");

        std::vector<double> sliced_bins;
        for(size_t n = 0; n < number_of_slices; ++n) {
            size_t k = n == 0 ? 0 : 1;
            const double shift = n * mtt_bins.back();
            for(; k < mtt_bins.size(); ++k) {
                const double value = mtt_bins.at(k) + shift;
                sliced_bins.push_back(value);
            }
        }
        return sliced_bins;
    }

    bool fill_all;
};

class FlatAnalyzerData_semileptonic : public FlatAnalyzerData {
public:
    TH1D_ENTRY_EX(mt_1, 50, 0, 50, "M_{T}[GeV]", "Events", false, 1.1, SaveAll)

    explicit FlatAnalyzerData_semileptonic(bool _fill_all) : FlatAnalyzerData(_fill_all) {}

    FlatAnalyzerData_semileptonic(TFile& outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData(outputFile, directoryName, fill_all) {}

    virtual void Fill(const FlatEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight);
        if(!fill_all) return;
        const ntuple::Flat& event = *eventInfo.event;
        mt_1().Fill(event.mt_1, weight);
    }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 200, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 500, 550,
                                                  600, 650, 700 };
        return bins;
    }
};

class FlatAnalyzerData_semileptonic_2tag : public FlatAnalyzerData_semileptonic {
public:
    explicit FlatAnalyzerData_semileptonic_2tag(bool _fill_all) : FlatAnalyzerData_semileptonic(_fill_all) {}

    FlatAnalyzerData_semileptonic_2tag(TFile& outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData_semileptonic(outputFile, directoryName, fill_all) {}

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }
};

class FlatAnalyzerData_tautau : public FlatAnalyzerData {
public:
    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}[GeV]", "Events", false, 1.1, SaveAll)
    TH1D_ENTRY_EX(iso_tau1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1, SaveAll)
    TH1D_ENTRY_EX(iso_tau2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1, SaveAll)

    explicit FlatAnalyzerData_tautau(bool _fill_all) : FlatAnalyzerData(_fill_all) {}

    FlatAnalyzerData_tautau(TFile& outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData(outputFile, directoryName, fill_all) {}

    virtual void Fill(const FlatEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight);
        if(!fill_all) return;
        const ntuple::Flat& event = *eventInfo.event;
        mt_1().Fill(event.mt_1, weight);
        iso_tau1().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1,weight);
        iso_tau2().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2,weight);
    }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 200, 250, 280, 310, 340, 370, 400, 500, 600, 700 };
        return bins;
    }
};

class FlatAnalyzerData_tautau_2tag : public FlatAnalyzerData_tautau {
public:
    explicit FlatAnalyzerData_tautau_2tag(bool _fill_all) : FlatAnalyzerData_tautau(_fill_all) {}

    FlatAnalyzerData_tautau_2tag(TFile& outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData_tautau(outputFile, directoryName, fill_all) {}

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }
};

} // namespace analysis
