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
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau}[GeV]", "dN/dm_{#tau#tau}[1/GeV]", false, 1.5)
    TH1D_ENTRY_CUSTOM_EX(m_vis, M_tt_Bins(), "M_{vis}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_CUSTOM_EX(m_bb_slice, M_tt_bbSlice_Bins(), "2DM_{sv}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1)

    TH1D_ENTRY_EX(pt_1, 20, 0, 200, "P_{T}(leading#tau_{h})[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 2)
    TH1D_ENTRY_EX(pt_2, 20, 0, 200, "P_{T}(subleading#tau_{h})[GeV]", "Events", false, 1)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 2)
    TH1D_ENTRY_EX(pt_b1, 20, 0, 200, "P_{T}[GeV](leading_jet)", "Events", false, 1.1)
    TH1D_ENTRY_EX(eta_b1, 25, -2.5, 2.5, "#eta(leading_jet)", "Events", false, 2)
    TH1D_ENTRY_EX(csv_b1, 25, 0, 1, "CSV(leading_jet)", "Events", false, 1.1)
    TH1D_ENTRY_EX(pt_b2, 20, 0, 200, "P_{T}[GeV](subleading_jet)", "Events", false, 1.1)
    TH1D_ENTRY_EX(eta_b2, 25, -2.5, 2.5, "#eta(subleading_jet)", "Events", false, 2)
    TH1D_ENTRY_EX(csv_b2, 25, 0, 1, "CSV(subleading_jet)", "Events", false, 1.1)
    TH1D_ENTRY_EX(pt_H_tt, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(pt_H_bb, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(pt_H_hh, 20, 0, 300, "P_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(m_bb, 30, 0, 600, "M_{jj}[GeV]", "dN/dm_{bb}[1/GeV]", false, 1.1)
    TH1D_ENTRY_EX(DeltaPhi_tt, 22, 0., 3.3, "#Delta#Phi_{#tau#tau}[rad]", "Events", false, 1.3)
    TH1D_ENTRY_EX(DeltaPhi_bb, 22, 0., 3.3, "#Delta#Phi_{bb}[rad]", "Events", false, 1.8)
    TH1D_ENTRY_EX(DeltaPhi_bb_MET, 22, 0., 3.3, "#Delta#Phi_{bb,MET}[rad]", "Events", false, 1.5)
    TH1D_ENTRY_EX(DeltaPhi_tt_MET, 22, 0., 3.3, "#Delta#Phi_{#tau#tau,MET}[rad]", "Events", false, 1.5)
    TH1D_ENTRY_EX(DeltaPhi_hh, 22, 0., 3.3, "#Delta#Phi_{#tau#taubb}[rad]", "Events", false, 1.5)
    TH1D_ENTRY_EX(DeltaR_tt, 40, 0, 6, "#DeltaR_{#tau#tau}", "Events", false, 1.1)
    TH1D_ENTRY_EX(DeltaR_bb, 40, 0, 6, "#DeltaR_{bb}[rad]", "Events", false, 1.7)
    TH1D_ENTRY_EX(DeltaR_hh, 40, 0, 6, "#DeltaR_{#tau#taubb}[rad]", "Events", false, 1.5)
    TH1D_ENTRY_EX(mt_2, 20, 0, 200, "M_{T}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_EX(pt_H_tt_MET, 20, 0, 300, "P_{T}[GeV]", "Evnets", false, 1.1)
    TH1D_ENTRY_EX(convergence, 10, -3.5, 6.5, "Fit_convergence", "Events", false, 1.1)
    TH1D_ENTRY_EX(chi2, 20, 0, 100, "#chi^{2}", "Events", false, 1.1)
    TH1D_ENTRY_EX(fit_probability, 20, 0, 1, "Fit_probability", "Events", false, 1.1)
    TH1D_ENTRY_EX(pull_balance, 20, -10, 10, "pull_balance", "Events", false, 2)
    TH1D_ENTRY_EX(pull_balance_1, 100, -10, 10, "pull_balance_1", "Events", false, 1.1)
    TH1D_ENTRY_EX(pull_balance_2, 100, -10, 10, "pull_balance_1", "Events", false, 1.1)
    TH1D_ENTRY_EX(MET, 20, 0, 100, "E_{T}^{miss}[GeV]", "Events", false, 1.1)

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

    static std::string HistogramSuffix(EventSubCategory subCategory, EventEnergyScale eventEnergyScale)
    {
        std::ostringstream ss;
        if(subCategory != EventSubCategory::NoCuts)
            ss << subCategory;
        if(eventEnergyScale != EventEnergyScale::Central) {
            if(ss.str().size())
                ss << "_";
            ss << eventEnergyScale;
        }
        return ss.str();
    }

    static std::string FullHistogramName(const std::string& hist_name, EventSubCategory subCategory,
                                         EventEnergyScale eventEnergyScale)
    {
        std::ostringstream ss;
        ss << hist_name;
        const std::string suffix = HistogramSuffix(subCategory, eventEnergyScale);
        if(suffix.size())
            ss << "_" << suffix;
        return ss.str();
    }

    virtual void Fill(const FlatEventInfo& eventInfo, double weight, EventEnergyScale eventEnergyScale,
                      EventSubCategory subCategory)
    {
        const ntuple::Flat& event = *eventInfo.event;
        const std::string key = HistogramSuffix(subCategory, eventEnergyScale);

        const double mass_tautau = event.m_sv_MC;
        m_sv(key).Fill(mass_tautau, weight);
        if (eventInfo.fitResults.has_valid_mass)
            m_ttbb_kinfit(key).Fill(eventInfo.fitResults.mass, weight);
        FillSlice(m_bb_slice(key), mass_tautau, eventInfo.Hbb.M(), weight);

        if (eventEnergyScale != EventEnergyScale::Central) return;

        pt_1(key).Fill(event.pt_1, weight);
        eta_1(key).Fill(event.eta_1, weight);
        pt_2(key).Fill(event.pt_2, weight);
        eta_2(key).Fill(event.eta_2, weight);
        DeltaPhi_tt(key).Fill(std::abs(eventInfo.lepton_momentums.at(0).DeltaPhi(eventInfo.lepton_momentums.at(1))), weight);
        DeltaR_tt(key).Fill(eventInfo.lepton_momentums.at(0).DeltaR(eventInfo.lepton_momentums.at(1)), weight);
        pt_H_tt(key).Fill(eventInfo.Htt.Pt(),weight);

        m_vis(key).Fill(eventInfo.Htt.M(),weight);
        pt_H_tt_MET(key).Fill(eventInfo.Htt_MET.Pt(), weight);
        DeltaPhi_tt_MET(key).Fill(std::abs(eventInfo.Htt.DeltaPhi(eventInfo.MET)), weight);
        mt_2(key).Fill(event.mt_2, weight);
        MET(key).Fill(eventInfo.MET.Pt(),weight);

        if(!eventInfo.has_bjet_pair) return;
        pt_b1(key).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(), weight);
        eta_b1(key).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Eta(), weight);
        csv_b1(key).Fill(eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.first), weight);
        pt_b2(key).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Pt(), weight);
        eta_b2(key).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Eta(), weight);
        csv_b2(key).Fill(eventInfo.event->csv_Bjets.at(eventInfo.selected_bjets.second), weight);
        DeltaPhi_bb(key).Fill(std::abs(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaPhi(
                                       eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second))), weight);
        DeltaR_bb(key).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaR(
                                     eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second)), weight);
        pt_H_bb(key).Fill(eventInfo.Hbb.Pt(),weight);
        m_bb(key).Fill(eventInfo.Hbb.M(), weight);
        DeltaPhi_bb_MET(key).Fill(std::abs(eventInfo.Hbb.DeltaPhi(eventInfo.MET)), weight);
        DeltaPhi_hh(key).Fill(std::abs(eventInfo.Htt.DeltaPhi(eventInfo.Hbb)), weight);
        DeltaR_hh(key).Fill(eventInfo.Htt.DeltaR(eventInfo.Hbb), weight);
        m_ttbb(key).Fill(eventInfo.resonance.M(), weight);
        pt_H_hh(key).Fill(eventInfo.resonance.Pt(), weight);

        convergence(key).Fill(eventInfo.fitResults.convergence,weight);
        chi2(key).Fill(eventInfo.fitResults.chi2,weight);
        fit_probability(key).Fill(eventInfo.fitResults.fit_probability,weight);
        pull_balance(key).Fill(eventInfo.fitResults.pull_balance,weight);
        pull_balance_1(key).Fill(eventInfo.fitResults.pull_balance_1,weight);
        pull_balance_2(key).Fill(eventInfo.fitResults.pull_balance_2,weight);
    }

    void FillSubCategories(const FlatEventInfo& eventInfo, double weight, EventEnergyScale eventEnergyScale)
    {
        using namespace cuts::massWindow;

        const double mass_tautau = eventInfo.event->m_sv_MC;
        Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::NoCuts);
        const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                && eventInfo.Hbb.M() > m_bb_low && eventInfo.Hbb.M() < m_bb_high;
        if(inside_mass_window)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::MassWindow);
        else
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::OutsideMassWindow);
        if(eventInfo.fitResults.has_valid_mass)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::KinematicFitConverged);
        if(eventInfo.fitResults.has_valid_mass && inside_mass_window)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::KinematicFitConvergedWithMassWindow);

        if(eventEnergyScale == EventEnergyScale::Central) {
            FillSubCategories(eventInfo, weight, EventEnergyScale::BtagEfficiencyUp);
            FillSubCategories(eventInfo, weight, EventEnergyScale::BtagEfficiencyDown);
        }
    }

    void FillEnergyScales(const FlatEventInfo& eventInfo, double weight, const std::set<EventEnergyScale>& energyScales)
    {
        for (EventEnergyScale energyScale : energyScales)
            FillSubCategories(eventInfo, weight, energyScale);
    }

    void FillCentralAndJetEnergyScales(const FlatEventInfo& eventInfo, double weight)
    {
        static const std::set<EventEnergyScale> energyScales =
            { EventEnergyScale::Central, EventEnergyScale::JetUp, EventEnergyScale::JetDown };
        FillEnergyScales(eventInfo, weight, energyScales);
    }

    void FillAllEnergyScales(const FlatEventInfo& eventInfo, double weight)
    {
        FillEnergyScales(eventInfo, weight, AllEventEnergyScales);
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

};

} // namespace analysis
