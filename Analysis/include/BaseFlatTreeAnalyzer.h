/*!
 * \file BaseFlatTreeAnalyzer.h
 * \brief Definition of BaseFlatTreeAnalyzer class, the base class for flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-03 created
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

#ifndef __APPLE__
#define override
#endif

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/AnalysisTypes.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

static const std::vector<double> mass_bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
                                               160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };
static const std::vector<double> mass_bins_2j2t = { 0,20,40,60,80,100,120,140,160,180,200,250,300,350};
//static const std::vector<double> mass_ttbb_bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280,
//                                                    300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
//                                                    650, 700, 750, 800, 850, 900, 950, 1000 };

static const std::vector<double> mass_ttbb_bins = { 200,250,270,290,310,330,350,370,390,410,430,450,500,550,600,650,700 };

static const std::vector<double> mass_bins_slice_2fette = { 
      0,	10,      20,	 30,	 40,	 50,	 60,	 70,	 80,	 90,	100,	110,	120,	130,
	140,	150,	160,	170,	180,	190,	200,	225,	250,	275,	300,	325,	350,
	360,	370,	380,	390,	400,	410,	420,	430,	440,	450,	460,	470,	480,
	490,	500,	510,	520,	530,	540,	550,	575,	600,	625,	650,	675,	700
};

static const std::vector<double> mass_bins_slice_5fette_lb = { 
       0,	  20,	  40,	  60,	  80,	 100,	 120,	 140,	 160,	 180,	 200,	 250,	 300,	350,
     370,	 390,	 410,	 430,	 450,	 470,	 490,	 510,	 530,	 550,	 600,	 650,	 700,
     720,	 740,	 760,	 780,	 800,	 820,	 840,	 860,	 880,	 900,	 950,	1000,	1050,
	1070,	1090,	1110,	1130,	1150,	1170,	1190,	1210,	1230,	1250,	1300,	1350,	1400,
    1420,	1440,	1460,	1480,	1500,	1520,	1540,	1560,	1580,	1600,	1650,	1700,	1750
};

static const std::vector<double> mass_bins_slice_5fette_fb = {
       0,     10,     20,     30,     40,     50,     60,     70,     80,     90,    100,    110,    120,    130,
     140,    150,    160,    170,    180,    190,    200,    225,    250,    275,    300,    325,    350,
     360,	 370,	 380,	 390,	 400,	 410,	 420,	 430,	 440,	 450,	 460,	 470,	 480,	 490,
     500,	 510,	 520,	 530,	 540,	 550,	 575,	 600,	 625,	 650,	 675,	 700,
     710,	 720,	 730,	 740,	 750,	 760,	 770,	 780,	 790,	 800,	 810,	 820,	 830,	 840,
     850,	 860,	 870,	 880,	 890,	 900,	 925,	 950,	 975,	1000,	1025,	1050,
    1060,	1070,	1080,	1090,	1100,	1110,	1120,	1130,	1140,	1150,	1160,	1170,	1180,	1190,
    1200,	1210,	1220,	1230,	1240,	1250,	1275,	1300,	1325,	1350,	1375,	1400,
    1410,	1420,	1430,	1440,	1450,	1460,	1470,	1480,	1490,	1500,	1510,	1520,	1530,	1540,
    1550,	1560,	1570,	1580,	1590,	1600,	1625,	1650,	1675,	1700,	1725,	1750
 };

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY(pt_1, 20, 0, 200)
    TH1D_ENTRY(eta_1, 25, -2.5, 2.5)
    TH1D_ENTRY(pt_2, 20, 0, 200)
    TH1D_ENTRY(eta_2, 25, -2.5, 2.5)
    TH1D_ENTRY(pt_b1, 20, 0, 200)
    TH1D_ENTRY(eta_b1, 25, -2.5, 2.5)
    TH1D_ENTRY(csv_b1, 25, 0, 1)
    TH1D_ENTRY(pt_b2, 20, 0, 200)
    TH1D_ENTRY(eta_b2, 25, -2.5, 2.5)
    TH1D_ENTRY(csv_b2, 25, 0, 1)
    TH1D_ENTRY(pt_H_tt, 20, 0, 300)
    TH1D_ENTRY(pt_H_bb, 20, 0, 300)
    TH1D_ENTRY(pt_H_hh, 20, 0, 300)
    TH1D_ENTRY_CUSTOM(m_sv, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_up, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_down, mass_bins)
    TH1D_ENTRY_CUSTOM(m_vis, mass_bins)
    TH1D_ENTRY(m_bb, 30, 0, 600)
    TH1D_ENTRY_CUSTOM(m_ttbb, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_up, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_down, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only_up, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only_down, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_massCut, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_up_massCut, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_down_massCut, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only_massCut, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only_up_massCut, mass_ttbb_bins)
    TH1D_ENTRY_CUSTOM(m_ttbb_kinfit_only_down_massCut, mass_ttbb_bins)
    TH1D_ENTRY(DeltaPhi_tt, 22, 0., 3.3)
    TH1D_ENTRY(DeltaPhi_bb, 22, 0., 3.3)
    TH1D_ENTRY(DeltaPhi_bb_MET, 22, 0., 3.3)
    TH1D_ENTRY(DeltaPhi_tt_MET, 22, 0., 3.3)
    TH1D_ENTRY(DeltaPhi_hh, 22, 0., 3.3)
    TH1D_ENTRY(DeltaR_tt, 40, 0, 6)
    TH1D_ENTRY(DeltaR_bb, 40, 0, 6)
    TH1D_ENTRY(DeltaR_hh, 40, 0, 6)
    TH1D_ENTRY_CUSTOM(m_bb_slice, mass_bins_slice_5fette_fb)
    TH1D_ENTRY_CUSTOM(m_bb_slice_up, mass_bins_slice_5fette_fb)
    TH1D_ENTRY_CUSTOM(m_bb_slice_down, mass_bins_slice_5fette_fb)
    TH1D_ENTRY(MVA_BDT, 40, -1, 1)
    TH1D_ENTRY(mt_2, 20, 0, 200)
    TH1D_ENTRY(pt_H_tt_MET, 20, 0, 300)
    TH1D_ENTRY(convergence, 10, -3.5, 6.5)
    TH1D_ENTRY(chi2, 20, 0, 100)
    TH1D_ENTRY(fit_probability, 20, 0, 1)
    TH1D_ENTRY(pull_balance, 20, -10, 10)
    TH1D_ENTRY(pull_balance_1, 100, -10, 10)
    TH1D_ENTRY(pull_balance_2, 100, -10, 10)
    TH1D_ENTRY(MET, 20, 0, 100)

    static std::string FullHistogramName(const std::string& hist_name, EventEnergyScale eventEnergyScale)
    {
        std::ostringstream ss;
        ss << hist_name << "_" << eventEnergyScale;
        return ss.str();
    }

    void CreateHistogramsWithCustomBinning(const std::string& hist_name, const std::vector<double>& binning)
    {
        for(EventEnergyScale eventEnergyScale : AllEventEnergyScales)
            Get((TH1D*)nullptr, hist_name, eventEnergyScale, binning);
    }

    virtual void Fill(const FlatEventInfo& eventInfo, double weight, bool fill_all, EventEnergyScale eventEnergyScale)
    {
        const ntuple::Flat& event = *eventInfo.event;
        double mass_tautau = event.m_sv_MC;
        m_sv(eventEnergyScale).Fill(mass_tautau, weight);
        if(!fill_all) return;

        const double m_ttbb_kinFit = eventInfo.fitResults.mass;
        if (eventInfo.fitResults.has_valid_mass)
            m_ttbb_kinfit_only(eventEnergyScale).Fill(m_ttbb_kinFit, weight);
        if (mass_tautau > cuts::massWindow::m_tautau_low && mass_tautau < cuts::massWindow::m_tautau_high &&
                eventInfo.Hbb.M() > cuts::massWindow::m_bb_low && eventInfo.Hbb.M() < cuts::massWindow::m_bb_high) {
            if (eventInfo.fitResults.has_valid_mass)
                m_ttbb_kinfit_only_massCut(eventEnergyScale).Fill(m_ttbb_kinFit, weight);
        }
        FillSlice(m_bb_slice(eventEnergyScale), mass_tautau, eventInfo.Hbb.M(), weight);

        if (eventEnergyScale != analysis::EventEnergyScale::Central) return;

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
//        MVA_BDT().Fill(eventInfo.mva_BDT, weight);
    }

    void FillAll(const FlatEventInfo& eventInfo, double weight, bool fill_all)
    {
        for (EventEnergyScale energyScale : AllEventEnergyScales){
            Fill(eventInfo, weight, fill_all, energyScale);
        }
    }

private:
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
};

class BaseFlatTreeAnalyzer {
public:
    typedef std::map<EventRegion, std::shared_ptr<FlatAnalyzerData> > RegionAnaData;

    typedef std::map<std::string, RegionAnaData> AnaDataForEventCategory;
    typedef std::map<EventCategory, AnaDataForEventCategory> FullAnaData;
    typedef std::map<EventRegion, PhysicalValue> PhysicalValueMap;

    static const std::string ReferenceHistogramName() { return "m_sv_Central"; }

    static const EventCategorySet& EventCategoriesToProcess() { return AllEventCategories; }


    virtual const EventCategorySet & CategoriesToRelax()
    {
        static const EventCategorySet categories= {EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag};
        //static const EventCategorySet categories;
        return categories;
    }

    virtual const EventRegionSet& EssentialEventRegions()
    {
        static const EventRegionSet regions = { EventRegion::OS_Isolated, EventRegion::SS_Isolated };
        return regions;
    }

    BaseFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& _outputFileName, Channel channel_id,
                         const std::string& signal_list)
        : inputPath(_inputPath), outputFileName(_outputFileName),
          dataCategoryCollection(source_cfg, signal_list, channel_id)
    {
        TH1::SetDefaultSumw2();
        histograms = HistogramDescriptor::ReadFromFile(hist_cfg);
    }

    void Run()
    {
        std::cout << "Processing data categories... " << std::endl;
        for(const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            if(!dataCategory->sources_sf.size()) continue;
            std::cout << *dataCategory << std::endl;
            for (const EventRegion& eventRegion : AllEventRegions){
                auto& anaData = GetAnaData(EventCategory::TwoJets_TwoBtag, dataCategory->name, eventRegion);
                anaData.CreateHistogramsWithCustomBinning("m_sv", mass_bins_2j2t);
                auto& anaData_loose = GetAnaData(EventCategory::TwoJets_TwoLooseBtag, dataCategory->name, eventRegion);
                anaData_loose.CreateHistogramsWithCustomBinning("m_sv", mass_bins_2j2t);
            }
            for(const auto& source_entry : dataCategory->sources_sf) {
                const std::string fullFileName = inputPath + "/" + source_entry.first;
                std::shared_ptr<TFile> file(new TFile(fullFileName.c_str(), "READ"));
                if(file->IsZombie())
                    throw exception("Input file '") << source_entry.first << "' not found.";
                std::shared_ptr<ntuple::FlatTree> tree(new ntuple::FlatTree(*file, "flatTree"));
                ProcessDataSource(*dataCategory, tree, source_entry.second);
            }
        }

        std::cout << "Calculating embedded scale factor... " << std::endl;
        const auto embeddedSF = CalculateEmbeddedScaleFactor(ReferenceHistogramName());
        //const auto TTembeddedSF = CalculateTTEmbeddedScaleFactor(ReferenceHistogramName());
//        const double embeddedSF = 1;
        std::cout << "Embedded SF = " << embeddedSF << std::endl;

        std::cout<<"Histrogram for ZTT - - - - - - - -"<<std::endl;
        for (EventCategory eventCategory : EventCategoriesToProcess()){
            for (const auto& hist : histograms) {
                //RescaleHistogramForTTembedded(eventCategory, hist.name, TTembeddedSF);
                CreateHistogramForZTT(eventCategory, hist.name, embeddedSF, false);
            }
        }

        std::cout << "Estimating QCD, Wjets and composit data categories... " << std::endl;
        for (EventCategory eventCategory : EventCategoriesToProcess()){
            std::cout<<"ScaleFactors for W - - - - - - - -"<<std::endl;
            const auto wjets_scale_factors = CalculateWjetsScaleFactors(eventCategory, ReferenceHistogramName());
            for (EventRegion eventRegion : QcdRegions){
                std::cout << eventCategory << ", " << eventRegion << ", scale factor = " <<
                             wjets_scale_factors.at(eventRegion) << std::endl;
            }
            std::cout<<"Estimate for W - - - - - - - -"<<std::endl;
            for (const auto& hist : histograms) {
                EstimateWjets(eventCategory, hist.name, wjets_scale_factors);
            }
        }

        for (EventCategory eventCategory : EventCategoriesToProcess()) {

            std::cout<<"ScaleFactor for QCD - - - - - - - -"<<std::endl;
            const auto qcd_scale_factor = CalculateQCDScaleFactor(ReferenceHistogramName(), eventCategory);
            //std::cout << eventCategory << " SS_Iso / SS_NotIso = " << qcd_scale_factor << std::endl;
            std::cout << eventCategory << " OS_NotIso / SS_NotIso = " << qcd_scale_factor << std::endl;

            for (const auto& hist : histograms) {
                EstimateQCD(hist.name, eventCategory, qcd_scale_factor);
                ProcessCompositDataCategories(eventCategory, hist.name);
            }
            std::cout << std::endl;
        }

        std::cout << "Saving tables... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");

        std::cout << "Saving datacards... " << std::endl;
        static const root_ext::SmartHistogram<TH1D> emptyDatacard_mSV("emptyDatacard_mSV",mass_bins);
        static const root_ext::SmartHistogram<TH1D> emptyDatacard_mttbb("emptyDatacard_mttbb", mass_ttbb_bins);
        static const root_ext::SmartHistogram<TH1D> emptyDatacard_slice("emptyDatacard_slice", mass_bins_slice_5fette_fb);

        ProduceFileForLimitsCalculation("m_sv", emptyDatacard_mSV);
        //ProduceFileForLimitsCalculation("m_ttbb_kinfit", emptyDatacard_mttbb);
        ProduceFileForLimitsCalculation("m_ttbb_kinfit_only", emptyDatacard_mttbb);
        //ProduceFileForLimitsCalculation("m_ttbb_kinfit_massCut", emptyDatacard_mttbb);
        ProduceFileForLimitsCalculation("m_ttbb_kinfit_only_massCut", emptyDatacard_mttbb);
        ProduceFileForLimitsCalculation("m_bb_slice", emptyDatacard_slice);

        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots(false, false);
        PrintStackedPlots(true, false);
        PrintStackedPlots(false, true);
        PrintStackedPlots(true, true);
    }

protected:
    virtual Channel ChannelId() const = 0;
    virtual EventRegion DetermineEventRegion(const ntuple::Flat& event, EventCategory eventCategory) = 0;
    virtual bool PassMvaCut(const FlatEventInfo& eventInfo, EventCategory eventCategory) { return true; }

    virtual PhysicalValue CalculateQCDScaleFactor(const std::string& hist_name, EventCategory eventCategory) = 0;
    virtual void EstimateQCD(const std::string& hist_name, EventCategory eventCategory,
                             const PhysicalValue& scale_factor) = 0;
    virtual PhysicalValueMap CalculateWjetsScaleFactors(EventCategory eventCategory, const std::string& hist_name) = 0;


    virtual PhysicalValue CalculateQCDScaleFactorEx(const std::string& hist_name,
                                                    EventCategory num_eventCategory, EventCategory den_eventCategory,
                                                    EventRegion num_eventRegion, EventRegion den_eventRegion)
    {
        const DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);
        auto hist_num_data = GetHistogram(num_eventCategory, data.name, num_eventRegion, hist_name);
        auto hist_den_data = GetHistogram(den_eventCategory, data.name, den_eventRegion, hist_name);

        if(!hist_num_data)
            throw analysis::exception("Unable to find num histograms for QCD scale factor estimation for ")
                << hist_name << ", num_eventCategory: " << num_eventCategory << ", num_eventRegion: " <<
                   num_eventRegion;

        if(!hist_den_data)
            throw analysis::exception("Unable to find den histograms for QCD scale factor estimation for ")
                << hist_name << ", den_eventCategory: " << den_eventCategory << ", den_eventRegion: " <<
                   den_eventRegion;

        TH1D* hist_num = GetHistogram(num_eventCategory, qcd.name, num_eventRegion, hist_name);
        if(!hist_num) {
            hist_num = &CloneHistogram(num_eventCategory, qcd.name, num_eventRegion, *hist_num_data);
            SubtractBackgroundHistograms(*hist_num, num_eventCategory, num_eventRegion, qcd.name, true);
        }

        TH1D* hist_den = GetHistogram(den_eventCategory, qcd.name, den_eventRegion, hist_name);
        if(!hist_den) {
            hist_den = &CloneHistogram(den_eventCategory, qcd.name, den_eventRegion, *hist_den_data);
            SubtractBackgroundHistograms(*hist_den, den_eventCategory, den_eventRegion, qcd.name, true);
        }

        const PhysicalValue n_num = Integral(*hist_num, false);
        const PhysicalValue n_den = Integral(*hist_den, false);
        if(n_num.value < 0 || n_den.value < 0)
            throw exception("Negative number of estimated events in QCD SF estimation for ") << num_eventCategory;
        return n_num / n_den;
    }

    virtual void EstimateQCDEx(const std::string& hist_name, EventCategory eventCategory, EventCategory shapeCategory,
                               EventRegion shapeRegion, const PhysicalValue& scale_factor)
    {
        const DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_shape_data = GetHistogram(shapeCategory, data.name, shapeRegion, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_Isolated, *hist_shape_data);
        SubtractBackgroundHistograms(histogram, shapeCategory, shapeRegion, qcd.name);
        histogram.Scale(scale_factor.value);
    }

    virtual void EstimateWjets(EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValueMap& scale_factor_map)
    {
        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);

        const EventCategory shapeEventCategory = CategoriesToRelax().count(eventCategory)
                ? analysis::MediumToLoose_EventCategoryMap.at(eventCategory) : eventCategory;

        for(EventRegion eventRegion : QcdRegions) {
            if(!scale_factor_map.count(eventRegion))
                throw exception("W-jet SF not found for ") << eventRegion;

            const PhysicalValue& scale_factor = scale_factor_map.at(eventRegion);
            TH1D* hist = nullptr;
            for (const analysis::DataCategory* wjets_category : wjets_mc_categories){
                if(auto hist_mc = GetHistogram(shapeEventCategory, wjets_category->name, eventRegion, hist_name)) {
                    if (!hist)
                        hist = &CloneHistogram(eventCategory, wjets.name, eventRegion, *hist_mc);
                    else
                        hist->Add(hist_mc);
                }
            }
            if (hist)
                hist->Scale(scale_factor.value);
        }
    }

    const std::string& ChannelName() const { return detail::ChannelNameMap.at(ChannelId()); }
    const std::string& ChannelNameLatex() const { return detail::ChannelNameMapLatex.at(ChannelId()); }

    virtual std::shared_ptr<FlatAnalyzerData> MakeAnaData()
    {
        return std::shared_ptr<FlatAnalyzerData>(new FlatAnalyzerData());
    }

    FlatAnalyzerData& GetAnaData(EventCategory eventCategory, const std::string& dataCategoryName,
                                 EventRegion eventRegion)
    {
        std::shared_ptr<FlatAnalyzerData>& anaData = fullAnaData[eventCategory][dataCategoryName][eventRegion];
        if(!anaData)
            anaData = MakeAnaData();
        return *anaData;
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                 EventRegion eventRegion, const std::string& histogramName)
    {
        return GetAnaData(eventCategory, dataCategoryName, eventRegion).GetPtr<TH1D>(histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                       const std::string& histogramName)
    {
        return GetHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, histogramName);
    }

    TH1D& CloneHistogram(EventCategory eventCategory, const std::string& dataCategoryName, EventRegion eventRegion,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return GetAnaData(eventCategory, dataCategoryName, eventRegion).Clone(originalHistogram);
    }

    TH1D& CloneSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, originalHistogram);
    }


    void ProcessDataSource(const DataCategory& dataCategory, std::shared_ptr<ntuple::FlatTree> tree, double scale_factor)
    {
        static const bool applyMVAcut = false;

        const analysis::DataCategory& DYJets_incl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_incl);
        const analysis::DataCategory& DYJets_excl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_excl);
        const analysis::DataCategory& WJets_incl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC_incl);

        for(Long64_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);
            const ntuple::Flat& event = tree->data;

            const EventCategoryVector eventCategories = DetermineEventCategories(event.csv_Bjets,
                                                                                 event.nBjets_retagged,
                                                                                 cuts::Htautau_Summer13::btag::CSVL,
                                                                                 cuts::Htautau_Summer13::btag::CSVM);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1), false);

            const double weight = dataCategory.IsData() ? 1 : event.weight * scale_factor;
            if(std::isnan(event.weight)) {
                std::cerr << "ERROR: event " << event.evt << " will not be processed because event weight is NaN."
                          << std::endl;
                continue;
            }

            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                if(eventRegion == EventRegion::Unknown) continue;
                const bool fill_all = EssentialEventRegions().count(eventRegion);


                UpdateMvaInfo(eventInfo, eventCategory, false, false, false);
                if(applyMVAcut && !PassMvaCut(eventInfo, eventCategory)) continue;

                const bool isMixedInclusiveSample = dataCategory.name == DYJets_incl.name
                                                 || dataCategory.name == WJets_incl.name;
                const bool haveExtraJets = eventInfo.event->n_extraJets_MC > 5 && eventInfo.event->n_extraJets_MC < 10;
                const double corrected_weight = isMixedInclusiveSample && haveExtraJets ? weight / 2 : weight;

                if(dataCategory.name == DYJets_excl.name || dataCategory.name == DYJets_incl.name)
                    FillDYjetHistograms(eventInfo, eventCategory, eventRegion, corrected_weight);
                if (dataCategory.IsData())
                    GetAnaData(eventCategory, dataCategory.name, eventRegion).FillAll(eventInfo,corrected_weight,fill_all);
                else
                    GetAnaData(eventCategory, dataCategory.name, eventRegion).Fill(eventInfo, corrected_weight,
                                                                               fill_all,eventInfo.eventEnergyScale);
            }

        }
    }

    void FillDYjetHistograms(const FlatEventInfo& eventInfo, EventCategory eventCategory, EventRegion eventRegion,
                             double weight)
    {
        const analysis::DataCategory& ZL = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZL);
        const analysis::DataCategory& ZJ = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZJ);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const analysis::DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);

        const std::map<ntuple::EventType, std::string> type_category_map = {
            { ntuple::EventType::ZL, ZL.name }, { ntuple::EventType::ZJ, ZJ.name },
            { ntuple::EventType::ZTT, ZTT_MC.name }, {ntuple::EventType::ZTT_L, ZTT_L.name}
        };

        if(type_category_map.count(eventInfo.eventType)) {
            const std::string& name = type_category_map.at(eventInfo.eventType);
            const bool fill_all = EssentialEventRegions().count(eventRegion);
            auto& anaData = GetAnaData(EventCategory::TwoJets_TwoBtag, name, eventRegion);
            anaData.CreateHistogramsWithCustomBinning("m_sv", mass_bins_2j2t);
            auto& anaData_loose = GetAnaData(EventCategory::TwoJets_TwoLooseBtag, name, eventRegion);
            anaData_loose.CreateHistogramsWithCustomBinning("m_sv", mass_bins_2j2t);
            GetAnaData(eventCategory, name, eventRegion).Fill(eventInfo, weight, fill_all, eventInfo.eventEnergyScale);
        }
    }

    PhysicalValue CalculateTTEmbeddedScaleFactor(const std::string& hist_name)
    {
        const analysis::DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);
        const analysis::DataCategory& TT_lept_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_lept);

        TH1D* hist_embedded = GetSignalHistogram(EventCategory::Inclusive, TTembedded.name, hist_name);
        TH1D* hist_tt = GetSignalHistogram(EventCategory::Inclusive, TT_lept_MC.name, hist_name);
        if(!hist_embedded || !hist_tt )
            throw std::runtime_error("TTembedded or tt_lept_MC hist not found");

        const PhysicalValue n_tt = Integral(*hist_tt, false);
        const PhysicalValue n_embedded = Integral(*hist_embedded, false);
        return n_tt / n_embedded;
    }

    void RescaleHistogramForTTembedded(EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValue& scale_factor)
    {
        const analysis::DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);
        const double embedded_scaleFactor = scale_factor.value;

        for(EventRegion eventRegion : AllEventRegions) {
            auto embedded_hist = GetHistogram(eventCategory, TTembedded.name, eventRegion, hist_name);
            if(!embedded_hist)
                throw std::runtime_error("TTembedded hist not found");
            TH1D* hist = static_cast<TH1D*>(embedded_hist->Clone());
            hist->Scale(embedded_scaleFactor);
        }
    }

    PhysicalValue CalculateEmbeddedScaleFactor(const std::string& hist_name)
    {
        const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        TH1D* hist_embedded = GetSignalHistogram(EventCategory::Inclusive, embedded.name, hist_name);
        TH1D* hist_ztautau = GetSignalHistogram(EventCategory::Inclusive, ZTT_MC.name, hist_name);
        if(!hist_embedded)
            throw exception("embedded hist not found");
        if(!hist_ztautau )
            throw exception("ztt hist not found");

        const PhysicalValue n_ztautau = Integral(*hist_ztautau, false);
        const PhysicalValue n_embedded = Integral(*hist_embedded, false);
        return n_ztautau / n_embedded;
    }

    void CreateHistogramForZTT(EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValue& scale_factor, bool useEmbedded)
    {
        const analysis::DataCategory& embedded = useEmbedded
                ? dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded)
                : dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const double embedded_scaleFactor = useEmbedded ? scale_factor.value : 1;
        const analysis::DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);
        const analysis::DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);
//        const analysis::DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);

        for(EventRegion eventRegion : AllEventRegions) {

            auto ztt_l_hist = GetHistogram(eventCategory, ZTT_L.name, eventRegion, hist_name);
            auto embedded_hist = GetHistogram(eventCategory, embedded.name, eventRegion, hist_name);
//            auto TTembedded_hist = GetHistogram(eventCategory, TTembedded.name, eventRegion, hist_name);

            if (embedded_hist){
                TH1D& ztt_hist = CloneHistogram(eventCategory, ZTT.name, eventRegion, *embedded_hist);
                ztt_hist.Scale(embedded_scaleFactor);
                if (ztt_l_hist){
                    ztt_hist.Add(ztt_l_hist);
                }
//                if (TTembedded_hist)
//                    ztt_hist.Add(TTembedded_hist, -1);
            }
            if (!embedded_hist && ztt_l_hist){
                CloneHistogram(eventCategory, ZTT.name, eventRegion, *ztt_l_hist);
            }
        }
    }

    void UpdateMvaInfo(FlatEventInfo& eventInfo, EventCategory eventCategory, bool calc_BDT, bool calc_BDTD,
                       bool calc_BDTMitFisher)
    {
        static double const default_value = std::numeric_limits<double>::lowest();
        std::ostringstream category_name;
        category_name << eventCategory;

        auto getMVA = [&](bool calc_MVA, MVA_Selections::MvaMethod method) -> double {
            if(calc_MVA) {
                auto mvaReader = MVA_Selections::MvaReader::Get(ChannelName(), category_name.str(), method);
                if(mvaReader)
                    return mvaReader->GetMva(eventInfo.lepton_momentums.at(0), eventInfo.lepton_momentums.at(1),
                                             eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1),
                                             eventInfo.MET);
            }
            return default_value;
        };

        eventInfo.mva_BDT = getMVA(calc_BDT, MVA_Selections::BDT);
        eventInfo.mva_BDTD = getMVA(calc_BDTD, MVA_Selections::BDTD);
        eventInfo.mva_BDTMitFisher = getMVA(calc_BDTMitFisher, MVA_Selections::BDTMitFisher);
    }

    void PrintStackedPlots(bool isBlind, bool drawRatio)
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        const std::string ratioCondition = drawRatio ? "_ratio" : "_noRatio";
        root_ext::PdfPrinter printer(outputFileName + blindCondition + ratioCondition + ".pdf");

        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            for (const HistogramDescriptor& hist : histograms) {
                //root_ext::PdfPrinter printer(outputFileName + blindCondition + "_" + hist.name + ratioCondition +".pdf");
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;
                StackedPlotDescriptor stackDescriptor(hist, ss_title.str(),false,ChannelNameLatex(),drawRatio);

                for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                    if(!category->draw) continue;

                    TH1D* histogram = GetSignalHistogram(eventCategory, category->name, hist.name);
                    if(!histogram) continue;

                    if(category->IsSignal())
                        stackDescriptor.AddSignalHistogram(histogram, category->title, category->color, category->draw_sf);
                    else if(category->IsBackground())
                        stackDescriptor.AddBackgroundHistogram(histogram, category->title, category->color);
                    else if(category->IsData())
                        stackDescriptor.AddDataHistogram(histogram, category->title, isBlind, GetBlindRegion(hist.name));
                }

                printer.PrintStack(stackDescriptor);
            }
        }
    }

    std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale) const
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_" << channel_name;
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else
            throw exception("Unsupported event energy scale ") << eventEnergyScale;
        full_name << "_8TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp)
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }


    void ProduceFileForLimitsCalculation(const std::string& hist_name,
                                         const root_ext::SmartHistogram<TH1D>& emptyDatacard)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" },
            { EventCategory::TwoJets_ZeroLooseBtag, "2jetloose0tag" },
            { EventCategory::TwoJets_OneLooseBtag, "2jetloose1tag" },
            { EventCategory::TwoJets_TwoLooseBtag, "2jetloose2tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        const std::string file_name = outputFileName + "_" + hist_name + ".root";
        std::shared_ptr<TFile> outputFile(new TFile(file_name.c_str(), "RECREATE"));
        outputFile->cd();
        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '") << dataCategory->name << "'.";
                for(const EventEnergyScale& eventEnergyScale : AllEventEnergyScales) {
                    const std::string full_hist_name = FlatAnalyzerData::FullHistogramName(hist_name, eventEnergyScale);
                    std::shared_ptr<TH1D> hist;
                    if(auto hist_orig = GetSignalHistogram(eventCategory, dataCategory->name, full_hist_name))
                        hist = std::shared_ptr<TH1D>(static_cast<TH1D*>(hist_orig->Clone()));
                    else {
                        std::cout << "Warning - Datacard histogram '" << full_hist_name
                                  << "' not found for data category '" << dataCategory->name << "' for eventCategory '"
                                  << categoryToDirectoryNameSuffix.at(eventCategory) << ".\n";

                        hist = std::shared_ptr<TH1D>(static_cast<TH1D*>(emptyDatacard.Clone()));
                    }
                    const std::string full_datacard_name = FullDataCardName(dataCategory->datacard, eventEnergyScale);
                    hist->Scale(dataCategory->limits_sf);
                    hist->Write(full_datacard_name.c_str());

                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "ZL") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = dataCategory->datacard + "_CMS_htt_" + dataCategory->datacard
                                + "Scale_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_up->Scale(1.02);
                        hist_syst_up->Write(name_syst_up.c_str());
                        std::shared_ptr<TH1D> hist_syst_down(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_down->Scale(0.98);
                        hist_syst_down->Write(name_syst_down.c_str());
                    }
                }
            }
        }
        outputFile->Close();
    }

    void SubtractBackgroundHistograms(TH1D& histogram, EventCategory eventCategory, EventRegion eventRegion,
                                      const std::string& current_category, bool verbose = false)
    {
        static const double correction_factor = 0.0000001;
        if(verbose)
            std::cout << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
                      << " for Event category '" << eventCategory << "'.\n"
                      << " for data category '" << current_category << "'.\n"
                      << "Initial integral: " << Integral(histogram, false) << ".\n";
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract ) continue;

            if(verbose)
                std::cout << "Sample '" << category->name << "': ";
            if(auto other_histogram = GetHistogram(eventCategory, category->name, eventRegion, histogram.GetName())) {
                histogram.Add(other_histogram, -1);
                if(verbose)
                    std::cout << Integral(*other_histogram, false) << ".\n";
            } else if(verbose)
                std::cout << "not found.\n";
        }

        const PhysicalValue original_Integral = Integral(histogram, false);
        if(verbose)
            std::cout << "Integral after bkg subtraction: " << original_Integral.value << ".\n" << std::endl;
        if (original_Integral.value < 0)
            throw exception("Integral after bkg subtraction is negative for histogram '")
                << histogram.GetName() << "' in event category " << eventCategory << " for event region " << eventRegion
                << ".";

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n){
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            std::cout << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                      << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                      << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                      << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(error*error + histogram.GetBinError(n)*histogram.GetBinError(n));
            histogram.SetBinContent(n,correction_factor);
            histogram.SetBinError(n,new_error);
        }
        const PhysicalValue new_Integral = Integral(histogram, false);
        const PhysicalValue correctionSF = original_Integral/new_Integral;
        histogram.Scale(correctionSF.value);
    }

    PhysicalValue CalculateFullIntegral(analysis::EventCategory eventCategory, analysis::EventRegion eventRegion,
                                        const std::string& hist_name, const DataCategoryPtrSet& dataCategories)
    {
        PhysicalValue integral;
        bool hist_found = false;

        for(const auto& dataCategory : dataCategories){
            if(auto hist_mc = GetHistogram(eventCategory, dataCategory->name, eventRegion, hist_name)) {
                hist_found = true;
                integral += Integral(*hist_mc, false);
            }
        }

        if(!hist_found)
            throw exception("No histogram with name '") << hist_name <<
                                                           "' was found in the given data category set, eventCategory: "
                                                        << eventCategory << " and eventRegion: " << eventRegion <<
                                                           "to calculate full integral.";

        return integral;
    }

private:

    static const std::vector< std::pair<double, double> >& GetBlindRegion(const std::string& hist_name)
    {
        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
            { { 100, 150 } },
            { { 200, 400 } },
            { { 100, 150 }, { 450, 500 }, { 800, 850 }, { 1150, 1200 }, { 1500, 1550 } }
        };

        static const std::map<std::string, size_t> histogramsToBlind = {
            { "m_sv", 1 }, { "m_vis", 1 }, { "m_bb", 1 }, { "m_ttbb", 2 }, { "m_ttbb_nomet", 2 },{ "m_ttbb_kinfit", 2 },
            { "m_ttbb_kinfit_only", 2 }, { "m_ttbb_kinfit_massCut", 2 }, { "m_ttbb_kinfit_only_massCut", 2 },
            { "m_bb_slice", 3 }
        };

        const auto findRegionId = [&]() -> size_t {
            if(histogramsToBlind.count(hist_name))
                return histogramsToBlind.at(hist_name);
            for(EventEnergyScale eventEnergyScale : AllEventEnergyScales) {
                for(const auto& entry : histogramsToBlind) {
                    const auto full_hist_name = FlatAnalyzerData::FullHistogramName(entry.first, eventEnergyScale);
                    if(full_hist_name == hist_name)
                        return entry.second;
                }
            }
            return 0;
        };

        const size_t regionId = findRegionId();
        if(regionId >= blindingRegions.size())
            throw analysis::exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");

        for(const HistogramDescriptor& hist : histograms) {
            if (hist.name != ReferenceHistogramName() && hist.name != "m_ttbb_kinfit_only_Central" &&
                    hist.name != "m_ttbb_kinfit_only_massCut_Central" ) continue;
            //PrintTables(of, sep, hist, false,false);
            //PrintTables(of, sep, hist, true, false);
            PrintTables(of, sep, hist, false,true);
            //PrintTables(of, sep, hist, true, true);
        }
        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const HistogramDescriptor& hist, bool includeOverflow,
                        bool includeError)
    {
        of << std::wstring(hist.title.begin(), hist.title.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (EventCategory eventCategory : EventCategoriesToProcess())
            of << eventCategory << sep;
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                if( TH1D* histogram = GetSignalHistogram(eventCategory, dataCategory->name, hist.name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    of << integral.ToString<wchar_t>(includeError) << sep;
                }
                else
                    of << "not found" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }

    void ProcessCompositDataCategories(EventCategory eventCategory, const std::string& hist_name)
    {
        for(const DataCategory* composit : dataCategoryCollection.GetCategories(DataCategoryType::Composit)) {
            for(const std::string& sub_name : composit->sub_categories) {
                const DataCategory& sub_category = dataCategoryCollection.FindCategory(sub_name);
                auto sub_hist = GetSignalHistogram(eventCategory, sub_category.name, hist_name);
                if(!sub_hist) continue;

                if(auto composit_hist = GetSignalHistogram(eventCategory, composit->name, hist_name))
                    composit_hist->Add(sub_hist);
                else
                    CloneSignalHistogram(eventCategory, composit->name, *sub_hist);
            }
        }
    }

protected:
    std::string inputPath;
    std::string outputFileName;
    DataCategoryCollection dataCategoryCollection;
    std::vector<HistogramDescriptor> histograms;
    FullAnaData fullAnaData;
};

} // namespace analysis
