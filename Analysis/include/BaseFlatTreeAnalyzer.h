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
#include "AnalysisBase/include/FlatTree.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"
#include "KinFit.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

static const std::vector<double> mass_bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140., 150,
                                               160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };

static const std::vector<double> mass_bins_slice_2fette = { 
0,	10,	20,	30,	40,	50,	60,	70,	80,	90,	100,	110,	120,	130,
	140,	150,	160,	170,	180,	190,	200,	225,	250,	275,	300,	325,	350,
	360,	370,	380,	390,	400,	410,	420,	430,	440,	450,	460,	470,	480,
	490,	500,	510,	520,	530,	540,	550,	575,	600,	625,	650,	675,	700
};

static const std::vector<double> mass_bins_slice_5fette_lb = { 
0,	20,	40,	60,	80,	100,	120,	140,	160,	180,	200,	250,	300,	350,
	370,	390,	410,	430,	450,	470,	490,	510,	530,	550,	600,	650,	700,
	720,	740,	760,	780,	800,	820,	840,	860,	880,	900,	950,	1000,	1050,
	1070,	1090,	1110,	1130,	1150,	1170,	1190,	1210,	1230,	1250,	1300,	1350,	1400,
	1420,	1440,	1460,	1480,	1500,	1520,	1540,	1560,	1580,	1600,	1650,	1700,	1750
};

static const std::vector<double> mass_bins_slice_5fette_fb = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140., 150,
                                               160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350,
360,	370,	380,	390,	400,	410,	420,	430,	440,	450,	460,	470,	480,	490,	500,	510,	520,	530,	540,	550,	575,	600,	625,	650,	675,	700,
710,	720,	730,	740,	750,	760,	770,	780,	790,	800,	810,	820,	830,	840,	850,	860,	870,	880,	890,	900,	925,	950,	975,	1000,	1025,	1050,
1060,	1070,	1080,	1090,	1100,	1110,	1120,	1130,	1140,	1150,	1160,	1170,	1180,	1190,	1200,	1210,	1220,	1230,	1240,	1250,	1275,	1300,	1325,	1350,	1375,	1400,
1410,	1420,	1430,	1440,	1450,	1460,	1470,	1480,	1490,	1500,	1510,	1520,	1530,	1540,	1550,	1560,	1570,	1580,	1590,	1600,	1625,	1650,	1675,	1700,	1725,	1750,
 };

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY(pt_1, 50, 0, 200)
    TH1D_ENTRY(eta_1, 60, -3, 3)
    TH1D_ENTRY(pt_2, 50, 0, 200)
    TH1D_ENTRY(eta_2, 60, -3, 3)
    TH1D_ENTRY(pt_b1, 50, 0, 200)
    TH1D_ENTRY(eta_b1, 60, -3, 3)
    TH1D_ENTRY(pt_b2, 50, 0, 200)
    TH1D_ENTRY(eta_b2, 60, -3, 3)
    TH1D_ENTRY(pt_H_tt, 35, 0, 350)
    TH1D_ENTRY(pt_H_bb, 35, 0, 350)
    TH1D_ENTRY(pt_H_hh, 35, 0, 350)
    TH1D_ENTRY_CUSTOM(m_sv, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_up, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_down, mass_bins)
    TH1D_ENTRY_CUSTOM(m_vis, mass_bins)
    TH1D_ENTRY(m_bb, 30, 0, 600)
    TH1D_ENTRY(m_ttbb, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit_up, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit_down, 50, 0, 1000)
    TH1D_ENTRY_CUSTOM(m_bb_slice, mass_bins_slice_5fette_fb)
    TH1D_ENTRY_CUSTOM(m_bb_slice_up, mass_bins_slice_5fette_fb)
    TH1D_ENTRY_CUSTOM(m_bb_slice_down, mass_bins_slice_5fette_fb)
    TH1D_ENTRY(m_ttbb_nomet, 100, 0, 1000)
    TH1D_ENTRY(pt_ttbb_nomet, 35, 0, 350)
    TH1D_ENTRY(DeltaPhi_tt, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_bb, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_bb_MET, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_tt_MET, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_hh, 80, -4, 4)
    TH1D_ENTRY(DeltaR_tt, 60, 0, 6)
    TH1D_ENTRY(DeltaR_bb, 60, 0, 6)
    TH1D_ENTRY(DeltaR_hh, 60, 0, 6)
    TH1D_ENTRY(MVA_BDT, 40, -1, 1)
    TH1D_ENTRY(MVA_BDTD, 40, -1, 1)
    TH1D_ENTRY(MVA_BDTMitFisher, 40, -1, 1)
    TH1D_ENTRY(mt_1, 50, 0, 50)
    TH1D_ENTRY(mt_2, 50, 0, 300)
    TH1D_ENTRY(pt_H_tt_MET, 35, 0, 350)
};

class BaseFlatTreeAnalyzer {
public:
    typedef std::map<EventType_QCD, FlatAnalyzerData> AnaDataQCD;
    typedef std::map<EventType_Wjets, FlatAnalyzerData> AnaDataWjets;

    struct AnaDataEventTypes {
        AnaDataQCD QCD;
        AnaDataWjets Wjets;
    };

    typedef std::map<std::string, AnaDataEventTypes> AnaDataForDataCategory;
    typedef std::map<EventCategory, AnaDataForDataCategory> FullAnaData;

    BaseFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& _outputFileName, const std::string& _signalName,
                         const std::string& _dataName, const std::string& _mvaXMLpath, bool _WjetsData = false,
                         bool _isBlind=false)
        : inputPath(_inputPath), signalName(_signalName), dataName(_dataName), outputFileName(_outputFileName),
          WjetsData(_WjetsData), isBlind(_isBlind)
    {
        TH1::SetDefaultSumw2();

        categories = DataCategory::ReadFromFile(source_cfg, dataName);
        histograms = HistogramDescriptor::ReadFromFile(hist_cfg);
    }

    void Run()
    {
        std::cout << "Processing sources... " << std::endl;
        for(DataCategory& category : categories) {
            std::cout << category << std::endl;
            //if (!category.IsSignal()) continue;
            for (DataSource& source : category.sources){
                const std::string fullFileName = inputPath + "/" + source.file_name;
                source.file = new TFile(fullFileName.c_str(), "READ");
                if(source.file->IsZombie()) {
                    std::ostringstream ss;
                    ss << "Input file '" << source.file_name << "' not found.";
                    throw std::runtime_error(ss.str());
                }
                source.tree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree(*source.file, "flatTree"));
                ProcessDataSource(category, source);
            }
        }

        std::cout << "Estimating QCD, Wjets and bkg sum... " << std::endl;
        for (auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const auto& hist : histograms) {
                EstimateQCD(eventCategory, anaData, hist);
                if (WjetsData) EstimateWjets(eventCategory, anaData, hist);

                EstimateSumBkg(eventCategory, anaData, hist);
            }
        }

        std::cout << "Saving tables and printing stacked plots... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");
        ProduceFileForLimitsCalculation();
        PrintStackedPlots();
    }

protected:
    void ProcessDataSource(const DataCategory& dataCategory, const DataSource& dataSource)
    {
        const analysis::DataCategory& Ztautau = FindCategory("LIMITS Ztautau");
        const analysis::DataCategory& DYJets = FindCategory("DYJets");

        for(size_t current_entry = 0; current_entry < dataSource.tree->GetEntries(); ++current_entry) {
            dataSource.tree->GetEntry(current_entry);
            const ntuple::Flat& event = dataSource.tree->data;
            const EventType_QCD eventTypeQCD = DetermineEventTypeForQCD(event);
            const EventType_Wjets eventTypeWjets = DetermineEventTypeForWjets(event);
            const EventCategoryVector eventCategories = DetermineEventCategories(event);
            const double weight = dataCategory.IsData() ? 1 : event.weight * dataSource.scale_factor;
            if(eventTypeQCD == EventType_QCD::Unknown) continue;
            for(auto eventCategory : eventCategories) {
//                if(!PassMvaCut(event, eventCategory)) continue;
                FillHistograms(fullAnaData[eventCategory][dataCategory.name].QCD[eventTypeQCD], event, weight, eventCategory);
                FillHistograms(fullAnaData[eventCategory][dataCategory.name].Wjets[eventTypeWjets], event, weight, eventCategory);
                if(dataCategory.name == DYJets.name && std::abs(event.pdgId_2_MC) == particles::tau.RawCode())
                    FillHistograms(fullAnaData[eventCategory][Ztautau.name].QCD[eventTypeQCD], event, weight, eventCategory);
            }
        }
    }

    virtual const std::string& ChannelName() = 0;

    virtual EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) = 0;
    virtual EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) = 0;
    virtual bool PassMvaCut(const ntuple::Flat& event, EventCategory eventCategory) = 0;

    virtual EventCategoryVector DetermineEventCategories(const ntuple::Flat& event)
    {
        EventCategoryVector categories;
        categories.push_back(EventCategory::Inclusive);

        std::vector<Float_t> goodCVSvalues;
        for (unsigned i = 0; i < event.eta_Bjets.size(); ++i){
            if ( std::abs(event.eta_Bjets.at(i)) >= cuts::Htautau_Summer13::btag::eta) continue;
            goodCVSvalues.push_back(event.csv_Bjets.at(i));
        }

        if (goodCVSvalues.size() == 1) {
            if (goodCVSvalues.at(0) <= cuts::Htautau_Summer13::btag::CSVT )
                categories.push_back(EventCategory::OneJet_ZeroBtag);
            else
                categories.push_back(EventCategory::OneJet_OneBtag);
        }

        if (goodCVSvalues.size() >= 2) {
            if (goodCVSvalues.at(0) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_ZeroBtag);
            else if ( goodCVSvalues.at(1) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_OneBtag);
            else
                categories.push_back(EventCategory::TwoJets_TwoBtag);
        }
        return categories;
    }

    virtual void EstimateQCD(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                             const HistogramDescriptor& hist) = 0;
    virtual void EstimateWjets(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                               const HistogramDescriptor& hist) = 0;

    virtual void FillHistograms(FlatAnalyzerData& anaData, const ntuple::Flat& event, double weight,
                                EventCategory eventCategory)
    {
        anaData.m_sv().Fill(event.m_sv, weight);
        anaData.m_sv_up().Fill(event.m_sv_Up, weight);
        anaData.m_sv_down().Fill(event.m_sv_Down, weight);
        anaData.pt_1().Fill(event.pt_1, weight);
        anaData.eta_1().Fill(event.eta_1, weight);
        anaData.pt_2().Fill(event.pt_2, weight);
        anaData.eta_2().Fill(event.eta_2, weight);
        TLorentzVector first_cand, second_cand;
        first_cand.SetPtEtaPhiE(event.pt_1,event.eta_1,event.phi_1,event.energy_1);
        second_cand.SetPtEtaPhiE(event.pt_2,event.eta_2,event.phi_2,event.energy_2);
        anaData.DeltaPhi_tt().Fill(first_cand.DeltaPhi(second_cand), weight);
        anaData.DeltaR_tt().Fill(first_cand.DeltaR(second_cand), weight);
        TLorentzVector Htt = first_cand + second_cand;
        anaData.pt_H_tt().Fill(Htt.Pt(),weight);
        anaData.m_vis().Fill(Htt.M(),weight);
        TLorentzVector MET;
        MET.SetPtEtaPhiM(event.mvamet,0,event.mvametphi,0);
        TMatrixD metcov(2,2);
        metcov(0,0)=event.mvacov00;
        metcov(1,0)=event.mvacov10;
        metcov(0,1)=event.mvacov01;
        metcov(1,1)=event.mvacov11;
        TLorentzVector Htt_MET = Htt + MET;
        anaData.pt_H_tt_MET().Fill(Htt_MET.Pt(), weight);
        anaData.DeltaPhi_tt_MET().Fill(Htt.DeltaPhi(MET), weight);
        anaData.mt_1().Fill(event.mt_1, weight);
        anaData.mt_2().Fill(event.mt_2, weight);
        if(event.mass_Bjets.size() >= 2) {
            std::vector<TLorentzVector> b_momentums(2);
            for(size_t n = 0; n < b_momentums.size(); ++n)
                b_momentums.at(n).SetPtEtaPhiM(event.pt_Bjets.at(n), event.eta_Bjets.at(n), event.phi_Bjets.at(n),
                                               event.mass_Bjets.at(n));
            anaData.pt_b1().Fill(b_momentums.at(0).Pt(), weight);
            anaData.eta_b1().Fill(b_momentums.at(0).Eta(), weight);
            anaData.pt_b2().Fill(b_momentums.at(1).Pt(), weight);
            anaData.eta_b2().Fill(b_momentums.at(1).Eta(), weight);
            anaData.DeltaPhi_bb().Fill(b_momentums.at(0).DeltaPhi(b_momentums.at(1)), weight);
            anaData.DeltaR_bb().Fill(b_momentums.at(0).DeltaR(b_momentums.at(1)), weight);
            const TLorentzVector Hbb = b_momentums.at(0) + b_momentums.at(1);
            anaData.pt_H_bb().Fill(Hbb.Pt(),weight);
            anaData.m_bb().Fill(Hbb.M(), weight);
            anaData.DeltaPhi_bb_MET().Fill(Hbb.DeltaPhi(MET), weight);
            anaData.DeltaPhi_hh().Fill(Htt.DeltaPhi(Hbb), weight);
            anaData.DeltaR_hh().Fill(Htt.DeltaR(Hbb), weight);
            const TLorentzVector Candidate_ttbb = Hbb + Htt + MET;
            anaData.m_ttbb().Fill(Candidate_ttbb.M(), weight);
            anaData.pt_H_hh().Fill(Candidate_ttbb.Pt(), weight);
            const TLorentzVector Candidate_ttbb_noMET = Hbb + Htt;
            anaData.pt_ttbb_nomet().Fill(Candidate_ttbb_noMET.Pt(), weight);
            anaData.m_ttbb_nomet().Fill(Candidate_ttbb_noMET.M(), weight);
//            const double m_ttbb_kinFit = analysis::CorrectMassByKinfit(b_momentums.at(0),b_momentums.at(1),first_cand,second_cand,MET,metcov);
            const double m_ttbb_kinFit = -1;
            anaData.m_ttbb_kinfit().Fill(m_ttbb_kinFit,weight);
            anaData.m_ttbb_kinfit_up().Fill(1.04*m_ttbb_kinFit,weight);
            anaData.m_ttbb_kinfit_down().Fill(0.96*m_ttbb_kinFit,weight);

		if(Hbb.M() < 60)
		{
			anaData.m_bb_slice().Fill(event.m_sv,weight);
			anaData.m_bb_slice_up().Fill(event.m_sv_Up,weight);
			anaData.m_bb_slice_down().Fill(event.m_sv_Down,weight);
		}
		else if (60 <= Hbb.M() && Hbb.M() < 100)
		{
			anaData.m_bb_slice().Fill(event.m_sv+350,weight);
			anaData.m_bb_slice_up().Fill(event.m_sv_Up+350,weight);
			anaData.m_bb_slice_down().Fill(event.m_sv_Down+350,weight);
		}
		else if (100 <= Hbb.M() && Hbb.M() < 140)
		{
			anaData.m_bb_slice().Fill(event.m_sv+700,weight);
			anaData.m_bb_slice_up().Fill(event.m_sv_Up++700,weight);
			anaData.m_bb_slice_down().Fill(event.m_sv_Down+700,weight);
		}
		else if (140 <= Hbb.M() && Hbb.M() < 200)
		{
			anaData.m_bb_slice().Fill(event.m_sv+1050,weight);
			anaData.m_bb_slice_up().Fill(event.m_sv_Up+1050,weight);
			anaData.m_bb_slice_down().Fill(event.m_sv_Down+1050,weight);
		}
		else if (200 <= Hbb.M() && Hbb.M() <= 600)
		{
			anaData.m_bb_slice().Fill(event.m_sv+1400,weight);
			anaData.m_bb_slice_up().Fill(event.m_sv_Up+1400,weight);
			anaData.m_bb_slice_down().Fill(event.m_sv_Down+1400,weight);
		}

        }
    }

    void PrintStackedPlots()
    {
        root_ext::PdfPrinter printer(outputFileName + ".pdf");
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const HistogramDescriptor& hist : histograms) {
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;
                StackedPlotDescriptor stackDescriptor(hist, ss_title.str());

                for(const DataCategory& category : categories) {
                    TH1D* histogram;
                    if(!(histogram = anaData[category.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) continue;
                    if(category.IsReference() || /*(category.IsSignal() && !category.NameContains(signalName)) || */
                            category.IsSumBkg() || category.IsForLimitsOnly()) continue;
                    if(category.IsData())
                        stackDescriptor.AddDataHistogram(histogram, category.title, isBlind, GetBlindRegion(hist.name));
                    else if(category.IsSignal()){
                        if (eventCategory == EventCategory::TwoJets_TwoBtag){
                            if (category.name == "SIGNAL ggHhh300")
                            stackDescriptor.AddSignalHistogram(histogram, "1x hh#rightarrow#tau#taubb(m_{H}=300 tan#beta=2.5)", category.color, 1);
                            if (category.name == "SIGNAL Radion500")
                            stackDescriptor.AddSignalHistogram(histogram, "1x Radion#rightarrowhh#rightarrow#tau#taubb(m_{H}=500)", category.color, 1);
                            if (category.name == "SIGNAL Graviton500")
                            stackDescriptor.AddSignalHistogram(histogram, "1x Graviton#rightarrowhh#rightarrow#tau#taubb(m_{H}=500)", category.color, 1);
                            if (category.name == "SIGNAL Radion1000")
                            stackDescriptor.AddSignalHistogram(histogram, "1x Radion#rightarrowhh#rightarrow#tau#taubb(m_{H}=1000)", category.color, 1);
                            if (category.name == "SIGNAL Graviton1000")
                            stackDescriptor.AddSignalHistogram(histogram, "1x Graviton#rightarrowhh#rightarrow#tau#taubb(m_{H}=1000)", category.color, 1);
                        }
                        else
                            stackDescriptor.AddSignalHistogram(histogram, category.title, category.color, 10);
                    }
                    else
                        stackDescriptor.AddBackgroundHistogram(histogram, category.title, category.color);
                }
                printer.PrintStack(stackDescriptor);
            }
        }
    }

    const analysis::DataCategory& FindCategory(const std::string& prefix) const
    {
        for (const analysis::DataCategory& category : categories){
            if (category.NameContains(prefix))
                return category;
        }
        throw analysis::exception("category not found : ") << prefix ;
    }

    void ProduceFileForLimitsCalculation()
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0tag" },
            { EventCategory::OneJet_OneBtag, "1jet1tag" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::vector< std::pair<std::string, std::string> > dataCategoriesForLimits = {
            { "VIRTUAL QCD", "QCD" }, { "TTbar", "TT" }, { "LIMITS VH125", "VH125" },
            { "LIMITS DiBoson", "VV" }, { "LIMITS Wjets", "W" }, { "ZJ", "ZJ" }, { "ZL", "ZL" },
            { "ZLL", "ZLL" }, { "LIMITS Ztautau", "ZTT" }, { "LIMITS bbH100", "bbH100" }, { "LIMITS bbH110", "bbH110" },
            { "LIMITS bbH120", "bbH120" }, { "LIMITS bbH130", "bbH130" }, { "LIMITS bbH140", "bbH140" }, { "LIMITS bbH160", "bbH160" },
            { "LIMITS bbH180", "bbH180" }, { "LIMITS bbH200", "bbH200" }, { "LIMITS bbH250", "bbH250" }, { "LIMITS bbH300", "bbH300" },
            { "LIMITS bbH350", "bbH350" }, { "LIMITS bbH400", "bbH400" }, { "LIMITS bbH90", "bbH90" },
            { "DATA Tau", "data_obs" }, { "DATA TauPlusX", "data_obs" },
            { "LIMITS ggAToZhToLLBB250", "ggAToZhToLLBB250" },
            { "LIMITS ggAToZhToLLBB260", "ggAToZhToLLBB260" }, { "LIMITS ggAToZhToLLBB270", "ggAToZhToLLBB270" },
            { "LIMITS ggAToZhToLLBB280", "ggAToZhToLLBB280" }, { "LIMITS ggAToZhToLLBB290", "ggAToZhToLLBB290" },
            { "LIMITS ggAToZhToLLBB300", "ggAToZhToLLBB300" }, { "LIMITS ggAToZhToLLBB310", "ggAToZhToLLBB310" },
            { "LIMITS ggAToZhToLLBB320", "ggAToZhToLLBB320" }, { "LIMITS ggAToZhToLLBB330", "ggAToZhToLLBB330" },
            { "LIMITS ggAToZhToLLBB340", "ggAToZhToLLBB340" }, { "LIMITS ggAToZhToLLBB350", "ggAToZhToLLBB350" },
            { "LIMITS ggAToZhToLLTauTau260", "ggAToZhToLLTauTau260" }, { "LIMITS ggAToZhToLLTauTau270", "ggAToZhToLLTauTau270" },
            { "LIMITS ggAToZhToLLTauTau280", "ggAToZhToLLTauTau280" }, { "LIMITS ggAToZhToLLTauTau290", "ggAToZhToLLTauTau290" },
            { "LIMITS ggAToZhToLLTauTau300", "ggAToZhToLLTauTau300" }, { "LIMITS ggAToZhToLLTauTau310", "ggAToZhToLLTauTau310" },
            { "LIMITS ggAToZhToLLTauTau320", "ggAToZhToLLTauTau320" }, { "LIMITS ggAToZhToLLTauTau330", "ggAToZhToLLTauTau330" },
            { "LIMITS ggAToZhToLLTauTau340", "ggAToZhToLLTauTau340" }, { "LIMITS ggAToZhToLLTauTau350", "ggAToZhToLLTauTau350" },
            { "LIMITS ggH125", "ggH125" }, { "LIMITS qqH125", "qqH125" },
            { "LIMITS ggHhh260", "ggHTohhTo2Tau2B260" }, { "LIMITS ggHhh270", "ggHTohhTo2Tau2B270" },
            { "LIMITS ggHhh280", "ggHTohhTo2Tau2B280" }, { "LIMITS ggHhh290", "ggHTohhTo2Tau2B290" },
            { "LIMITS ggHhh300", "ggHTohhTo2Tau2B300" }, { "LIMITS ggHhh310", "ggHTohhTo2Tau2B310" },
            { "LIMITS ggHhh320", "ggHTohhTo2Tau2B320" }, { "LIMITS ggHhh330", "ggHTohhTo2Tau2B330" },
            { "LIMITS ggHhh340", "ggHTohhTo2Tau2B340" }, { "LIMITS ggHhh350", "ggHTohhTo2Tau2B350" },
            { "LIMITS Radion300", "ggRadionTohhTo2Tau2B300" }, { "LIMITS Radion500", "ggRadionTohhTo2Tau2B500" },
            { "LIMITS Radion700", "ggRadionTohhTo2Tau2B700" }, { "LIMITS Radion1000", "ggRadionTohhTo2Tau2B1000" },
            { "LIMITS Graviton270", "ggGravitonTohhTo2Tau2B270" }, { "LIMITS Graviton300", "ggGravitonTohhTo2Tau2B300" },
            { "LIMITS Graviton500", "ggGravitonTohhTo2Tau2B500" }, { "LIMITS Graviton700", "ggGravitonTohhTo2Tau2B700" },
            { "LIMITS Graviton1000", "ggGravitonTohhTo2Tau2B1000" },
            // Debug
            { "LIMITS TTbar_hadr", "TT_hadr" }, { "LIMITS TTbar_lept", "TT_lept" }, { "LIMITS TTbar_semi", "TT_semi" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);

        std::shared_ptr<TFile> outputFile(new TFile((outputFileName + ".root").c_str(), "RECREATE"));
        outputFile->cd();
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaDataForCategory = fullAnaDataEntry.second;
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const auto& limitDataCategory : dataCategoriesForLimits) {
                if(!anaDataForCategory.count(limitDataCategory.first)) {
                    std::cerr << "WARNING: Dataset '" << limitDataCategory.first << "' for event category '"
                              << eventCategory << "' not found." << std::endl;
                    continue;
                }
                FlatAnalyzerData& anaData = anaDataForCategory[limitDataCategory.first].QCD[EventType_QCD::OS_Isolated];
//		anaData.m_sv().Write(limitDataCategory.second.c_str());
//		anaData.m_ttbb_kinfit().Write(limitDataCategory.second.c_str());
		anaData.m_bb_slice().Write(limitDataCategory.second.c_str());
                const std::string namePrefix = limitDataCategory.second + "_CMS_scale_t_" + channel_name + "_8TeV";
                const std::string nameDown = namePrefix + "Down";
                const std::string nameUp = namePrefix + "Up";

//		anaData.m_ttbb_kinfit_up().Write(nameUp.c_str());
//		anaData.m_ttbb_kinfit_down().Write(nameDown.c_str());
//		anaData.m_sv_down().Write(nameDown.c_str());
//		anaData.m_sv_up().Write(nameUp.c_str());
		anaData.m_bb_slice_up().Write(nameUp.c_str());
		anaData.m_bb_slice_down().Write(nameDown.c_str());
            }
        }
        outputFile->Close();
    }

private:

    static const std::pair<double, double>& GetBlindRegion(const std::string& hist_name)
    {
        static const std::vector< std::pair<double, double> > blindingRegions = {
            { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() }, { 100, 150 }, { 250, 350 }
        };
        static const std::map<std::string, size_t> histogramsToBlind = {
            { "m_sv", 1 }, { "m_sv_up", 1 }, { "m_sv_down", 1 }, { "m_vis", 1 }, { "m_bb", 1 },
            { "m_ttbb", 2 }, { "m_ttbb_nomet", 2 },
	    { "m_ttbb_kinfit", 2 }, { "m_ttbb_kinfit_up", 2 }, { "m_ttbb_kinfit_down", 2 },
	    { "m_bb_slice", 2}, { "m_bb_slice_up", 2}, { "m_bb_slice_down", 2}
        };

        if(!histogramsToBlind.count(hist_name)) return blindingRegions.at(0);
        const size_t regionId = histogramsToBlind.at(hist_name);
        if(regionId >= blindingRegions.size())
            throw analysis::exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");

        for(const HistogramDescriptor& hist : histograms) {
            if (hist.name != "m_sv") continue;
            PrintHistogram(of, sep, hist, false,false);
            PrintHistogram(of, sep, hist, true, false);
            PrintHistogram(of, sep, hist, false,true);
            PrintHistogram(of, sep, hist, true, true);
        }
        of.flush();
        of.close();
    }

    void PrintHistogram(std::wostream& of, const std::wstring& sep, const HistogramDescriptor& hist, bool includeOverflow,
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

        for (const auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            of << eventCategory << sep;
        }
        of << std::endl;

        for (DataCategory& dataCategory : categories) {
            if (dataCategory.IsReference()) continue;
            of << std::wstring(dataCategory.title.begin(), dataCategory.title.end()) << sep;
            for (auto& fullAnaDataEntry : fullAnaData) {
                AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
                if( TH1D* histogram = anaData[dataCategory.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name) ) {
                    typedef std::pair<Int_t, Int_t> limit_pair;
                    const limit_pair limits = includeOverflow ? limit_pair(0, histogram->GetNbinsX() + 1)
                                                              : limit_pair(1, histogram->GetNbinsX());
                    double error = 0.;
                    const double integral = histogram->IntegralAndError(limits.first, limits.second, error);
                    of << MakeStringRepresentationForValueWithError(integral, error, includeError) << sep;
                }
                else
                    of << "NaN" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;

    }

    void EstimateSumBkg(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                             const analysis::HistogramDescriptor& hist)
    {
        const analysis::DataCategory& sumBkg = FindCategory("SUM");
        const analysis::DataCategory& ewk = FindCategory("EWK");
        root_ext::SmartHistogram<TH1D>* hist_ewk;
        if(!(hist_ewk = anaData[ewk.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) return;
        TH1D& histogram = anaData[sumBkg.name].QCD[EventType_QCD::OS_Isolated].Clone(*hist_ewk);
        for (const analysis::DataCategory& category : categories){
            if (category.IsReference() || category.IsData() || category.IsSignal() || category.name == ewk.name
                    || category.IsSumBkg() || category.IsForLimitsOnly()) continue;

            TH1D* otherBkg_hist;
            if(!(otherBkg_hist = anaData[category.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) continue;
            histogram.Add(otherBkg_hist);
        }
    }

protected:
    std::string inputPath;
    std::string signalName, dataName;
    std::string outputFileName;
    DataCategoryCollection categories;
    std::vector<HistogramDescriptor> histograms;
    FullAnaData fullAnaData;
    bool WjetsData;
    bool isBlind;
};

} // namespace analysis
