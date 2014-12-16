/*!
 * \file HH_SyncTreeProducer.C
 * \brief Sync Tree Producer from FlatTree.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \author Maria Agnese Ciocci (Siena University, INFN Pisa)
 * \date 2014-11-13 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>,
 *                Maria Agnese Ciocci <mariaagnese.ciocci@pi.infn.it>
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

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"
#include "AnalysisBase/include/SyncTree.h"

class SyncTreeProducer : public analysis::LightBaseFlatTreeAnalyzer {
public:
    typedef std::map<analysis::EventEnergyScale, ntuple::Flat> ES_toEvent_Map;
    typedef std::map<analysis::EventId, ES_toEvent_Map> EventId_ToES_Map;

    SyncTreeProducer(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), inclusive(0), passed(0)
    {
        GetOutputFile().cd();
        syncTree = std::shared_ptr<ntuple::SyncTree>(new ntuple::SyncTree("syncTree"));
        recalc_kinfit = false;
    }

    virtual ~SyncTreeProducer() {
        syncTree->Write();
    }

protected:
    static bool PassSyncTreeSelection(const analysis::FlatEventInfo& eventInfo)
    {
        using analysis::EventRegion;
        const ntuple::Flat& event = *eventInfo.event;
        if (eventInfo.channel == analysis::Channel::MuTau){
            using namespace cuts::Htautau_Summer13::MuTau;

            return !(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= muonID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::ETau){
            using namespace cuts::Htautau_Summer13::ETau;

            return !(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= electronID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::TauTau){
            using namespace cuts::Htautau_Summer13::TauTau::tauID;

            return !(event.againstElectronLooseMVA_2 <= againstElectronLooseMVA3
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits);
//            if(!event.againstElectronLooseMVA_2
//                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
//                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
//                    || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
//                        && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
//                return false;

//            const bool os = event.q_1 * event.q_2 == -1;
//            const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
//                             event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

//            if(eventInfo.bjet_momentums.size() < 2 || !eventInfo.fitResults.has_valid_mass)
//                return false;

//            using namespace cuts::massWindow;
//            const double mass_tautau = eventInfo.event->m_sv_MC;
//            const double mass_bb = eventInfo.Hbb.M();
//            const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
//                    && mass_bb > m_bb_low && mass_bb < m_bb_high;

//            return os && iso && inside_mass_window;
        }
        throw analysis::exception("unsupported channel ") << eventInfo.channel;
    }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;

        if (category != EventCategory::Inclusive) return;
        ++inclusive;
        if (!PassSyncTreeSelection(eventInfo)) return;
        ++passed;

        const ntuple::Flat& event = *eventInfo.event;
        const analysis::EventId eventId(event.run,event.lumi,event.evt);
        eventId_ToES_Map[eventId][static_cast<analysis::EventEnergyScale>(event.eventEnergyScale)] = event;
    }

    virtual void EndOfRun() override
    {
        std::cout << "inclusive evt = " << inclusive << ", passed = " << passed << std::endl;
        std::cout << "eventId_ToES_Map.size = " << eventId_ToES_Map.size() << std::endl;
        for (const auto& event_iter : eventId_ToES_Map){
            const ES_toEvent_Map& es_toEventMap = event_iter.second;

            if (!es_toEventMap.count(analysis::EventEnergyScale::Central)) continue;
            const ntuple::Flat& event = es_toEventMap.at(analysis::EventEnergyScale::Central);
            analysis::FlatEventInfo eventInfo(event, analysis::FlatEventInfo::BjetPair(0, 1), false);
            static const double default_value = ntuple::DefaultFillValueForSyncTree();
            syncTree->run() = event.run;
            syncTree->lumi() = event.lumi;
            syncTree->evt() = event.evt;

            syncTree->npv() = event.npv;
            syncTree->npu() = event.npu;

            syncTree->puweight() = event.puweight;
            syncTree->trigweight_1() = event.trigweight_1;
            syncTree->trigweight_2() = event.trigweight_2;
            syncTree->idweight_1() = event.idweight_1;
            syncTree->idweight_2() = event.idweight_2;
            syncTree->isoweight_1() = event.isoweight_1;
            syncTree->isoweight_2() = event.isoweight_2;
            syncTree->fakeweight() = event.fakeweight_2;

            syncTree->weight() = event.weight;
            syncTree->embeddedWeight() = event.embeddedWeight;

            syncTree->mvis() = eventInfo.Htt.M();
            syncTree->m_sv() = event.m_sv_MC;
            syncTree->pt_sv() = event.pt_sv_MC;
            syncTree->eta_sv() = event.eta_sv_MC;
            syncTree->phi_sv() = event.phi_sv_MC;
            if (es_toEventMap.count(analysis::EventEnergyScale::TauUp)){
                const ntuple::Flat& event_up = es_toEventMap.at(analysis::EventEnergyScale::TauUp);
                syncTree->m_sv_Up() = event_up.m_sv_MC;
            }
            else
                syncTree->m_sv_Up() = default_value;
            if (es_toEventMap.count(analysis::EventEnergyScale::TauDown)){
                const ntuple::Flat& event_down = es_toEventMap.at(analysis::EventEnergyScale::TauDown);
                syncTree->m_sv_Down() = event_down.m_sv_MC;
            }
            else
                syncTree->m_sv_Down() = default_value;

            syncTree->pt_1() = event.pt_1;
            syncTree->phi_1() = event.phi_1;
            syncTree->eta_1() = event.eta_1;
            syncTree->m_1() = event.m_1;
            syncTree->q_1() = event.q_1;
            syncTree->mt_1() = event.mt_1;
            syncTree->d0_1() = event.d0_1;
            syncTree->dZ_1() = event.dZ_1;

            // MVA iso for hadronic Tau, Delta Beta for muon and electron
            syncTree->iso_1() = event.iso_1;
            syncTree->mva_1() = event.mva_1;
            syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
            syncTree->againstElectronMVA3raw_1() = event.againstElectronMVA3raw_1;
            syncTree->byIsolationMVA2raw_1() = event.byIsolationMVA2raw_1;
            syncTree->againstMuonLoose2_1() = event.againstMuonLoose_1;
            syncTree->againstMuonMedium2_1() = event.againstMuonMedium_1;
            syncTree->againstMuonTight2_1() = event.againstMuonTight_1;

            syncTree->pt_2() = event.pt_2;
            syncTree->phi_2() = event.phi_2;
            syncTree->eta_2() = event.eta_2;
            syncTree->m_2() = event.m_2;
            syncTree->q_2() = event.q_2;
            syncTree->mt_2() = event.mt_2;
            syncTree->d0_2() = event.d0_2;
            syncTree->dZ_2() = event.dZ_2;

            syncTree->iso_2() = event.iso_2;
            syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
            syncTree->againstElectronMVA3raw_2() = event.againstElectronMVA3raw_2;
            syncTree->byIsolationMVA2raw_2() = event.byIsolationMVA2raw_2;
            syncTree->againstMuonLoose2_2() = event.againstMuonLoose_2;
            syncTree->againstMuonMedium2_2() = event.againstMuonMedium_2;
            syncTree->againstMuonTight2_2() = event.againstMuonTight_2;

            syncTree->pt_tt() = eventInfo.Htt_MET.Pt();

            syncTree->met() = event.met;
            syncTree->metphi() = event.metphi;
            syncTree->metcov00() = event.metcov00;
            syncTree->metcov01() = event.metcov01;
            syncTree->metcov10() = event.metcov10;
            syncTree->metcov11() = event.metcov11;

            syncTree->mvamet() = event.mvamet;
            syncTree->mvametphi() = event.mvametphi;
            syncTree->mvacov00() = event.mvacov00;
            syncTree->mvacov01() = event.mvacov01;
            syncTree->mvacov10() = event.mvacov10;
            syncTree->mvacov11() = event.mvacov11;

            syncTree->njets() = event.njets;
            syncTree->njetspt20() = event.njetspt20;

            syncTree->jpt_1() = default_value;
            syncTree->jeta_1() = default_value;
            syncTree->jphi_1() = default_value;
            syncTree->jptraw_1() = default_value;
            //syncTree->jptunc_1() = default_value;
            syncTree->jmva_1() = default_value;
            syncTree->jpass_1() = default_value;
            syncTree->jpt_2() = default_value;
            syncTree->jeta_2() = default_value;
            syncTree->jphi_2() = default_value;
            syncTree->jptraw_2() = default_value;
    //            syncTree->jptunc_2() = default_value;
            syncTree->jmva_2() = default_value;
            syncTree->jpass_2() = default_value;

            if (event.njets >= 1) {
                syncTree->jpt_1() = event.pt_jets.at(0);
                syncTree->jeta_1() = event.eta_jets.at(0);
                syncTree->jphi_1() = event.phi_jets.at(0);
                syncTree->jptraw_1() = event.ptraw_jets.at(0);
                //syncTree->jptunc_1();
                syncTree->jmva_1() = event.mva_jets.at(0);
                syncTree->jpass_1() = event.passPU_jets.at(0);
            }

            if (event.njets >= 2) {
                syncTree->jpt_2() = event.pt_jets.at(1);
                syncTree->jeta_2() = event.eta_jets.at(1);
                syncTree->jphi_2() = event.phi_jets.at(1);
                syncTree->jptraw_2() = event.ptraw_jets.at(1);
                //syncTree->jptunc_2();
                syncTree->jmva_2() = event.mva_jets.at(1);
                syncTree->jpass_2() = event.passPU_jets.at(1);
            }

            syncTree->nbtag() = event.nBjets_retagged;
            //syncTree->nbtag() = event.nBjets;

            syncTree->bpt_1() = default_value;
            syncTree->beta_1() = default_value;
            syncTree->bphi_1() = default_value;
            syncTree->bcsv_1() = default_value;
            syncTree->bpt_2() = default_value;
            syncTree->beta_2() = default_value;
            syncTree->bphi_2() = default_value;
            syncTree->bcsv_2() = default_value;
            syncTree->m_bb() = default_value;
            syncTree->m_ttbb() = default_value;
            syncTree->bpt_3() = default_value;
            syncTree->beta_3() = default_value;
            syncTree->bphi_3() = default_value;
            syncTree->bcsv_3() = default_value;

            if (event.nBjets >= 1 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM) {
                syncTree->bpt_1() = event.pt_Bjets.at(0);
                syncTree->beta_1() = event.eta_Bjets.at(0);
                syncTree->bphi_1() = event.phi_Bjets.at(0);
                syncTree->bcsv_1() = event.csv_Bjets.at(0);
            }

            if (event.nBjets >= 2 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(1) > cuts::Htautau_Summer13::btag::CSVM) {
                syncTree->bpt_2() = event.pt_Bjets.at(1);
                syncTree->beta_2() = event.eta_Bjets.at(1);
                syncTree->bphi_2() = event.phi_Bjets.at(1);
                syncTree->bcsv_2() = event.csv_Bjets.at(1);
                syncTree->m_bb() = eventInfo.Hbb.M();
                syncTree->m_ttbb() = eventInfo.resonance.M();
            }

            if (event.nBjets >= 3 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(1) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(2) > cuts::Htautau_Summer13::btag::CSVM){
                syncTree->bpt_3() = event.pt_Bjets.at(2);
                syncTree->beta_3() = event.eta_Bjets.at(2);
                syncTree->bphi_3() = event.phi_Bjets.at(2);
                syncTree->bcsv_3() = event.csv_Bjets.at(2);
            }

            syncTree->Fill();
        }
    }

private:
    std::shared_ptr<ntuple::SyncTree> syncTree;
    EventId_ToES_Map eventId_ToES_Map;
    unsigned inclusive, passed;
};
