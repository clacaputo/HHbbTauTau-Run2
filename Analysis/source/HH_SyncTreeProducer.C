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

class HH_SyncTreeProducerData : public analysis::LightFlatAnalyzerData {
public:
    HH_SyncTreeProducerData(TFile& outputFile) : LightFlatAnalyzerData(outputFile) {}


};

class HH_SyncTreeProducer : public analysis::LightBaseFlatTreeAnalyzer {
public:
    HH_SyncTreeProducer(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), anaData(GetOutputFile())
    {
        anaData.getOutputFile().cd();
        syncTree = std::shared_ptr<ntuple::SyncTree>(new ntuple::SyncTree("syncTree"));
    }

    virtual ~HH_SyncTreeProducer() {
        syncTree->Write();
    }
protected:
    virtual analysis::LightFlatAnalyzerData& GetAnaData() override { return anaData; }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;

        const analysis::EventRegionSet interestedEventRegion = {analysis::EventRegion::OS_Isolated,
                                                     analysis::EventRegion::SS_Isolated};

        if (category != EventCategory::Inclusive) return;
        const analysis::EventRegion eventRegion = DetermineEventRegion(eventInfo);
        if (!interestedEventRegion.count(eventRegion)) return;

        const ntuple::Flat& event = *eventInfo.event;
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
        syncTree->m_sv() = event.m_sv_vegas;
        syncTree->pt_sv() = default_value;
        syncTree->eta_sv() = default_value;
        syncTree->phi_sv() = default_value;
        syncTree->m_sv_Up() = event.m_sv_up_vegas;
        syncTree->m_sv_Down() = event.m_sv_down_vegas;

        syncTree->pt_1() = event.pt_1;
        syncTree->phi_1() = event.phi_1;
        syncTree->eta_1() = event.eta_1;
        syncTree->m_1() = event.m_1;
        syncTree->q_1() = event.q_1;
        syncTree->mt_1() = event.mt_1;
        syncTree->d0_1() = event.d0_1;
        syncTree->dZ_1() = event.dZ_1;

        // leg1 lepton specific variable should be filled outside. Here all them set to the default value.
        syncTree->iso_1() = event.pfRelIso_1;
        syncTree->mva_1() = event.mva_1;
        syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
        syncTree->againstElectronMVA3raw_1() = default_value;
        syncTree->byIsolationMVA2raw_1() = default_value;
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

        syncTree->iso_2() = default_value;
        syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
        syncTree->againstElectronMVA3raw_2() = default_value;
        syncTree->byIsolationMVA2raw_2() = default_value;
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

//        if (event.njets >= 1) {
//            const Candidate& jet = jets.at(0);
//            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);
//            syncTree->jpt_1() = event
//            syncTree->jeta_1() = jet.momentum.Eta();
//            syncTree->jphi_1() = jet.momentum.Phi();
//            syncTree->jptraw_1() = ntuple_jet.pt_raw;
//            //syncTree->jptunc_1();
//            syncTree->jmva_1() = ntuple_jet.puIdMVA;
//            //syncTree->jlrm_1();
//            //syncTree->jctm_1();
//            syncTree->jpass_1() = ntuple::JetID_MVA::PassLooseId(ntuple_jet.puIdBits);
//        } else {
//            syncTree->jpt_1() = default_value;
//            syncTree->jeta_1() = default_value;
//            syncTree->jphi_1() = default_value;
//            syncTree->jptraw_1() = default_value;
//            //syncTree->jptunc_1() = default_value;
//            syncTree->jmva_1() = default_value;
////            syncTree->jlrm_1() = default_value;
////            syncTree->jctm_1() = default_value;
//            syncTree->jpass_1() = default_value;
//        }

//        if (jets.size() >= 2) {
//            const Candidate& jet = jets.at(1);
//            const ntuple::Jet& ntuple_jet = event->jets().at(jet.index);
//            syncTree->jpt_2() = jet.momentum.Pt();
//            syncTree->jeta_2() = jet.momentum.Eta();
//            syncTree->jphi_2() = jet.momentum.Phi();
//            syncTree->jptraw_2() = ntuple_jet.pt_raw;
//            //syncTree->jptunc_2();
//            syncTree->jmva_2() = ntuple_jet.puIdMVA;
//            //syncTree->jlrm_2();
//            //syncTree->jctm_2();
//            syncTree->jpass_2() = ntuple::JetID_MVA::PassLooseId(ntuple_jet.puIdBits);
//        } else {
//            syncTree->jpt_2() = default_value;
//            syncTree->jeta_2() = default_value;
//            syncTree->jphi_2() = default_value;
//            syncTree->jptraw_2() = default_value;
////            syncTree->jptunc_2() = default_value;
//            syncTree->jmva_2() = default_value;
////            syncTree->jlrm_2() = default_value;
////            syncTree->jctm_2() = default_value;
//            syncTree->jpass_2() = default_value;
//        }

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

private:
    HH_SyncTreeProducerData anaData;
    std::shared_ptr<ntuple::SyncTree> syncTree;
};

