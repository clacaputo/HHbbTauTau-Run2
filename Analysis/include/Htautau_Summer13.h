/*!
 * \file Htautau_Summer13.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-13 created
 */

#pragma once
#include <string>
#include <vector>

namespace cuts {
namespace Htautau_Summer13 {

const double DeltaR_betweenSignalObjects = 0.5; // >

namespace tauID {
    namespace ETau {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-10-82
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronMediumMVA3 = 0.5; //  > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                      // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                      // MVA3 is recommended, but it does not exists any more
        const double LooseCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                      // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
    namespace MuTau {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-10-82
        const double againstMuonTight = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double LooseCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                      // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
    namespace TauTau {
        //const double pt = 45; // > AN-13-189
        const double pt_lead = 25; // > AN-2013-187
        const double pt_sublead = 20; // > AN-2013-187
        const double eta = 2.1; // < AN-13-189
        const double decayModeFinding = 0.5; // > AN-10-82
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA3 = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more
        const double MediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }

    namespace veto {
        const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double eta = 2.5; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double decayModeFinding = 0.5;  // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double LooseCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                      // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double dz = 0.2; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double deltaR_signalObjects = 0.4; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
    }
}

namespace muonID {
    namespace MuTau {
        const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double eta = 2.1; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state

        const bool isTightMuon = true; // = HiggsToTauTauWorkingSummer2013#Muon_ID

        //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
        const bool isGlobalMuonPromptTight = true;
              // = https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc#567 and 590 definition
              // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
        const bool isPFMuon = true; // = def of isTightMuon
        const int nMatched_Stations = 1; // > def of isTightMuon
        const double dB = 0.045; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID
        const double dz = 0.2; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID
        const int pixHits = 0; // > def of isTightMuon
        const int trackerLayersWithMeasurement = 5; // > def of isTightMuon

        const double pFRelIso = 0.1; //  < twiki SWGuideMuonId#Muon_Isolation_AN1
    }

    namespace veto {
        const double pt = 5; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const bool isTightMuon = true; // = HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double pFRelIso = 0.15; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double dz = 0.2; // < cm, see HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double deltaR_signalObjects = 0.4; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
    }
}

namespace electronID{
    const double eta_CrackVeto_low = 1.4442;
    const double eta_CrackVeto_high = 1.566;

    namespace ETau {
        const double pt = 24; // >  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double eta_high = 2.1; // <  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double dz = 0.2; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const int missingHits = 1; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const bool hasMatchedConversion = false; // =  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double dB = 0.045; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double pFRelIso = 0.1; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Muon_Isolation
        const double scEta_min[2] = { 0.8, 1.479 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID
        const double MVApogNonTrig[3] = { 0.925, 0.975, 0.985 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID
    }

    namespace veto {
        const double pt = 10; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double eta_high = 2.5; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double dz = 0.2; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double pFRelIso = 0.3; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double ref_pt = 20; // twiki HiggsToTauTauWorkingSummer2013#Electron_ID
        const double scEta_min[2] = {0.8, 1.479}; // loose HiggsToTauTauWorkingSummer2013#Electron_ID
        const double MVApogNonTrig[2][3] = {{0.925, 0.915, 0.965},{0.905,0.955, 0.975}};
                                                  // loose HiggsToTauTauWorkingSummer2013#Electron_ID
        const double deltaR_signalObjects = 0.4; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
    }
}

namespace btag {
    const double CSVL = 0.244; // loose twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
    const double CSVM = 0.679; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
    const double CSVT = 0.898; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5

    namespace signal {
        const double pt = 20; // > AN-13-075 HH bb gamma gamma
        const double eta = 2.4; // < twiki HiggsToTauTauWorkingSummer2013#bjets
        const double CSV = CSVM; // recommended twiki HiggsToTauTauWorkingSummer2013#bjets
    }
    namespace veto {
        const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double eta = 2.4; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double CSV = CSVT; // twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double passLooseID = true; // = twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        const double deltaR_signalObjects = 0.4; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
    }
}

namespace vertex {
    const double ndf = 4; // >= twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_mu_mu_tau
    const double z = 24.0; // < cm, twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_mu_mu_tau
    const double r = 2.0; // < cm, twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_mu_mu_tau
    const bool chooseHighestSumPt2 = true; // = AN-13-188
}

namespace trigger {
    namespace MuTau {
        const std::vector<std::string> hltPaths =
                {"HLT_IsoMu17_eta2p1_LooseIsoPFTau20","HLT_IsoMu18_eta2p1_LooseIsoPFTau20"};

    }
    namespace ETau {
        const std::vector<std::string> hltPaths =
            {"HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20",
             "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20"};
    }

    namespace TauTau {
        const std::vector<std::string> hltPaths =
                {"HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1",
                 "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1" };
    }
}

} // Htautau_Summer13
} // cuts
