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
// AN-2013/178 H->etau,mutau
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/analyzers/DiLeptonAnalyzer.py
const double DeltaR_betweenSignalObjects = 0.5; // >

namespace MuTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const std::vector<std::string> hltPaths =
                {"HLT_IsoMu17_eta2p1_LooseIsoPFTau20","HLT_IsoMu18_eta2p1_LooseIsoPFTau20"};
    }

    namespace muonID {
        const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double eta = 2.1; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state

//        const bool isTightMuon = true; // = HiggsToTauTauWorkingSummer2013#Muon_ID
        //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
        const bool isGlobalMuonPromptTight = true;
              // = https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc#567 and 590 definition
              // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
        const bool isPFMuon = true; // = def of isTightMuon
        const int nMatched_Stations = 1; // > def of isTightMuon
        const int pixHits = 0; // > def of isTightMuon
        const int trackerLayersWithMeasurement = 5; // > def of isTightMuon
        const double dz = 0.2; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID
        const double dB = 0.045; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID

        const double pFRelIso = 0.1; // < twiki SWGuideMuonId#Muon_Isolation_AN1
    }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
    namespace ZmumuVeto {
        const double pt = 15; // >
        const double eta = 2.4; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // < same definition as db
        const bool isGlobalMuon = true; // = already applied at Tree selection level as skim - MuonBlock.cc
        const bool isTrackerMuon = true; // =
        const bool isPFMuon = true; // =
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =
    }
}

namespace ETau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const std::vector<std::string> hltPaths =
            {"HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20",
             "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20"};
    }

    namespace electronID{
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

    namespace tauID {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronMediumMVA3 = 0.5; //  > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                      // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                      // MVA3 is recommended, but it does not exists any more
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
    namespace ZeeVeto {
        const double pt = 15; // >
        const double eta = 2.5; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // <
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =

        // WP95 definition for Veto twiki EgammaCutBasedIdentification
        const double barrel_eta_high = 1.479; // <=
        const double endcap_eta_low = 1.479; // >
        const double endcap_eta_high = 2.5; // <
        const size_t barrel_index = 0, endcap_index = 1;
        const double sigma_ieta_ieta[] = { 0.01, 0.03 }; // <
        const double delta_eta[] = { 0.007, 0.01 }; // <
        const double delta_phi[] = { 0.8, 0.7 }; // <
        const double HoverE[] = { 0.15, std::numeric_limits<double>::max() }; // <
        const double dZ_vtx[] = { 0.2, 0.2 }; // <
    }
}

namespace TauTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Tau_Tau_Final_state
        const std::vector<std::string> hltPaths =
                {"HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30",
                 "HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1",
                 "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1" };
    }

    namespace tauID {
        const double pt = 45; // > AN-2013/189 H->tautau full hadronic
        const double eta = 2.1; // < AN-2013/189 H->tautau full hadronic
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA3 = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading pt tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more for new tauID
        const double byMediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/physicsobjects/HTauTauElectron.py
namespace electronVeto {
    const double pt = 10; // >
    const double eta_high = 2.5; // <
    const double dz = 0.2; // <
    const double d0 = 0.045; // <
    const double pFRelIso = 0.3; // <
    const double ref_pt = 20; // twiki HiggsToTauTauWorkingSummer2013#Electron_ID
    const double scEta_min[2] = {0.8, 1.479}; // loose HiggsToTauTauWorkingSummer2013#Electron_ID
    const double MVApogNonTrig[2][3] = {{0.925, 0.915, 0.965},{0.905,0.955, 0.975}};
                                              // loose HiggsToTauTauWorkingSummer2013#Electron_ID
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py
namespace muonVeto {
    const double pt = 10; // >
    const double eta = 2.4; // <
    const double dz = 0.2; // <
    const double d0 = 0.045; // <

    //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
    const bool isGlobalMuonPromptTight = true;
          // = https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc#567 and 590 definition
          // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
    const bool isPFMuon = true; // = def of isTightMuon
    const int nMatched_Stations = 1; // > def of isTightMuon
    const int pixHits = 0; // > def of isTightMuon
    const int trackerLayersWithMeasurement = 5; // > def of isTightMuon

    const double pfRelIso = 0.3; // <
}

namespace jetID {
    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 30; // >
    const double eta = 4.7; // <
    const double puLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >

    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const double pfLooseID = true; // =
}

namespace btag {
    // twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
    const double CSVL = 0.244; // > loose
    const double CSVM = 0.679; // > medium
    const double CSVT = 0.898; // > tight

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 20; // >
    const double eta = 2.4; // <
    const double CSV = CSVM; // >
    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const double puLooseID = true; // =
    const double pfLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >
}

// AN-2013/188 H->tautau physics objects
namespace vertex {
    const double ndf = 4; // >
    const double z = 24.0; // < cm
    const double r = 2.0; // < cm
    const bool chooseHighestSumPt2 = true; // =
}

} // Htautau_Summer13
} // cuts
