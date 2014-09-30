/*!
 * \file Htautau_Summer13.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-02-13 created
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
#include <string>
#include <vector>
#include <map>

#include <TMath.h>

#include "AnalysisBase/include/Tools.h"
#include "TreeProduction/interface/Tau.h"
#include "TreeProduction/interface/Jet.h"

namespace cuts {
namespace Htautau_Summer13 {
// AN-2013/178 H->etau,mutau
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/analyzers/DiLeptonAnalyzer.py
const double DeltaR_betweenSignalObjects = 0.5; // >

// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauMuAnalyzer.py#L272
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L665
const double DeltaR_triggerMatch = 0.5; // <

namespace MuTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const std::vector<std::string> hltPaths =
            { "HLT_IsoMu17_eta2p1_LooseIsoPFTau20", "HLT_IsoMu18_eta2p1_LooseIsoPFTau20" };
    }

    // twiki HiggsToTauTauWorkingSummer2013
    namespace muonIDscaleFactor {
        const std::vector<double> eta = { 0.8, 1.2, 2.1 };
        const std::vector<double> pt = { 20, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.9818, 0.9829, 0.9869 },
                                                                    { 0.9852, 0.9852, 0.9884 } };
    }

    // twiki HiggsToTauTauWorkingSummer2013
    namespace muonISOscaleFactor {
        const std::vector<double> eta = { 0.8, 1.2, 2.1 };
        const std::vector<double> pt = { 20, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.9494, 0.9835, 0.9923 },
                                                                    { 0.9883, 0.9937, 0.9996 } };
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

        const double mt = 30; // < not used for sync, only for the final selection.
    }

    namespace tauID {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonTight = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation

        // https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsHTohh/test/HiggsHTohh.cpp#L1086
        const double dz = 0.2; // <
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

    namespace BackgroundEstimation {
        const double HighMtRegion = 70; // > For W-jets data driven estimation
    }
}

namespace ETau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const std::vector<std::string> hltPaths =
            { "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20",
              "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20" };
    }

    // twiki HiggsToTauTauWorkingSummer2013
    namespace electronIDscaleFactor {
        const std::vector<double> eta = { 1.479, 2.1 };
        const std::vector<double> pt = { 24, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.8999, 0.7945 },
                                                                    { 0.9486, 0.8866 } };
    }

    // twiki HiggsToTauTauWorkingSummer2013
    namespace electronISOscaleFactor {
        const std::vector<double> eta = { 1.479, 2.1 };
        const std::vector<double> pt = { 24, 30 };
        const std::vector< std::vector< double > > scaleFactors = { { 0.9417, 0.9471 },
                                                                    { 0.9804, 0.9900 } };
    }

    namespace electronID{
        const double pt = 24; // >  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double eta_high = 2.1; // <  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double dz = 0.2; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const int missingHits = 1; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const bool hasMatchedConversion = false; // =  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double d0 = 0.045; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double pFRelIso = 0.1; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Muon_Isolation
        const double scEta_min[2] = { 0.8, 1.479 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID
        const double MVApogNonTrig[3] = { 0.925, 0.975, 0.985 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID

        const double mt = 30; // < not used for sync, only for the final selection.
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

        // > https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/physicsobjects/Tau.py#L62
        const std::vector<double> againstElectronMediumMVA3_customValues = {
            0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95,0.897,0.958,0.955,0.942
        };
        // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauEleAnalyzer.py#L187
        // https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsHTohh/test/HiggsHTohh.cpp#L1060 (no dB)
        const double dz = 0.2; // <
        const double dB = 0.045; // <
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

    namespace BackgroundEstimation {
        const double HighMtRegion = 70; // > For W-jets data driven estimation
    }
}

namespace TauTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Tau_Tau_Final_state
        const std::map< std::string, bool > hltPathsMap =
        {{"HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1", false},
         {"HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1",false} };

        const std::vector<std::string> hltPaths = tools::collect_map_keys(hltPathsMap);
    }

    namespace tauID {
        const double pt = 45; // > AN-2013/189 H->tautau full hadronic
        const double eta = 2.1; // < AN-2013/189 H->tautau full hadronic
        const double dz = 0.2; // < AN-2013/189 H->tautau full hadronic
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA3 = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading pt tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more for new tauID
        const double byMediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1; // not equivalent to the medium working point (1.5)
                                                     // custom value used for QCD estimation
                                                     // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation

        namespace BackgroundEstimation {
            const double Isolation_upperLimit = 4; // < upper value of isolation for QCD data driven estimation
        }

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
    const int missingHits = 1; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
    const bool hasMatchedConversion = false; // =  HiggsToTauTauWorkingSummer2013#Electron_ID
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
    const bool puLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >

    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const bool pfLooseID = true; // =

    const double pt_loose = 20; // >

    //https://github.com/ajgilbert/ICHiggsTauTau/blob/38f0bedebe1e7ff432bdcbd7753f38cfaf95405f/plugins/MVAMETPairProducer.cc#L410
    inline bool passPFLooseId(const ntuple::Jet& jet)
    {
        TLorentzVector momentum;
        momentum.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass);
        if(momentum.E() == 0)                                  return false;
        if(jet.neutralHadronEnergyFraction > 0.99)   return false;
        if(jet.neutralEmEnergyFraction     > 0.99)   return false;
        if(jet.nConstituents <  2)                          return false;
        if(jet.chargedHadronEnergyFraction <= 0 && std::abs(jet.eta) < 2.4 ) return false;
        if(jet.chargedEmEnergyFraction >  0.99  && std::abs(jet.eta) < 2.4 ) return false;
        if(jet.chargedMultiplicity     < 1      && std::abs(jet.eta) < 2.4 ) return false;
        return true;
    }
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
    const bool puLooseID = true; // =
    const bool pfLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >
}

// AN-2013/188 H->tautau physics objects
namespace vertex {
    const double ndf = 4; // >
    const double z = 24.0; // < cm
    const double r = 2.0; // < cm
    const bool chooseHighestSumPt2 = true; // =
}

namespace tauCorrections {

    const double DecayModeWeight = 0.88; // = HiggsToTauTauWorkingSummer2013#TauTau_scale_factors
                                         // for 1-prong no pi 0 taus

    const double deltaR = 0.3; // < Updated to be compatible with H->tautau code

    const double energyUncertainty = 0.03;

    // For taus that matched MC truth.
    // Original corrections from HiggsToTauTauWorkingSummer2013. Updated to be compatible with H->tautau code.
    inline double MomentumScaleFactor(bool hasMCmatch, double pt, ntuple::tau_id::hadronicDecayMode decayMode,
                                      bool useLegacyCorrections)
    {
        if(!hasMCmatch) return 1.0;
        if(decayMode == ntuple::tau_id::kOneProng1PiZero || decayMode == ntuple::tau_id::kOneProng2PiZero) {
            if(useLegacyCorrections)
                return 1.025 + 0.001 * std::min(std::max(pt - 45.0, 0.0), 10.0);
            return 1.012;
        }
        if(decayMode == ntuple::tau_id::kThreeProng0PiZero) {
            if(useLegacyCorrections)
                return 1.012 + 0.001 * std::min(std::max(pt - 32.0, 0.0), 18.0);
            return 1.012;
        }
        return 1.0;
    }
}

namespace customTauMVA {
    bool ComputeAntiElectronMVA3New(int category, float raw, int WP)
    {
      if (category < 0 ) return false ;
      if (category > 15) return true  ;

      float cutsLoose    [16] = {0.835,0.831,0.849,0.859,0.873,0.823,0.85 ,0.855,0.816,0.861,0.862,0.847,0.893,0.82 ,0.845,0.851} ;
      float cutsMedium   [16] = {0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95 ,0.897,0.958,0.955,0.942} ;
      float cutsTight    [16] = { 0.96,0.968,0.971,0.972,0.969,0.959,0.981,0.965,0.975,0.972,0.974,0.971,0.897,0.971,0.961,0.97 } ;
      float cutsVeryTight[16] = {0.978,0.98 ,0.982,0.985,0.977,0.974,0.989,0.977,0.986,0.983,0.984,0.983,0.971,0.987,0.977,0.981} ;

      switch (WP)
      {
        case 0 : return (raw > cutsLoose    [category]) ;
        case 1 : return (raw > cutsMedium   [category]) ;
        case 2 : return (raw > cutsTight    [category]) ;
        case 3 : return (raw > cutsVeryTight[category]) ;
      }

      return false ; // to avoid warnings, smarter solutions exist
    }

}

} // Htautau_Summer13
} // cuts
