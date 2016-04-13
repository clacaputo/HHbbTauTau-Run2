/*!
 * \file Htautau_2015.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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
#include "../../AnalysisBase/include/AnalysisMath.h"
#include "../../AnalysisBase/include/Tools.h"
#include "../../AnalysisBase/include/exception.h"
#include "../../TreeProduction/interface/Tau.h"
#include "../../TreeProduction/interface/Jet.h"

namespace analysis{

enum class DataSourceType { Spring15MC=0, Fall15MC=1, Run2015B=2, Run2015C , Run2015D };

static const std::map<DataSourceType,std::string> dataSourceTypeMap = {{DataSourceType::Spring15MC , "Spring15MC"},
                                                                       {DataSourceType::Fall15MC , "Fall15MC"},
                                                                       {DataSourceType::Run2015B , "Run2015B"},
                                                                       {DataSourceType::Run2015C , "Run2015C"},
                                                                       {DataSourceType::Run2015D , "Run2015D"}};

static const std::map<std::string, DataSourceType> stringToDataSourceTypeMap = {{ "Spring15MC", DataSourceType::Spring15MC},
                                                                       { "Fall15MC", DataSourceType::Fall15MC},
                                                                       { "Run2015B" , DataSourceType::Run2015B},
                                                                       { "Run2015C" , DataSourceType::Run2015C},
                                                                       { "Run2015D" , DataSourceType::Run2015D}};


inline std::ostream& operator<< (std::ostream& s, const DataSourceType& dataSourceType) {
    s << dataSourceTypeMap.at(dataSourceType);
    return s;
}

}//namespace analysis

namespace cuts {
namespace Htautau_2015 {
// AN-2013/178 H->etau,mutau
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/analyzers/DiLeptonAnalyzer.py
const double DeltaR_betweenSignalObjects = 0.5; // >

// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauMuAnalyzer.py#L272
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L665
const double DeltaR_triggerMatch = 0.5; // <

namespace MuTau {
    namespace trigger {
        // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015  - Trigger Session
        const std::map<analysis::DataSourceType , const std::set<std::string>> hltPathMaps =
                                    {{analysis::DataSourceType::Spring15MC,{"HLT_IsoMu17_eta2p1_v1"}},
                                     {analysis::DataSourceType::Fall15MC,{"HLT_IsoMu18_v2"}},
                                     {analysis::DataSourceType::Run2015B,{"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",
                                                                            "HLT_IsoMu24_eta2p1_v2"}},
                                     {analysis::DataSourceType::Run2015C,{"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",
                                                                          "HLT_IsoMu24_eta2p1_v2"}},
                                     {analysis::DataSourceType::Run2015D,{"HLT_IsoMu18_v"}}};
    }

    namespace muonID {
        const double pt = 19; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.1; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dB = 0.045; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        const bool isMediumMuon = true; // = twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
                                        // def : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Medium_Muon

        //After Synctuple

        const double pFRelIso = 0.1; // < twiki HiggsToTauTauWorking2015#Muons

     /*   //const bool isTightMuon = true; // = HiggsToTauTauWorkingSummer2013#Muon_ID
                                         //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
        const bool isGlobalMuonPromptTight = true;
              // = https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/DataFormats/MuonReco/src/MuonSelectors.cc#L534
              // def: https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/DataFormats/MuonReco/src/MuonSelectors.cc#L534
              // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
        const bool isPFMuon = true; // = def of isTightMuon
        const int nMatched_Stations = 1; // > def of isTightMuon
        const int pixHits = 0; // > def of isTightMuon
        const int trackerLayersWithMeasurement = 5; // > def of isTightMuon

        const double mt = 30; // < not used for sync, only for the final selection.
        */
    }

    namespace tauID {
        const double pt = 20; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.3; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double decayModeFinding = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        //After Synctuple

        const double againstMuonTight3 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronVLooseMVA5 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

   }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
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
//        const std::set<std::string> hltPaths =
//            { "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20",
//              "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20" };

        const std::map<analysis::DataSourceType , const std::set<std::string>> hltPathMaps =
                                    {{analysis::DataSourceType::Spring15MC,{"HLT_Ele22_eta2p1_WP75_Gsf"}},
                                     {analysis::DataSourceType::Fall15MC,{"HLT_Ele23_WPLoose_Gsf_v3"}},
                                     {analysis::DataSourceType::Run2015C,{"HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1",
                                                                          "HLT_Ele32_eta2p1_WPTight_Gsf_v1"}},
                                     {analysis::DataSourceType::Run2015D,{"HLT_Ele23_WPLoose_Gsf_v"}}};
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

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#Baseline_e_tau_h_AN1
    namespace electronID{
        const double pt = 24; // >
        const double eta = 2.1; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // <
        const bool MVApogNonTrig = true; // tight Use Value map compute by POG people
        const int missingHits = 1; // <
        const bool passConversionVeto = true; // =

        const double pFRelIso = 0.1; // < twiki
        const double scEta_min[2] = { 0.8, 1.479 }; // tight


        const double mt = 30; // < not used for sync, only for the final selection.
    }

    namespace tauID {
        const double pt = 20; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double eta = 2.3; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double decayModeFinding = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double dz = 0.2; // < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

        //After Synctuple

        const double againstMuonTight3 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double againstElectronVLooseMVA5 = 0.5; // > twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorking2015#Baseline_mu_tau_h

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
        const std::map< std::string, bool > hltPathsMap =
        {{"HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1", false},
         {"HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1",false} };

        const std::set<std::string> hltPaths = analysis::tools::collect_map_keys(hltPathsMap);

        const double jet_pt = 50; // > GeV twiki HiggsToTauTauWorkingSummer2013#Tau_Tau_Final_state
        const double jet_eta = 3.0; // < twiki HiggsToTauTauWorkingSummer2013#Tau_Tau_Final_state
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
                                                     // MVA3 is recommended, but it does not exists any more
                                                     // for new tauID
        const double byMediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1; // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                //https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/RecoTauTag/Configuration/python/HPSPFTaus_cff.py#L183

        namespace BackgroundEstimation {
            const double Isolation_upperLimit = 4; // < upper value of isolation for QCD data driven estimation
        }

    }
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorking2015#Third_lepton_vetoes
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
    const double dB = 0.045; // <

    //def of isMediumMuon: twiki SWGuideMuonId#Medium_Muon
    const bool isMediumMuon = true; // =

    const double pfRelIso = 0.3; // <
}

namespace jetID {
    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 30; // > , njets is filled with a 30 GeV cuts. All the other variables use jets with 20 pt cut
    const double eta = 4.7; // <
    const bool puLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >

    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const bool pfLooseID = true; // =

    const double pt_loose = 20; // >

    //  https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    inline bool passPFLooseId(const pat::Jet& jet)
    {
        TLorentzVector momentum;
        momentum.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
        if(momentum.E() == 0)                                  return false;
        if(jet.neutralHadronEnergyFraction() > 0.99)   return false;
        if(jet.neutralEmEnergyFraction()     > 0.99)   return false;
        if(jet.nConstituents() <  1)                   return false;
        if(jet.chargedHadronEnergyFraction() <= 0 && std::abs(jet.eta()) < 2.4 ) return false;
        if(jet.chargedEmEnergyFraction() >  0.99  && std::abs(jet.eta()) < 2.4 ) return false;
        if(jet.chargedMultiplicity()     <= 0      && std::abs(jet.eta()) < 2.4 ) return false;
        return true;
    }
}

namespace btag {
    // pfCombinedInclusiveSecondaryVertexV2BJetTags, twiki
    const double CSVL = 0.605; // > loose
    const double CSVM = 0.89; // > medium
    const double CSVT = 0.97; // > tight

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 20; // >
    const double eta = 2.4; // <
    const double CSV = 0.8; // > HTauTau Twiki
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

    const double deltaR_matchGenParticle = 0.5; // gen Particle match
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

namespace jetToTauFakeRateWeight {

    inline double CalculateJetToTauFakeWeight(double tauPt)
    {
        const double tau_pt = tauPt < 200 ? tauPt : 200 ;
        return (1.15743)-(0.00736136*tau_pt)+(4.3699e-05*tau_pt*tau_pt)-(1.188e-07*tau_pt*tau_pt*tau_pt);
    }
}

// twiki HiggsToTauTauWorkingSummer2013
namespace electronEtoTauFakeRateWeight {

    //inline double CalculateEtoTauFakeWeight(const ntuple::Tau& tau_leg)
    inline double CalculateEtoTauFakeWeight(double tau_eta, ntuple::tau_id::hadronicDecayMode tau_decayMode)
    {
        static const double eta = 1.5 ;
        static const std::map<ntuple::tau_id::hadronicDecayMode, std::vector<double>> decayModeMap= {
        { ntuple::tau_id::kOneProng0PiZero, {1.37 , 1.11} }, { ntuple::tau_id::kOneProng1PiZero, {2.18 , 0.47} } };

        if (!decayModeMap.count(ntuple::tau_id::ConvertToHadronicDecayMode(tau_decayMode)))
            return 1;
        const size_t eta_bin = std::abs(tau_eta) < eta ? 0 : 1;
        return decayModeMap.at(ntuple::tau_id::ConvertToHadronicDecayMode(tau_decayMode)).at(eta_bin);
    }

}

namespace DrellYannCategorization {
    const double minimal_genParticle_pt = 8; // > GeV
    const double deltaR_matchGenParticle = 0.5; // < twiki HiggsToTauTauWorkingSummer2013#E_MU_TAU_channel
    // < 0.3 https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L757

    const double minimal_visible_momentum = 18; // > GeV
}

namespace DYEmbedded {
    namespace trigger {
        // Riccardo repository: https://github.com/rmanzoni/HTT/tree/
        //a03b227073b2d4d8a2abe95367c014694588bf98/CMGTools/H2TauTau/python/proto/samples/run2012/trigger*.py
        const std::set<std::string> hltPaths = { "HLT_Mu17_Mu8"};
    }

    const double invariantMassCut = 50.; // > GeV
           // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/EmbedWeighter.py#L132

    const double deltaR_matchGenParticle = DrellYannCategorization::deltaR_matchGenParticle;
}

namespace customTauMVA {
    // https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/FnPredicates.cc#L319
    inline bool ComputeAntiElectronMVA3New(const ntuple::Tau& tau, size_t WP, bool applyDzCut = false)
    {
        static size_t n_categories = 16;
        static const std::vector< std::vector<float> > cuts {
            { 0.835, 0.831, 0.849, 0.859, 0.873, 0.823, 0.85, 0.855, 0.816, 0.861, 0.862, 0.847, 0.893, 0.82,
                0.845, 0.851 }, // loose
            { 0.933, 0.921, 0.944, 0.945, 0.918, 0.941, 0.981, 0.943, 0.956, 0.947, 0.951, 0.95, 0.897, 0.958,
                0.955, 0.942 }, // medium
            { 0.96, 0.968, 0.971, 0.972, 0.969, 0.959, 0.981, 0.965, 0.975, 0.972, 0.974, 0.971, 0.897, 0.971,
                0.961, 0.97 }, // tight
            { 0.978, 0.98, 0.982, 0.985, 0.977, 0.974, 0.989, 0.977, 0.986, 0.983, 0.984, 0.983, 0.971, 0.987,
                0.977, 0.981 } // very tight
        };

        if(applyDzCut) {
            const TLorentzVector momentum = analysis::MakeLorentzVectorPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
            const double z_2 = tau.vz + (130. / std::tan(momentum.Theta()));
            if (z_2 > -1.5 && z_2 < 0.5) return false;
        }

        const int category = std::round(tau.againstElectronMVA3category);
        const float raw = tau.againstElectronMVA3raw;

        if(category < 0) return false;
        const size_t u_category = static_cast<size_t>(category);
        if(u_category >= n_categories) return true;

        if(WP >= cuts.size())
            throw std::runtime_error("Bad working point");

        return raw > cuts.at(WP).at(u_category);
    }
}

} // Htautau_2015
} // cuts
