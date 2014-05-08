/*!
 * \file custom_cuts.h
 * \brief Cuts which are customized for HHbbtautau analysis
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-01 created
 */

#pragma once

namespace cuts {
    const double minDeltaR_betweenHiggses = 1; // >

    namespace btag {
        const double CSVL = 0.244; // loose twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
        const double CSVM = 0.679; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
        const double CSVT = 0.898; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5

        namespace signal {
            const double pt = 20; // > AN-2013/075 HH->bb gamma gamma
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

}
