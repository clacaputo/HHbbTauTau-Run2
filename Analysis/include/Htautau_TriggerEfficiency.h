/*!
 * \file Htautau_TriggerEfficiency.h
 * \brief Definition of functions to calculate trigger efficiencies.Higgs.
 *
 * Code copied from CMG Tools:
 *   https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/interface/TriggerEfficiency.h
 *   https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/src/TriggerEfficiency.cc
 *
 * For description see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 */

#pragma once

#include <vector>
#include <TMath.h>
#include <TLorentzVector.h>

namespace analysis {
namespace Htautau_Summer13 {

namespace trigger {
namespace detail {
inline double efficiency(double m, double m0, double sigma, double alpha, double n, double norm)
{
    if(m<1. || 1000.<m) return 0.; //safety check

    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) sigma);
    double t = (m - m0)/sig;
    if(alpha < 0) t = -t;
    double absAlpha = fabs(alpha/sig);
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = absAlpha - n/absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
    double area = leftArea + rightArea;
    if( t <= absAlpha ) {
        arg = t / sqrt2;
        if(arg > 5.) ApproxErf = 1;
        else if (arg < -5.) ApproxErf = -1;
        else ApproxErf = TMath::Erf(arg);
        return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
    }
    else {
        return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -  1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
    }
}

template<typename EffFunction>
std::vector<double> CalculateWeights(const TLorentzVector& momentum_leg1, const TLorentzVector& momentum_leg2,
                                     const EffFunction& eff_leg1_data_fn, const EffFunction& eff_leg2_data_fn,
                                     const EffFunction& eff_leg1_mc_fn, const EffFunction& eff_leg2_mc_fn)
{
    std::vector<double> weights;
    const double eff_leg1_data = eff_leg1_data_fn(momentum_leg1.Pt(), momentum_leg1.Eta());
    const double eff_leg2_data = eff_leg2_data_fn(momentum_leg2.Pt(), momentum_leg2.Eta());
    const double eff_leg1_mc = eff_leg1_mc_fn(momentum_leg1.Pt(), momentum_leg1.Eta());
    const double eff_leg2_mc = eff_leg2_mc_fn(momentum_leg2.Pt(), momentum_leg2.Eta());
    weights.push_back(eff_leg1_data/eff_leg1_mc);
    weights.push_back(eff_leg2_data/eff_leg2_mc);
    return weights;
}

} // namespace detail

namespace Run2012ABCD {

namespace ETau {
    namespace Data {
        inline double electronEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.479) return detail::efficiency(pt, 22.9704, 1.0258, 1.26889,   1.31024, 1.06409);
            else                  return detail::efficiency(pt, 21.9816, 1.40993, 0.978597, 2.33144, 0.937552);
        }

        inline double tauEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.5) return detail::efficiency(pt, 18.538229, 0.651562, 0.324869, 13.099048, 0.902365);
            else                return detail::efficiency(pt, 18.756548, 0.230732, 0.142859, 3.358497,  0.851919);
        }
    }
    namespace MC {
        inline double electronEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.479) return detail::efficiency(pt, 21.7243, 0.619015, 0.739301, 1.34903, 1.02594);
            else                  return detail::efficiency(pt, 22.1217, 1.34054,  1.8885,   1.01855, 4.7241);
        }

        inline double tauEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.5) return detail::efficiency(pt, 18.605055, 0.264062, 0.139561, 4.792849,  0.915035);
            else                return detail::efficiency(pt, 18.557810, 0.280908, 0.119282, 17.749043, 0.865756);
        }
    }

    inline std::vector<double> CalculateWeights(const TLorentzVector& ele_momentum, const TLorentzVector& tau_momentum)
    {
        return detail::CalculateWeights(ele_momentum, tau_momentum,
                                        &Data::electronEfficiency, &Data::tauEfficiency,
                                        &MC::electronEfficiency, &MC::tauEfficiency);
    }
} // namespace ETau


namespace MuTau {
    namespace Data {
        inline double muonEfficiency(double pt, double eta) {
            if      (eta < -1.2) return detail::efficiency(pt, 15.9977, 7.64004e-05, 6.4951e-08,  1.57403, 0.865325);
            else if (eta < -0.8) return detail::efficiency(pt, 17.3974, 0.804001,    1.47145,     1.24295, 0.928198);
            else if (eta < 0.0)  return detail::efficiency(pt, 16.4307, 0.226312,    0.265553,    1.55756, 0.974462);
            else if (eta < 0.8)  return detail::efficiency(pt, 17.313,  0.662731,    1.3412,      1.05778, 1.26624);
            else if (eta < 1.2)  return detail::efficiency(pt, 16.9966, 0.550532,    0.807863,    1.55402, 0.885134);
            else                 return detail::efficiency(pt, 15.9962, 0.000106195, 4.95058e-08, 1.9991,  0.851294);
        }

        inline double tauEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.5) return detail::efficiency(pt, 18.604910, 0.276042, 0.137039, 2.698437, 0.940721);
            else                return detail::efficiency(pt, 18.701715, 0.216523, 0.148111, 2.245081, 0.895320);
        }
    }
    namespace MC {
        inline double muonEfficiency(double pt, double eta) {
            if      (eta < -1.2) return detail::efficiency(pt, 16.0051, 2.45144e-05, 4.3335e-09,  1.66134, 0.87045);
            else if (eta < -0.8) return detail::efficiency(pt, 17.3135, 0.747636,    1.21803,     1.40611, 0.934983);
            else if (eta < 0.0)  return detail::efficiency(pt, 15.9556, 0.0236127,   0.00589832,  1.75409, 0.981338);
            else if (eta < 0.8)  return detail::efficiency(pt, 15.9289, 0.0271317,   0.00448573,  1.92101, 0.978625);
            else if (eta < 1.2)  return detail::efficiency(pt, 16.5678, 0.328333,    0.354533,    1.67085, 0.91699);
            else                 return detail::efficiency(pt, 15.997,  7.90069e-05, 4.40036e-08, 1.66272, 0.884502);
        }

        inline double tauEfficiency(double pt, double eta) {
            if(fabs(eta) < 1.5)  return detail::efficiency(pt, 18.532997, 1.027880, 2.262950, 1.003322,  5.297292);
            else                 return detail::efficiency(pt, 18.212782, 0.338119, 0.122828, 12.577926, 0.893975);
        }
    }

    inline std::vector<double> CalculateWeights(const TLorentzVector& mu_momentum, const TLorentzVector& tau_momentum)
    {
        return detail::CalculateWeights(mu_momentum, tau_momentum,
                                        &Data::muonEfficiency, &Data::tauEfficiency,
                                        &MC::muonEfficiency, &MC::tauEfficiency);
    }
} // namespace MuTau

namespace TauTau {
    namespace Data {
        inline double tauEfficiency(double pt, double /*eta*/) {
            return ( 0.826969 * 0.5 * (TMath::Erf((pt - 42.2274) / 2. / 0.783258 /sqrt(pt)) + 1.) ) ;
        }
    }
    namespace MC {
        inline double tauEfficiency(double pt, double /*eta*/) {
            const double data_plateau = 0.826969 ;
            if (pt < 140.)     return ( 0.813769 * 0.5 * (TMath::Erf((pt - 39.9322) / 2. / 0.819354  /sqrt(pt)) + 1.) );
            else if (pt > 400) return data_plateau / 2.03467;
            else if (pt > 300) return data_plateau / 1.31593;
            else if (pt > 250) return data_plateau / 1.25698;
            else if (pt > 200) return data_plateau / 1.18941;
            else if (pt > 180) return data_plateau / 1.17448;
            else if (pt > 160) return data_plateau / 1.0964 ;
            else               return data_plateau / 1.09279;
        }
    }

    inline std::vector<double> CalculateWeights(const TLorentzVector& lead_tau_momentum,
                                                const TLorentzVector& sublead_tau_momentum)
    {
        return detail::CalculateWeights(lead_tau_momentum, sublead_tau_momentum,
                                        &Data::tauEfficiency, &Data::tauEfficiency,
                                        &MC::tauEfficiency, &MC::tauEfficiency);
    }
} // namespace TauTau

} // namesapce Run2012ABCD
} // namespace trigger
} // namespace Htautau_Summer13
} // namespace analysis
