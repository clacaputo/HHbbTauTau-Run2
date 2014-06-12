#ifndef TauAnalysis_SVfitStandAlone_svFitStandAloneAuxFunctions_h
#define TauAnalysis_SVfitStandAlone_svFitStandAloneAuxFunctions_h

#include <TH1.h>
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

namespace svFitStandalone
{
  //-----------------------------------------------------------------------------
  // define masses, widths and lifetimes of particles
  // relevant for computing values of likelihood functions in SVfit algorithm
  //
  // NOTE: the values are taken from
  //        K. Nakamura et al. (Particle Data Group),
  //        J. Phys. G 37, 075021 (2010)
  //
  const double electronMass = 0.51100e-3; // GeV
  const double electronMass2 = electronMass*electronMass;
  const double muonMass = 0.10566; // GeV
  const double muonMass2 = muonMass*muonMass; 
  
  const double chargedPionMass = 0.13957; // GeV
  const double chargedPionMass2 = chargedPionMass*chargedPionMass;
  const double neutralPionMass = 0.13498; // GeV
  const double neutralPionMass2 = neutralPionMass*neutralPionMass;

  const double tauLeptonMass = 1.77685; // GeV
  const double tauLeptonMass2 = tauLeptonMass*tauLeptonMass;
  const double tauLeptonMass3 = tauLeptonMass2*tauLeptonMass;
  const double tauLeptonMass4 = tauLeptonMass3*tauLeptonMass;
  const double cTauLifetime = 8.711e-3; // centimeters
  //-----------------------------------------------------------------------------

  inline double square(double x)
  {
    return x*x;
  }

  inline double cube(double x)
  {
    return x*x*x;
  }

  inline double fourth(double x)
  {
    return x*x*x*x;
  }

  /**
     \typedef SVfitStandalone::Vector
     \brief   spacial momentum vector (equivalent to Vector)
  */
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
  /**
     \typedef SVfitStandalone::LorentzVector
     \brief   lorentz vector (equivalent to reco::Candidate::LorentzVector)
  */
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  /// Determine Gottfried-Jackson angle from visible energy fraction X
  double gjAngleFromX(double, double, double, double, double);

  /// Determine visible tau rest frame energy given visible mass and neutrino mass
  double pVisRestFrame(double, double, double);

  /// Convert the Gottfried-Jackson angle into a lab frame opening angle
  double gjAngleToLabFrame(double, double, double);

  /// Determine the tau momentum in the lab frame given the rest frame assumptions
  double motherMomentumLabFrame(double, double, double, double, double);

  /// Determine the tau direction given our parameterization
  Vector motherDirection(const Vector&, double, double);

  /// Compute the tau four vector given the tau direction and momentum
  LorentzVector motherP4(const Vector&, double, double);

  /// Extract maximum, mean and { 0.84, 0.50, 0.16 } quantiles of distribution
  void extractHistogramProperties(const TH1*, const TH1*, double&, double&, double&, double&, double&, double&, double&, double&, int = 0);
}

#endif
