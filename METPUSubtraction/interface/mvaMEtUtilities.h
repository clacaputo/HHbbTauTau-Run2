#ifndef RecoMET_METPUSubtraction_mvaMEtUtilities_h
#define RecoMET_METPUSubtraction_mvaMEtUtilities_h

#include <vector>
#include <utility>
#include <TLorentzVector.h>
#include "CommonMETData.h"

class mvaMEtUtilities 
{
 public:

  struct pfCandInfo 
  {
    pfCandInfo() : p4_(0.,0.,0.,0.), dZ_(0.) {}
    TLorentzVector p4_;
    double dZ_;
  };

  struct leptonInfo 
  {
    leptonInfo() : p4_(0.,0.,0.,0.), chargedFrac_(0.) {}
    TLorentzVector p4_;
    double chargedFrac_;
  };

  struct JetInfo 
  {
    JetInfo() : p4_(0.,0.,0.,0.), mva_(0.), neutralEnFrac_(0.) {}
    TLorentzVector p4_;
    double mva_;
    double neutralEnFrac_;  
  };

  mvaMEtUtilities();

  friend bool operator<(const JetInfo&, const JetInfo&);

  bool passesMVA(const TLorentzVector&, double);

  TLorentzVector leadJetP4(const std::vector<JetInfo>&);
  TLorentzVector subleadJetP4(const std::vector<JetInfo>&);
  unsigned numJetsAboveThreshold(const std::vector<JetInfo>&, double);

  std::vector<JetInfo> cleanJets(const std::vector<JetInfo>&, 
				 const std::vector<leptonInfo>&, double, double);

  std::vector<pfCandInfo> cleanPFCands(const std::vector<pfCandInfo>&, 
				       const std::vector<leptonInfo>&, double, bool);

  CommonMETData computePFCandSum(const std::vector<pfCandInfo>&, double, int);
  CommonMETData computeJetSum_neutral(const std::vector<JetInfo>&, bool);

  CommonMETData computePUMEt(const std::vector<pfCandInfo>&, const std::vector<JetInfo>&, double);
  
  CommonMETData computeNegPFRecoil   (const CommonMETData&, const std::vector<pfCandInfo>&, double);
  CommonMETData computeNegTrackRecoil(const CommonMETData&, const std::vector<pfCandInfo>&, double);
  CommonMETData computeNegNoPURecoil (const CommonMETData&, const std::vector<pfCandInfo>&, const std::vector<JetInfo>&, double);
  CommonMETData computeNegPUCRecoil  (const CommonMETData&, const std::vector<pfCandInfo>&, const std::vector<JetInfo>&, double);
  CommonMETData computeSumLeptons    (const std::vector<leptonInfo>& leptons, bool iCharged);
  void finalize(CommonMETData& metData);
 protected:

  TLorentzVector jetP4(const std::vector<JetInfo>&, unsigned);

  // cuts on jet Id. MVA output in bins of jet Pt and eta
  double mvaCut_[3][4][4]; 
};

#endif
