#ifndef RecoMET_METPUSubtraction_PFMETAlgorithmMVA_h
#define RecoMET_METPUSubtraction_PFMETAlgorithmMVA_h

/** \class PFMETAlgorithmMVA
 *
 * MVA based algorithm for computing the particle-flow missing Et
 *
 * \authors Phil Harris, CERN
 *          Christian Veelken, LLR
 *
 */

#include <string>
#include <vector>
#include <ostream>

#include <TMatrixD.h>

//#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "GBRForest.h"

//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"

#include "mvaMEtUtilities.h"



class PFMETAlgorithmMVA 
{
 public:

  PFMETAlgorithmMVA(double dZcut);
  ~PFMETAlgorithmMVA();

  void initialize(const std::string& inputFileNameU, const std::string& inputFileNameDPhi,
                  const std::string& inputFileNameCovU1, const std::string& inputFileNameCovU2);

  void setHasPhotons(bool hasPhotons) { hasPhotons_ = hasPhotons; }

  void setInput(const std::vector<mvaMEtUtilities::leptonInfo>&,
		const std::vector<mvaMEtUtilities::JetInfo>&,
		const std::vector<mvaMEtUtilities::pfCandInfo>&,
        const std::vector<TVector3>&);

  void evaluateMVA();

  TLorentzVector getMEt()    const { return mvaMEt_;    }
  const TMatrixD&                getMEtCov() const { return mvaMEtCov_; }

  double getU()     const { return mvaOutputU_;    }
  double getDPhi()  const { return mvaOutputDPhi_;  }
  double getCovU1() const { return mvaOutputCovU1_; }
  double getCovU2() const { return mvaOutputCovU2_; }
  
  void print(std::ostream&) const;

 protected:

  void setInput(double, double, double,
		double, double, double,
		double, double, double,
		double, double, double,
		double, double, double,
		double, double, double,
		double, double, double,
		double, double, 
		double);

  void evaluateU();
  void evaluateDPhi();
  void evaluateCovU1();
  void evaluateCovU2();

  mvaMEtUtilities utils_;
    
  std::string mvaNameU_;
  std::string mvaNameDPhi_;
  std::string mvaNameCovU1_;
  std::string mvaNameCovU2_;

//  int    mvaType_;
//  bool   is42_;
//  bool   isOld42_;
  bool   hasPhotons_;
  double dZcut_;

  Float_t pfSumEt_;
  Float_t pfU_;
  Float_t pfPhi_;
  Float_t tkSumEt_;
  Float_t tkU_;
  Float_t tkPhi_;
  Float_t npuSumEt_;
  Float_t npuU_;
  Float_t npuPhi_;
  Float_t puSumEt_;
  Float_t puMEt_;
  Float_t puPhi_;
  Float_t pucSumEt_;
  Float_t pucU_;
  Float_t pucPhi_;
  Float_t jet1Pt_;
  Float_t jet1Eta_;
  Float_t jet1Phi_;
  Float_t jet2Pt_;
  Float_t jet2Eta_;
  Float_t jet2Phi_;
  Float_t numJetsPtGt30_;
  Float_t numJets_;
  Float_t numVertices_;

  Float_t* mvaInputU_;
  Float_t* mvaInputDPhi_;
  Float_t* mvaInputCovU1_;
  Float_t* mvaInputCovU2_;
  
  Float_t mvaOutputU_;
  Float_t mvaOutputDPhi_;
  Float_t mvaOutputCovU1_;
  Float_t mvaOutputCovU2_;

  double sumLeptonPx_;
  double sumLeptonPy_;
  double chargedSumLeptonPx_;
  double chargedSumLeptonPy_;

  TLorentzVector mvaMEt_;
  TMatrixD mvaMEtCov_;

  const GBRForest* mvaReaderU_;
  const GBRForest* mvaReaderDPhi_;
  const GBRForest* mvaReaderCovU1_;
  const GBRForest* mvaReaderCovU2_;

//  bool loadMVAfromDB_;

//  edm::ParameterSet cfg_;
};
#endif
