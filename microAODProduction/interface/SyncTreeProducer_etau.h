// system include files
#include <memory>
#include <vector>

// user include files

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "BaseEDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TTree.h"
#include "Math/VectorUtil.h"


namespace analysis {

class SyncAnalyzerData_eTau : public analysis::BaseEDAnalyzerData {
public:
    SyncAnalyzerData_eTau(std::shared_ptr<TFile> outputFile) : BaseEDAnalyzerData(outputFile) {}
    SyncAnalyzerData_eTau(const std::string& outputFileName) : BaseEDAnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};
}

//
// class declaration
//

class SyncTreeProducer_etau: public BaseEDAnalyzer {
   public:
      explicit SyncTreeProducer_etau(const edm::ParameterSet&);
      ~SyncTreeProducer_etau();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0,
              TRUE_PROMPT_ELECTRON,
              TRUE_ELECTRON_FROM_TAU,
              TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      virtual analysis::SyncAnalyzerData_eTau& GetAnaData() override { return anaData; }
      virtual analysis::CandidateV2Ptr SelectHiggs(analysis::CandidateV2PtrVector& higgses) override;

      // ----------member data --------------------------

      Run2::SyncTree syncTree;
      analysis::SyncAnalyzerData_eTau anaData;
};