#ifndef MVASELECTIONS_H
#define MVASELECTIONS_H

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatTree.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"

#include <TLorentzVector.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class MVASelections
{
public:

    MVASelections(const std::string& _mvaMethodName, const std::string& _mvaXMLfile):mvaMethodName(_mvaMethodName),
    mvaXMLfile(_mvaXMLfile)
    {
        TMVA::Tools::Instance();

        reader = new TMVA::Reader("Silent");
        reader->AddVariable("pt_mu", &var1);
        reader->AddVariable("pt_tau", &var2);
        reader->AddVariable("pt_b1", &var3);
        reader->AddVariable("pt_b2", &var4);
        reader->AddVariable("DR_bb", &var5);
        reader->AddVariable("DPhi_BBMET", &var6);
        reader->AddVariable("DR_ll", &var7);
        reader->AddVariable("Pt_Htt", &var8);
        reader->AddVariable("DR_HBBHTT", &var9);
        reader->AddVariable("Pt_Hbb", &var10);
        reader->AddVariable("DeltaPhi_METTT", &var11 );
        reader->AddVariable("PtH", &var12);
        reader->AddVariable("mT2", &var13);
        reader->BookMVA(mvaMethodName.c_str(),mvaXMLfile.c_str());
    }

    double GetMVA (const TLorentzVector& leg1_momentum, const TLorentzVector& leg2_momentum,
                        const TLorentzVector& b1_momentum, const TLorentzVector& b2_momentum, const TLorentzVector& MET){

        TLorentzVector h_tt   = leg1_momentum + leg2_momentum;
        TLorentzVector h_bb   = b1_momentum   + b2_momentum;
        TLorentzVector H      = h_tt + h_bb;

        var1  = leg1_momentum.Pt();
        var2  = leg2_momentum.Pt();
        var3  = b1_momentum.Pt();
        var4  = b2_momentum.Pt();
        var5  = b1_momentum.DeltaR(b2_momentum);
        var6  = h_bb.DeltaPhi(MET);
        var7  = leg1_momentum.DeltaR(leg2_momentum);
        var8  = h_tt.Pt();
        var9  = h_bb.DeltaR(h_tt);
        var10 = h_bb.Pt();
        var11 = MET.DeltaPhi(h_tt);
        var12 = H.Pt();
        var13 = analysis::Calculate_MT(leg1_momentum,MET.Pt(),MET.Phi());

        return reader->EvaluateMVA( mvaMethodName.c_str() );

    }


private:
    TMVA::Reader* reader;
    std::string mvaMethodName;
    std::string mvaXMLfile;
    float var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13;
};

#endif // MVASELECTIONS_H
