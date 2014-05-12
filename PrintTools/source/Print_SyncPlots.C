/*!
 * \file Print_SyncPlots.C
 * \brief Print control plots that were selected to synchronize produced tree-toople.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-28 created
 */

#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TText.h>
#include <TLine.h>
#include <TPad.h>

class Print_SyncPlots {
public:
    Print_SyncPlots(const std::string& _channel, const std::string& _myGroup, const std::string& _myPath,
                    const std::string& _myRootFile, const std::string& _myTree, const std::string& _group,
                    const std::string& _groupPath, const std::string& _groupRootFile, const std::string& _groupTree,
                    bool _plotOnlyMatchedEvents, const std::string& _mySel="", const std::string& _groupSel="")
        : channel(_channel), myGroup(_myGroup), myPath(_myPath), myRootFile(_myRootFile), myTree(_myTree),
          group(_group), groupPath(_groupPath), groupRootFile(_groupRootFile), groupTree(_groupTree), mySel(_mySel),
          groupSel(_groupSel), plotOnlyMatchedEvents(_plotOnlyMatchedEvents)
    {}

    void Run()
    {
        using namespace std;

        cout<<channel<<endl;
        cout<<myGroup<<"  "<<myPath<<"  "<<myRootFile<<"  "<<myTree<<"  "<<mySel<<endl;
        cout<<group<<"  "<<groupPath<<"  "<<groupRootFile<<"  "<<groupTree<<" "<<groupSel<<endl;

        if(mySel.CompareTo("")==0)mySel="1";
        if(groupSel.CompareTo("")==0)groupSel="1";


        TFile Fmine(myPath+"/"+myRootFile+".root");
        TTree*Tmine=(TTree*)Fmine.Get(myTree.Data());
        TFile Fother(groupPath+"/"+groupRootFile+".root");
        TTree*Tother=(TTree*)Fother.Get(groupTree.Data());
        if(!Tmine) {
            std::ostringstream ss;
            ss << "File " << Fmine.GetName() << " is empty.";
            throw std::runtime_error(ss.str());

        }
        if(!Tother) {
            std::ostringstream ss;
            ss << "File " << Fother.GetName() << " is empty.";
            throw std::runtime_error(ss.str());
        }

        if(plotOnlyMatchedEvents) {
            Tmine->SetBranchStatus("*", 0);
            Tother->SetBranchStatus("*", 0);
            CollectEvents(Tmine, Tother);
        }
        ////////////////

        cout<<"Mine: "<<Fmine.GetName()<<endl;
        cout<<"Other: "<<Fother.GetName()<<endl;
        TCanvas C;

        //TString filename=TString("PlotsDiff_")+myGroup+"_"+myRootFile+"_"+group+"_"+groupRootFile+".pdf";
        TString filename=TString("PlotsDiff_")+channel+"_"+myGroup+"_"+group;
        if(plotOnlyMatchedEvents)
            filename += "_sameEvents";
        filename += ".pdf";
        C.Print(filename+"[");


        //inclusive
        TString selection="1";
        //drawHistos(&C,filename,"inclusive",Tmine,Tother,"run",200,197000,199000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,"inclusive",Tmine,Tother,"run",300,170000,173000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,"inclusive",Tmine,Tother,"run",200,190000,194000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,"inclusive",Tmine,Tother,"run",200,190000,210000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Int_t>(&C,filename,"inclusive",Tmine,Tother,"npu",50,0,50,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"pt_1",100,0,100,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"eta_1",60,-3,3,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"iso_1",60,-.02,0.12,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"pt_2",100,0,100,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"eta_2",60,-3,3,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"met",20,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"metphi",30,-3.5,3.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvamet",30,0,150,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvametphi",35,-3.5,3.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvacov00",40,0,1000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvacov01",40,0,1000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvacov10",40,0,1000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mvacov11",40,0,1000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"mt_1",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"m_sv",60,0,300,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);


        ///Jets
        drawHistos<Int_t>(&C,filename,"inclusive",Tmine,Tother,"njets",5,-0.5,4.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        selection="(njets>=1)";
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"jpt_1",50,0,300,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"jeta_1",50,-5,5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        selection="(njets>=2)";
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"jpt_2",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"jeta_2",50,-5,5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"mjj",50,0,3000,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"jdeta",50,0,10,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Int_t>(&C,filename,selection,Tmine,Tother,"njetingap",4,-0.5,3.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //drawHistos(&C,filename,selection,Tmine,Tother,"visjeteta",100,0,10,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,selection,Tmine,Tother,"ptvis",100,0,500,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,selection,Tmine,Tother,"jdphi",100,0,3.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,selection,Tmine,Tother,"dijetpt",100,0,500,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,selection,Tmine,Tother,"hdijetphi",100,0,3.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        //drawHistos(&C,filename,selection,Tmine,Tother,"mva",40,-1,1.001,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);


        /////b-jets
        selection="1";
        drawHistos<Int_t>(&C,filename,"inclusive",Tmine,Tother,"nbtag",5,-0.5,4.5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        selection="(nbtag>=1)";
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"bpt",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,selection,Tmine,Tother,"beta",50,-5,5,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        selection="1";
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"puweight",25,-.1,0.9,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"effweight",30,.8,1.1,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"embeddedWeight",50,0,1,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,"inclusive",Tmine,Tother,"weight",30,0,2.0,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with puweight
        TString weight="(puweight)";
        selection=weight;
        drawHistos<Int_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"npv",50,0,50,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with trigger weight
        weight="(effweight)";
        selection=weight;
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mt_1",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        weight="(weight)";
        selection=weight;
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mt_1",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with 0jet
        weight="(njets==0)";
        selection=weight;
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with mTcut
        weight="(mt_1<50)";
        selection=weight;
        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        C.Print(filename+"]");
    }

private:
    template<typename VarType>
    void drawHistos(TCanvas * C, TString filename, TString category, TTree* Tmine, TTree* Tother, TString var,
                    int nbins, float xmin, float xmax, TString selection, TString myGroup, TString myRootFile,
                    TString group, TString groupRootFile,TString mySel,TString groupSel)
    {
        try {

            //   cout<<Tmine->GetName()<<" "<<myGroup<<" "<<myRootFile<<" "<<mySel<<" "<<endl;
            //   cout<<Tother->GetName()<<" "<<group<<" "<<groupRootFile<<" "<<groupSel<<" "<<endl;
            //   cout<<var<<" "<<nbins<<" "<<xmin<<" "<<xmax<<" "<<selection<<endl;


            std::unique_ptr<TH1F> Hmine(new TH1F(TString("Hmine")+var,"",nbins,xmin,xmax));
            Hmine->GetYaxis()->SetTitle(category);
            Hmine->GetXaxis()->SetTitle(var);
            Hmine->SetLineColor(1);
            Hmine->SetMarkerColor(1);
            Hmine->SetStats(0);

            std::unique_ptr<TH1F> Hother(new TH1F(TString("Hother")+var,"",nbins,xmin,xmax));
            Hother->GetYaxis()->SetTitle(category);
            Hother->GetXaxis()->SetTitle(var);
            Hother->SetLineColor(2);
            Hother->SetMarkerColor(2);
            Hother->SetStats(0);

            std::unique_ptr<TH2F> Hmine_vs_other(new TH2F(TString("Hmine_vs_other") + var, "", nbins, xmin, xmax,
                                                          nbins, xmin, xmax));
            Hmine_vs_other->GetXaxis()->SetTitle(var + "_mine");
            Hmine_vs_other->GetYaxis()->SetTitle(var + "_other");

            TText TXmine;
            TXmine.SetTextColor(1);
            TXmine.SetTextSize(.04);
            TText TXother;
            TXother.SetTextColor(2);
            TXother.SetTextSize(.04);


            if(plotOnlyMatchedEvents) {
                const std::vector<VarType> my_values = CollectValues<VarType>(Tmine, var.Data());
                const std::vector<VarType> other_values = CollectValues<VarType>(Tother, var.Data());
                FillHistograms<VarType>(my_values, other_values, *Hmine, *Hother, *Hmine_vs_other);
            } else {
                Tmine->Draw(var+">>"+Hmine->GetName(),selection+"*("+mySel+")");
                Tother->Draw(var+">>"+Hother->GetName(),selection+"*("+groupSel+")");
            }

            TPad pad1("pad1","",0,0.2,1,1);
            TPad pad2("pad2","",0,0,1,0.2);

            ////////////////////////////////////////////
            pad1.cd();

            ////Draw one histogram on top of the other
            if(Hmine->GetMaximum()>Hother->GetMaximum())
            Hmine->GetYaxis()->SetRangeUser(0,Hmine->GetMaximum()*1.1);
            else Hmine->GetYaxis()->SetRangeUser(0,Hother->GetMaximum()*1.1);
            Hmine->SetTitle(selection);
            Hmine->Draw("hist");
            Hother->Draw("histsame");

            //Print the integrals of the histograms a the top
            //TXmine.DrawTextNDC(.2,.965,myGroup+"_"+myRootFile+": "+(long)(Hmine->Integral(0,Hmine->GetNbinsX()+1)));
            //TXother.DrawTextNDC(.2,.93,group+"_"+groupRootFile+": "+(long)(Hother->Integral(0,Hother->GetNbinsX()+1)));
            TXmine.DrawTextNDC(.23,.84,myGroup+" : "+(long)(Hmine->Integral(0,Hmine->GetNbinsX()+1)));
            TXother.DrawTextNDC(.53,.84,group+": "+(long)(Hother->Integral(0,Hother->GetNbinsX()+1)));

            ////////////////////////////////////////////
            pad2.cd();

            //   ///Draw the difference of the historgrams
            //   TH1F*HDiff=(TH1F*)Hmine->Clone("HDiff");
            //   HDiff->Add(Hother,-1);
            //   int max= abs(HDiff->GetMaximum())>abs( HDiff->GetMinimum()) ?   abs(HDiff->GetMaximum()): abs( HDiff->GetMinimum());
            //   HDiff->GetYaxis()->SetRangeUser(-2*(max>0?max:1),2*(max>0?max:1));
            //   HDiff->Draw("hist");
            //   TLine line;
            //   line.DrawLine(HDiff->GetXaxis()->GetXmin(),0,HDiff->GetXaxis()->GetXmax(),0);

            ///Draw the ratio of the historgrams
            std::unique_ptr<TH1F> HDiff((TH1F*)Hother->Clone("HDiff"));
            HDiff->Divide(Hmine.get());
            ///HDiff->GetYaxis()->SetRangeUser(0.9,1.1);
            HDiff->GetYaxis()->SetRangeUser(0.9,1.1);
            //HDiff->GetYaxis()->SetRangeUser(0.98,1.02);
            //HDiff->GetYaxis()->SetRangeUser(0.,2.0);
            HDiff->GetYaxis()->SetNdivisions(3);
            HDiff->GetYaxis()->SetLabelSize(0.1);
            HDiff->GetYaxis()->SetTitleSize(0.1);
            HDiff->GetYaxis()->SetTitleOffset(0.5);
            //HDiff->GetYaxis()->SetTitle(myGroup + " / " + group);
            HDiff->GetYaxis()->SetTitle("Ratio");
            HDiff->GetXaxis()->SetNdivisions(-1);
            HDiff->GetXaxis()->SetTitle("");
            HDiff->GetXaxis()->SetLabelSize(0.0001);
            HDiff->SetMarkerColor(2);
            HDiff->Draw("histp");
            TLine line;
            line.DrawLine(HDiff->GetXaxis()->GetXmin(),1,HDiff->GetXaxis()->GetXmax(),1);


            C->Clear();
            pad1.Draw();
            pad2.Draw();

            C->Print(filename);

            if(plotOnlyMatchedEvents) {
                pad1.Clear();
                pad1.cd();
                Hmine_vs_other->Draw();
                TXmine.DrawTextNDC(.23,.84,myGroup+" : "+(long)(Hmine->Integral(0,Hmine->GetNbinsX()+1)));
                TXother.DrawTextNDC(.53,.84,group+": "+(long)(Hother->Integral(0,Hother->GetNbinsX()+1)));
                C->Clear();
                pad1.Draw();
                C->Print(filename);
            }

        } catch(std::runtime_error& e){
            std::cerr << "ERROR: " << e.what() << std::endl;
        }
    }

    template<typename VarType, typename Histogram, typename Histogram2D>
    void FillHistograms(const std::vector<VarType>& my_values, const std::vector<VarType>& other_values,
                        Histogram& my_histogram, Histogram& other_histogram, Histogram2D& histogram2D)
    {
        for(const auto& event_entry_pair : event_to_entry_pair_map) {
            const size_t my_entry = event_entry_pair.second.first;
            const size_t other_entry = event_entry_pair.second.second;
            const VarType& my_value = my_values.at(my_entry);
            const VarType& other_value = other_values.at(other_entry);
            my_histogram.Fill(my_value);
            other_histogram.Fill(other_value);
            histogram2D.Fill(my_value, other_value);
        }
    }

    void CollectEvents(TTree* my_tree, TTree* other_tree)
    {
        my_events = CollectValues<Int_t>(my_tree, "evt");
        std::set<Int_t> my_events_set(my_events.begin(), my_events.end());
        other_events = CollectValues<Int_t>(other_tree, "evt");
        std::set<Int_t> other_events_set(other_events.begin(), other_events.end());
        const size_t max_intersection_size = std::max(my_events_set.size(), other_events_set.size());
        std::vector<Int_t> intersection_vector(max_intersection_size);
        const auto iter = std::set_intersection(my_events_set.begin(), my_events_set.end(),
                                                other_events_set.begin(), other_events_set.end(),
                                                intersection_vector.begin());
        intersection_vector.resize(iter - intersection_vector.begin());
        std::cout << "# my events = " << my_events.size() << ", " << "# my unique events = " << my_events_set.size()
                  << "\n# other events = " << other_events.size()
                  << ", # other unique events = " << other_events_set.size()
                  << "\n # common events = " << intersection_vector.size() << std::endl;
        intersection.clear();
        intersection.insert(intersection_vector.begin(), intersection_vector.end());
        std::map<Int_t, size_t> my_event_to_entry_map, other_event_to_entry_map;
        FillEventToEntryMap(my_events, intersection, my_event_to_entry_map);
        FillEventToEntryMap(other_events, intersection, other_event_to_entry_map);
        for(const auto& event_entry : my_event_to_entry_map) {
            event_to_entry_pair_map[event_entry.first] =
                    std::pair<size_t, size_t>(event_entry.second, other_event_to_entry_map.at(event_entry.first));
        }
    }

    void FillEventToEntryMap(const std::vector<Int_t>& events, const std::set<Int_t>& intersection,
                             std::map<Int_t, size_t>& event_to_entry_map)
    {
        for(size_t n = 0; n < events.size(); ++n) {
            if(intersection.count(events[n]))
                event_to_entry_map[events[n]] = n;
        }
    }

    template<typename VarType>
    std::vector<VarType> CollectValues(TTree* tree, const std::string& name)
    {
        EnableBranch(tree, name, true);
        std::vector<VarType> result;
        VarType value;
        tree->SetBranchAddress(name.c_str(), &value);
        const Long64_t N = tree->GetEntries();
        for(Long64_t n = 0; n < N;++n) {
            if(tree->GetEntry(n) < 0)
                throw std::runtime_error("error while reading tree.");
            result.push_back(value);
        }
        EnableBranch(tree, name, false);
        return result;
    }

    void EnableBranch(TTree* tree, const std::string& name, bool enable)
    {
        UInt_t n_found = 0;
        tree->SetBranchStatus(name.c_str(), enable, &n_found);
        if(n_found != 1) {
            std::ostringstream ss;
            ss << "Branch '" << name << "' is not found.";
            throw std::runtime_error(ss.str());
        }
    }


private:
    TString channel, myGroup, myPath, myRootFile, myTree, group, groupPath, groupRootFile, groupTree, mySel, groupSel;
    bool plotOnlyMatchedEvents;
    std::vector<Int_t> my_events, other_events;
    std::set<Int_t> intersection;
    std::map<Int_t, std::pair<size_t, size_t> > event_to_entry_pair_map;
};
