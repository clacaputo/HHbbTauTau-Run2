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
#include <Rtypes.h>

class Print_SyncPlots {
public:
    typedef std::set<Int_t> EventSet;
    typedef std::vector<Int_t> EventVector;
    typedef std::map<Int_t, size_t> EventToEntryMap;
    typedef std::pair<size_t, size_t> EntryPair;
    typedef std::map<Int_t, EntryPair> EventToEntryPairMap;

    Print_SyncPlots(const std::string& _channel,
                    const std::string& _myGroup, const std::string& _myRootFile, const std::string& _myTree,
                    const std::string& _group, const std::string& _groupRootFile, const std::string& _groupTree)
        : channel(_channel), myGroup(_myGroup), myRootFile(_myRootFile), myTree(_myTree),
          group(_group), groupRootFile(_groupRootFile), groupTree(_groupTree)
    {
        std::cout << channel << std::endl;
        std::cout << myGroup << "  " << myRootFile << "  " << myTree << std::endl;
        std::cout << group << "  " << groupRootFile << "  " << groupTree << std::endl;

        Tmine = LoadTree(Fmine, myRootFile, myTree);
        Tother = LoadTree(Fother, groupRootFile, groupTree);

        CollectEvents(Tmine, Tother);

        std::cout << "Mine: " << Fmine->GetName() << std::endl;
        std::cout << "Other: " << Fother->GetName() << std::endl;

        file_name = std::string("PlotsDiff_") + channel + "_" + myGroup + "_" + group + ".pdf";
    }

    void Run()
    {
        using namespace std;

        canvas.Print((file_name+"[").c_str());

        //drawHistos("run", 200, 197000, 199000);
        //drawHistos("run", 300, 170000, 173000);
        //drawHistos("run", 200, 190000, 194000);
        //drawHistos("run", 200, 190000, 210000);
        drawHistos<Int_t>("npu", 50, 0, 50);
        drawHistos<Double_t>("pt_1", 100, 0, 100);
        drawHistos<Double_t>("eta_1", 60, -3, 3);
        drawHistos<Double_t>("iso_1", 200, -1, 1);
        drawHistos<Double_t>("pt_2", 100, 0, 100);
        drawHistos<Double_t>("eta_2", 60, -3, 3);
        drawHistos<Double_t>("iso_2", 200, -1, 1);
        drawHistos<Double_t>("mvis", 50, 0, 200);
        drawHistos<Double_t>("met", 20, 0, 200);
        drawHistos<Double_t>("metphi", 30, -3.5, 3.5);
        drawHistos<Double_t>("mvamet", 30, 0, 150);
        drawHistos<Double_t>("mvametphi", 35, -3.5, 3.5);
        drawHistos<Double_t>("mvacov00", 40, 0, 1000);
        drawHistos<Double_t>("mvacov01", 40, 0, 1000);
        drawHistos<Double_t>("mvacov10", 40, 0, 1000);
        drawHistos<Double_t>("mvacov11", 40, 0, 1000);
        //drawHistos<Double_t>("mt_1", 50, 0, 200);
        drawHistos<Double_t>("m_sv", 60, 0, 300);


        // Jets
        const std::vector<Int_t> my_njets = CollectValues<Int_t>(Tmine, "njets");
        const std::vector<Int_t> other_njets = CollectValues<Int_t>(Tother, "njets");
        const auto noJets_my = [&](size_t entry_id) -> bool {
            return my_njets.at(entry_id) == 0;
        };
        const auto atLeast1jet_my = [&](size_t entry_id) -> bool {
            return my_njets.at(entry_id) >= 1;
        };
        const auto atLeast2jets_my = [&](size_t entry_id) -> bool {
            return my_njets.at(entry_id) >= 2;
        };
        const auto noJets_other = [&](size_t entry_id) -> bool {
            return other_njets.at(entry_id) == 0;
        };
        const auto atLeast1jet_other = [&](size_t entry_id) -> bool {
            return other_njets.at(entry_id) >= 1;
        };
        const auto atLeast2jets_other = [&](size_t entry_id) -> bool {
            return other_njets.at(entry_id) >= 2;
        };

        drawHistos<Int_t>("njets", 5, -0.5, 4.5);
        drawHistos<Double_t>("jpt_1", 50, 0, 300, atLeast1jet_my, atLeast1jet_other, "njets>=1");
        drawHistos<Double_t>("jeta_1", 50, -5, 5, atLeast1jet_my, atLeast1jet_other, "njets>=1");
        drawHistos<Double_t>("jpt_2", 50, 0, 200, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        drawHistos<Double_t>("jeta_2", 50, -5, 5, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("mjj", 50, 0, 3000, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("jdeta", 50, 0, 10, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Int_t>("njetingap", 4, -0.5, 3.5, atLeast2jets_my, atLeast2jets_other, "njets>=2");

        //drawHistos<Double_t>("visjeteta", 100, 0, 10, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("ptvis", 100, 0, 500, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("jdphi", 100, 0, 3.5, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("dijetpt", 100, 0, 500, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("hdijetphi", 100, 0, 3.5, atLeast2jets_my, atLeast2jets_other, "njets>=2");
        //drawHistos<Double_t>("mva", 40, -1, 1.001, atLeast2jets_my, atLeast2jets_other, "njets>=2");

        drawHistos<Double_t>("mvis", 50, 0, 200, noJets_my, noJets_other, "njets==0");


        // b-jets
        const std::vector<Int_t> my_nbjets = CollectValues<Int_t>(Tmine, "nbtag");
        const std::vector<Int_t> other_nbjets = CollectValues<Int_t>(Tother, "nbtag");
        const auto atLeast1bjet_my = [&](size_t entry_id) -> bool {
            return my_nbjets.at(entry_id) >= 1;
        };
        const auto atLeast1bjet_other = [&](size_t entry_id) -> bool {
            return other_nbjets.at(entry_id) >= 1;
        };

        drawHistos<Int_t>("nbtag", 5, -0.5, 4.5);
        drawHistos<Double_t>("bpt", 50, 0, 200, atLeast1bjet_my, atLeast1bjet_other, "nbtag>=1");
        drawHistos<Double_t>("beta", 50, -5, 5, atLeast1bjet_my, atLeast1bjet_other, "nbtag>=1");

        drawHistos<Double_t>("puweight", 25, -.1, 0.9);
        //drawHistos<Double_t>("effweight", 30, .8, 1.1);
        //drawHistos<Double_t>("embeddedWeight", 50, 0, 1);
        //drawHistos<Double_t>("weight", 30, 0, 2.0);

        //with puweight
//        TString weight="(puweight)";
//        selection=weight;
//        drawHistos<Int_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"npv",50,0,50,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
//        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with trigger weight
//        weight="(effweight)";
//        selection=weight;
//        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mt_1",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);
//        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

//        weight="(weight)";
//        selection=weight;
//        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mt_1",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        //with mTcut
//        weight="(mt_1<50)";
//        selection=weight;
//        drawHistos<Double_t>(&C,filename,TString("inclusive * ")+weight,Tmine,Tother,"mvis",50,0,200,selection,myGroup,myRootFile,group,groupRootFile,mySel,groupSel);

        canvas.Print((file_name+"]").c_str());
    }

private:
    static TTree* LoadTree(std::shared_ptr<TFile>& file, const std::string& fileName, const std::string& treeName)
    {
        file = std::shared_ptr<TFile>(new TFile(fileName.c_str(), "READ"));
        TTree* tree = (TTree*)file->Get(treeName.c_str());
        if(!tree) {
            std::ostringstream ss;
            ss << "File " << fileName << " is empty.";
            throw std::runtime_error(ss.str());

        }
        tree->SetBranchStatus("*", 0);
        return tree;
    }

    template<typename VarType>
    void drawHistos(const std::string& var, int nbins, float xmin, float xmax)
    {
        const auto SelectAll = [](size_t entry_id) -> bool { return true; };
        drawHistos<VarType>(var, nbins, xmin, xmax, SelectAll, SelectAll, "All");
    }

    template<typename VarType, typename MySelector, typename OtherSelector>
    void drawHistos(const std::string& var, int nbins, float xmin, float xmax, const MySelector& my_selector,
                    const OtherSelector& other_selector, const std::string& selection_label)
    {
        try {

            std::shared_ptr<TH1F> Hmine_all(new TH1F(TString("Hmine") + var + "all","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_all(new TH1F(TString("Hother")+var + "all","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hmine_common(new TH1F(TString("Hmine")+var + "common","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_common(new TH1F(TString("Hother")+var + "common","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hmine_diff(new TH1F(TString("Hmine")+var + "diff","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_diff(new TH1F(TString("Hother")+var + "diff","",nbins,xmin,xmax));

            std::shared_ptr<TH2F> Hmine_vs_other(new TH2F(TString("Hmine_vs_other") + var, "", nbins, xmin, xmax,
                                                          nbins, xmin, xmax));


            const std::vector<VarType> my_values = CollectValues<VarType>(Tmine, var);
            const std::vector<VarType> other_values = CollectValues<VarType>(Tother, var);
            FillCommonHistograms<VarType>(my_values, other_values, my_selector, other_selector,
                                          *Hmine_common, *Hother_common, *Hmine_vs_other);
            FillInclusiveHistogram(my_values, my_selector, *Hmine_all);
            FillInclusiveHistogram(other_values, other_selector, *Hother_all);
            FillExclusiveHistogram(my_values, my_selector, *Hmine_diff, my_events_only_map);
            FillExclusiveHistogram(other_values, other_selector, *Hother_diff, other_events_only_map);

            DrawSuperimposedHistograms(Hmine_all, Hother_all, selection_label + " (all)", var);
            DrawSuperimposedHistograms(Hmine_common, Hother_common, selection_label + " (common)", var);
            DrawSuperimposedHistograms(Hmine_diff, Hother_diff, selection_label + " (diff)", var);
            Draw2DHistogram(Hmine_vs_other, selection_label, var);

        } catch(std::runtime_error& e){
            std::cerr << "ERROR: " << e.what() << std::endl;
        }
    }

    void DrawSuperimposedHistograms(std::shared_ptr<TH1F> Hmine, std::shared_ptr<TH1F> Hother,
                                    const std::string& selection_label, const std::string& var)
    {
        Hmine->SetTitle(selection_label.c_str());
        Hmine->GetYaxis()->SetTitle(selection_label.c_str());
        Hmine->GetXaxis()->SetTitle(var.c_str());
        Hmine->SetLineColor(1);
        Hmine->SetMarkerColor(1);
        Hmine->SetStats(0);

        Hother->GetYaxis()->SetTitle(selection_label.c_str());
        Hother->GetXaxis()->SetTitle(var.c_str());
        Hother->SetLineColor(2);
        Hother->SetMarkerColor(2);
        Hother->SetStats(0);

        TPad pad1("pad1","",0,0.2,1,1);
        TPad pad2("pad2","",0,0,1,0.2);

        pad1.cd();

        // Draw one histogram on top of the other
        if(Hmine->GetMaximum()>Hother->GetMaximum())
            Hmine->GetYaxis()->SetRangeUser(0,Hmine->GetMaximum()*1.1);
        else
            Hmine->GetYaxis()->SetRangeUser(0,Hother->GetMaximum()*1.1);
        Hmine->Draw("hist");
        Hother->Draw("histsame");
        DrawTextLabels(Hmine->Integral(0,Hmine->GetNbinsX()+1), Hother->Integral(0,Hother->GetNbinsX()+1));

        pad2.cd();

        // Draw the ratio of the historgrams
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
        HDiff->SetMarkerStyle(7);
        HDiff->SetMarkerColor(2);
        HDiff->Draw("histp");
        TLine line;
        line.DrawLine(HDiff->GetXaxis()->GetXmin(),1,HDiff->GetXaxis()->GetXmax(),1);

        canvas.Clear();
        pad1.Draw();
        pad2.Draw();

        PrintCanvas();
    }

    void Draw2DHistogram(std::shared_ptr<TH2F> Hmine_vs_other, const std::string& selection_label,
                         const std::string& var)
    {
        Hmine_vs_other->SetTitle(selection_label.c_str());
        Hmine_vs_other->GetXaxis()->SetTitle((var + "_mine").c_str());
        Hmine_vs_other->GetYaxis()->SetTitle((var + "_other").c_str());

        TPad pad1("pad1","", 0, 0, 1, 1);
        pad1.cd();
        Hmine_vs_other->Draw("colz");
        const size_t n_events = Hmine_vs_other->Integral(0, Hmine_vs_other->GetNbinsX() + 1,
                                                         0, Hmine_vs_other->GetNbinsY() + 1);
        DrawTextLabels(n_events, n_events);
        canvas.Clear();
        pad1.Draw();
        PrintCanvas();
    }

    void DrawTextLabels(size_t n_events_mine, size_t n_events_other)
    {
        TText TXmine;
        TXmine.SetTextColor(1);
        TXmine.SetTextSize(.04);
        TText TXother;
        TXother.SetTextColor(2);
        TXother.SetTextSize(.04);

        //Print the integrals of the histograms a the top
        //TXmine.DrawTextNDC(.2,.965,myGroup+"_"+myRootFile+": "+(long)(Hmine->Integral(0,Hmine->GetNbinsX()+1)));
        //TXother.DrawTextNDC(.2,.93,group+"_"+groupRootFile+": "+(long)(Hother->Integral(0,Hother->GetNbinsX()+1)));
        TXmine.DrawTextNDC(.23,.84,myGroup+" : " + n_events_mine);
        TXother.DrawTextNDC(.53,.84,group+": " + n_events_other);
    }

    template<typename VarType, typename MySelector, typename OtherSelector, typename Histogram, typename Histogram2D>
    void FillCommonHistograms(const std::vector<VarType>& my_values, const std::vector<VarType>& other_values,
                              const MySelector& my_selector, const OtherSelector& other_selector,
                              Histogram& my_histogram, Histogram& other_histogram, Histogram2D& histogram2D)
    {
        for(const auto& event_entry_pair : common_event_to_entry_pair_map) {
            const size_t my_entry = event_entry_pair.second.first;
            const size_t other_entry = event_entry_pair.second.second;
            if(!my_selector(my_entry) || !other_selector(other_entry))
                continue;
            const VarType& my_value = my_values.at(my_entry);
            const VarType& other_value = other_values.at(other_entry);
            my_histogram.Fill(my_value);
            other_histogram.Fill(other_value);
            histogram2D.Fill(my_value, other_value);
        }
    }

    template<typename VarType, typename Histogram, typename Selector>
    void FillInclusiveHistogram(const std::vector<VarType>& values, const Selector& selector, Histogram& histogram)
    {
        for(size_t n = 0; n < values.size(); ++n) {
            if(selector(n))
                histogram.Fill(values.at(n));
        }
    }

    template<typename VarType, typename Histogram, typename Selector>
    void FillExclusiveHistogram(const std::vector<VarType>& values, const Selector& selector,
                                Histogram& histogram, const EventToEntryMap& exclusive_event_to_entry_map)
    {
        for(const auto& event_entry : exclusive_event_to_entry_map) {
            if(selector(event_entry.second)) {
                const VarType& value = values.at(event_entry.second);
                histogram.Fill(value);
            }
        }
    }

    void CollectEvents(TTree* my_tree, TTree* other_tree)
    {
        my_events = CollectValues<Int_t>(my_tree, "evt");
        const EventSet my_events_set(my_events.begin(), my_events.end());
        other_events = CollectValues<Int_t>(other_tree, "evt");
        const EventSet other_events_set(other_events.begin(), other_events.end());

        EventSet intersection, my_events_only, other_events_only;
        EventSetIntersection(my_events_set, other_events_set, intersection);
        EventSetDifference(my_events_set, other_events_set, my_events_only);
        EventSetDifference(other_events_set, my_events_set, other_events_only);

        EventToEntryMap my_event_to_entry_map, other_event_to_entry_map;
        FillEventToEntryMap(my_events, my_event_to_entry_map);
        FillEventToEntryMap(other_events, other_event_to_entry_map);

        std::cout << "Our events" << std::endl;
        for(const auto& event_entry : my_event_to_entry_map) {
            if(intersection.count(event_entry.first)) {
                common_event_to_entry_pair_map[event_entry.first] =
                        std::pair<size_t, size_t>(event_entry.second, other_event_to_entry_map.at(event_entry.first));
            }
            if(my_events_only.count(event_entry.first)) {
                my_events_only_map[event_entry.first] = event_entry.second;
                std::cout << "eventId = " << event_entry.first << std::endl;
            }
        }

        std::cout << "Riccardo events" << std::endl;
        for(const auto& event_entry : other_event_to_entry_map) {
            if(other_events_only.count(event_entry.first)){
                other_events_only_map[event_entry.first] = event_entry.second;
                std::cout << "eventId = " << event_entry.first << std::endl;
            }
        }

        std::cout << "# my events = " << my_events.size() << ", " << "# my unique events = " << my_events_set.size()
                  << "\n# other events = " << other_events.size()
                  << ", # other unique events = " << other_events_set.size()
                  << "\n# common events = " << intersection.size() << std::endl;
    }

    static void EventSetIntersection(const EventSet& first_set, const EventSet& second_set, EventSet& intersection_set)
    {
        const size_t max_intersection_size = std::max(first_set.size(), second_set.size());
        EventVector intersection_vector(max_intersection_size);
        const auto iter = std::set_intersection(first_set.begin(), first_set.end(),
                                                second_set.begin(), second_set.end(),
                                                intersection_vector.begin());
        intersection_vector.resize(iter - intersection_vector.begin());
        intersection_set.clear();
        intersection_set.insert(intersection_vector.begin(), intersection_vector.end());
    }

    static void EventSetDifference(const EventSet& first_set, const EventSet& second_set, EventSet& diff_set)
    {
        EventVector diff_vector(first_set.size());
        const auto iter = std::set_difference(first_set.begin(), first_set.end(),
                                              second_set.begin(), second_set.end(),
                                              diff_vector.begin());
        diff_vector.resize(iter - diff_vector.begin());
        diff_set.clear();
        diff_set.insert(diff_vector.begin(), diff_vector.end());
    }

    static void FillEventToEntryMap(const EventVector& events, EventToEntryMap& event_to_entry_map)
    {
        for(size_t n = 0; n < events.size(); ++n)
            event_to_entry_map[events[n]] = n;
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

    void PrintCanvas()
    {
        canvas.Print(file_name.c_str());
    }

private:
    std::string channel, myGroup, myRootFile, myTree, group, groupRootFile, groupTree;
    std::shared_ptr<TFile> Fmine, Fother;
    TTree *Tmine, *Tother;

    EventVector my_events, other_events;
    EventToEntryMap my_events_only_map, other_events_only_map;
    EventToEntryPairMap common_event_to_entry_pair_map;

    TCanvas canvas;
    std::string file_name;
};
