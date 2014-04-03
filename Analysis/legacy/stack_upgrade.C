#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TProfile.h"
#include "THStack.h"
#include <string>
#include <sstream>

using namespace std;


int stack_upgrade() {

        gROOT->ProcessLine(".L ./tdrstyle.C");
        TH1::SetDefaultSumw2();
        setTDRStyle();


  const int nFiles = 5;

  Int_t r = 4;
  TString fileNames[nFiles] = {"data", "TTT","WZto3lNu","ZZto4L","DY"};

  Double_t weights[nFiles] = {1.00, 0.0074359927, 0.003447122, 0.00042929, 0.6319085498};

  TFile * f[nFiles];
  std::stringstream indexes;

  for(int iFile = 0; iFile < nFiles; iFile++)

  {
          indexes.str("");
          indexes << fileNames[iFile];
          std::string input_file="/nfs/scratch/fynu/lperrini/ZHtautauAnalysis/config/output_HISTOS_PU_OBSERVED/"+indexes.str()+".root";

          f[iFile] = TFile::Open(input_file.c_str());

          if(!f[iFile]) {
                  std::cerr << "Error: file " << input_file << " could not be opened." << std::endl;
                  return 1;
          }
          
          else std::cout << "File " << input_file << " succesfully opened!" << std::endl;
  }



  const int nHist1 = 29;


  TString histNames1[nHist1] = {"h_mu1Z_pt", "h_mu2Z_pt", "h_Z_lep1_eta", "h_Z_lep2_eta", "h_Z_lep1_phi", "h_Z_lep2_phi", "h_Zmass_mumu", "h_Zmass_ee","h_Zpt_mumu", "h_Zpt_ee","h_Zmass", "h_Zpt","h_H_mass", "h_H_pt","h_H_eta", "h_H_phi", "h_H_mass_type_1", "h_H_mass_type_2", "h_H_mass_type_3", "h_H_mass_type_4", "h_H_mass_type_5", "h_H_mass_type_6", "h_H_mass_type_7", "h_H_mass_type_8", "h_Tmass","h_H_lep1_eta","h_H_lep2_eta","h_H_lep1_phi","h_H_lep2_phi"};

  TH1F *                h_1d[nHist1][nFiles];

  for(int iFile = 0; iFile < nFiles; iFile++)

  {
          for(int iHist = 0; iHist < nHist1; iHist++)
          {
                  h_1d[iHist][iFile] = (TH1F*)f[iFile]->Get(histNames1[iHist]);
                  h_1d[iHist][iFile]->Rebin(r);
                  h_1d[iHist][iFile]->Scale(weights[iFile]);
          }
  }

  TCanvas* c1 = new TCanvas("c1","c1", 800,600);

  for(int iHist = 0; iHist < nHist1; iHist++)

  {
          THStack *hs = new THStack("hs","Stacked MC histograms");


          for(int iFile=1; iFile < nFiles; iFile++)

          {

                  h_1d[iHist][iFile]->SetLineWidth(0);

                  if(iFile == 1){ h_1d[iHist][iFile]->SetFillColor(33); h_1d[iHist][iFile]->SetLineColor(33);  }
                  else if(iFile == 2){  h_1d[iHist][iFile]->SetFillColor(2); h_1d[iHist][iFile]->SetLineColor(2); }
                  else if(iFile == 3){ h_1d[iHist][iFile]->SetFillColor(3); h_1d[iHist][iFile]->SetLineColor(3); }
                  else if(iFile == 4){ h_1d[iHist][iFile]->SetFillColor(5); h_1d[iHist][iFile]->SetLineColor(5); }
                  hs->Add(h_1d[iHist][iFile],"hist");

          }

          //h_1d[iHist][0]->SetMarkerStyle(21);
          //h_1d[iHist][0]->SetMarkerSize(0.7);

 TLegend* leg = new TLegend( 0.57, 0.65, 0.77, 0.92);
 leg->SetFillColor(0);
 leg->SetTextSize(0.035);

          leg->AddEntry(h_1d[iHist][0],"data 2011","p");
          leg->AddEntry(h_1d[iHist][1],"t#bar{t}","f");
          leg->AddEntry(h_1d[iHist][2],"WZ","f");
          leg->AddEntry(h_1d[iHist][3],"ZZ","f");
          leg->AddEntry(h_1d[iHist][4],"drell-yan","f");

TString lumist="4.4 fb^{-1}";
  TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
  ll->SetTextSize(0.03);
  ll->SetTextFont(42);
  ll->SetFillColor(0);
  ll->SetBorderSize(0);
  ll->SetMargin(0.01);
  ll->SetTextAlign(12); // align left
  TString text = "CMS Preliminary";
  ll->AddText(0.01,0.5,text);
  text = "#sqrt{s} = 7 TeV  L = ";
  text = text + lumist;
  //  ll->SetTextAlign(32); // align right
  ll->AddText(0.65, 0.6, text);

          hs->Draw();
          h_1d[iHist][0]->Draw("samePE01");
          leg->Draw("same");
          ll->Draw("same");

          c1->Print(histNames1[iHist]+"_all.png");
          c1->Print(histNames1[iHist]+"_all.eps");

          leg->Clear();
          hs->Clear();

  }

  return 0;

}

