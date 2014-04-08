#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include "TMath.h"

TH1 *plot_MC;
TH1 *plot_data;
TH1 *plot_data_div_MC;
TFile *_i1, *_i2, *_o;
TCanvas *c,*c1,*c2;

void Fit_nPU(char _input1_filename[]  = "Pileup_2011_to_173692_LPLumiScale_68mb.root", 
	 char _input2_filename[] = "/cmshome/grippo/CMSSW_4_2_6/src/VHTauTau/TreeMaker/test/codice/batch_Analysis/18_10_2011/histo_DYmumu.root")
  
{
  _i1 = new TFile(_input1_filename,"read");
  
  _i2 = new TFile(_input2_filename,"read");
  
  char plotnametoy[60];
  char plotnamedata[60];
  sprintf (plotnametoy, "hNPU");
  sprintf (plotnamedata, "pileup");
  
  _i2->cd();
  plot_MC = (TH1F*)(_i2 -> FindObjectAny(plotnametoy));
  plot_MC->Sumw2();
  //plot_MC->Scale(1./plot_MC->Integral(1,50));
  //plot_MC->Scale(1./plot_MC->Integral(1,24));
  plot_MC->Scale(1./plot_MC->Integral());
  c=new TCanvas("c","c");
  c->cd();
  plot_MC->Draw();
    
  _i1->cd();
  plot_data = (TH1F*)(_i1 -> FindObjectAny(plotnamedata));
  plot_data->Sumw2();
  //plot_data->Scale(1./plot_data->Integral(1,50));
  plot_data->Scale(1./plot_data->Integral(1,24));
  //plot_data->Scale(1./plot_data->Integral());
  c1=new TCanvas("c1","c1");
  c1->cd();
  plot_data->Draw();

  plot_data_div_MC = (TH1F*)plot_data->Clone("plot_data_div_MC");
  plot_data_div_MC->Divide(plot_MC);
  c2=new TCanvas("c2","c2");
  c2->cd();
  plot_data_div_MC->Draw();
  
  _o = new TFile("reweight_21_10.root","RECREATE");
  plot_data_div_MC->Write();
  plot_MC->Write();
  plot_data->Write();
  _o->Close();
}
