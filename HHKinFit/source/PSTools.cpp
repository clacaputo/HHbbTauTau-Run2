#include "../include/PSTools.h"

// #include <RooEllipse.h>
#include <Rtypes.h>
#include <TH1.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

//  Tools by Peter Schleper


void
PSTools::coutf(int w, TString s)
{    //  tools for formated std::cout
  std::string s1 = std::string(s.Data());
  coutf(w, s1);
}

void
PSTools::coutf(int w, int i)
{
  std::cout << std::setw(w)   << i ;
}

//void std::coutf(int w,string s) {std::cout << std::setw(w) << s ;}

void
PSTools::coutf(int w, std::string s)
{    // emulate "left" (bug in ROOT?
  int l = s.length();
  const std::string blank = "                       ";
  if (l <= w) {
    std::cout << s.substr(0, l);
    std::cout << blank.substr(0, w - l);
  }
  else {
    std::cout << s.substr(0, w);
  }
}

void
PSTools::coutf(int w, double x)
{
  std::cout << std::setw(w) << x;
}

void
PSTools::coutf(int w, int p, double x)
{
  printf("%*.*f ", w, p, x);
  //std::cout << std::setiosflags(std::ios::fixed) << std::setw(w) << std::setprecision(p) << x ; 
  //std::cout << std::setprecision(6) << std::resetiosflags(std::ios::fixed); 
}

