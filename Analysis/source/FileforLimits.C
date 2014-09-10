/*!
 * \class FileForLimits FileForLimits.C
 * \brief Class to copy Rebecca's value to compute limits.
 *
 *
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2013-09-10 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "TH1D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TKey.h"
#include "TSystem.h"
#include "TTree.h"
#include "../include/BaseFlatTreeAnalyzer.h"

class FileForLimits{
public:
    FileForLimits(const std::string& _our_File_name , const std::string& _Rebecca_File_name,
                  const std::string& limits_fileName)
        : our_File_name(_our_File_name), Rebecca_File_name(_Rebecca_File_name),
          outputFile(new TFile(limits_fileName.c_str(),"RECREATE"))
    {

    }


    void Run()
    {
        outputFile->cd();
        CopyFile(Rebecca_File_name.c_str(),true);
        std::cout<< "Rebecca copied" << std::endl;
        CopyFile(our_File_name.c_str(),true);
        std::cout<< "mine copied" << std::endl;
        outputFile->ls();
        delete outputFile;
    }

private:
    void CopyDir(TDirectory *source, bool createDir = true) {
       //copy all objects and subdirs of directory source as a subdir of the current directory
       source->ls();
       TDirectory *savdir = gDirectory;
       TDirectory *adir;
       std::cout<< "created dir from mine file" << std::endl;
       if (createDir)
           adir = savdir->mkdir(source->GetName());
       else{
           adir = static_cast<TDirectory*>(savdir->Get(source->GetName()));
           std::cout<< "got Dir; name: " << source->GetName() <<std::endl;
       }
       adir->cd();
       std::cout<< "entered in adir " <<std::endl;
       //loop on all entries of this directory
       TKey *key;
       TIter nextkey(source->GetListOfKeys());
       while ((key = (TKey*)nextkey())) {
          const char *classname = key->GetClassName();
          TClass *cl = gROOT->GetClass(classname);
          if (!cl) continue;
          if (cl->InheritsFrom("TDirectory")) {
             source->cd(key->GetName());
             TDirectory *subdir = gDirectory;
             adir->cd();
             CopyDir(subdir);
             adir->cd();
          } else if (cl->InheritsFrom("TTree")) {
             TTree *T = (TTree*)source->Get(key->GetName());
             adir->cd();
             TTree *newT = T->CloneTree();
             newT->Write();
          } else {
             source->cd();
             TObject *obj = key->ReadObj();
             adir->cd();
             obj->Write();
             delete obj;
         }
      }
      adir->SaveSelf(kTRUE);
      savdir->cd();
    }
    void CopyFile(const char *fname, bool createDir = true) {
       //Copy all objects and subdirs of file fname as a subdir of the current directory
       TDirectory *target = gDirectory;
       TFile *f = TFile::Open(fname);
       if (!f || f->IsZombie()) {
          printf("Cannot copy file: %s\n",fname);
          target->cd();
          return;
       }
       target->cd();
       std::cout<< "I am in CopyFile for mine file" << std::endl;
       CopyDir(f,createDir);
       delete f;
       target->cd();
    }

    std::string our_File_name;
    std::string Rebecca_File_name;
    TFile* outputFile;

};
