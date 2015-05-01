/*!
 * \file MergeRootFiles.C
 * \brief Code to merge two root files taking first file as a reference and adding missing info from the second.
 *
 * Some parts of the code are taken from copyFile.C written by Rene Brun.
 *
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

#include <iostream>

#include <TROOT.h>
#include <TKey.h>
#include <TSystem.h>
#include <TTree.h>
#include <memory>
#include "AnalysisBase/include/RootExt.h"

class MergeRootFiles {
public:
    MergeRootFiles(const std::string& originalFileName , const std::string& referenceFileName,
                  const std::string& outputFileName)
        : originalFile(root_ext::OpenRootFile(originalFileName)),
          referenceFile(root_ext::OpenRootFile(referenceFileName)),
          outputFile(root_ext::CreateRootFile(outputFileName)) {}

    void Run()
    {
        std::cout << "Copying original file..." << std::endl;
        CopyDirectory(originalFile.get(), outputFile.get(), false);
        std::cout << "Copying reference file..." << std::endl;
        CopyDirectory(referenceFile.get(), outputFile.get(), true);
        std::cout << "Original and reference files has been merged." << std::endl;
    }

private:
    /// Copy all objects and subdirs of the source directory to the destination directory.
    static void CopyDirectory(TDirectory *source, TDirectory *destination, bool isReference)
    {
        std::cout<< "CopyDir. Current direcotry: " << destination->GetName() << std::endl;

        TIter nextkey(source->GetListOfKeys());
        for(TKey* key; (key = (TKey*)nextkey());) {
            //std::cout << "Processing key: " << key->GetName() << std::endl;
            const char *classname = key->GetClassName();
            TClass *cl = gROOT->GetClass(classname);
            if (!cl) continue;
            bool objectWritten = false;
            if (cl->InheritsFrom("TDirectory")) {
                TDirectory *subdir_source = static_cast<TDirectory*>(source->Get(key->GetName()));
                TDirectory *subdir_destination;
                if(isReference) {
                    subdir_destination = static_cast<TDirectory*>(destination->Get(subdir_source->GetName()));
                    if(!subdir_destination) {
                        std::cout << "Skipping reference directory '" << subdir_source->GetName()
                                  << "', which is missing in the origin file.";
                        continue;
                    }

                } else
                    subdir_destination = destination->mkdir(subdir_source->GetName());

                CopyDirectory(subdir_source, subdir_destination, isReference);
            } else if(destination->Get(key->GetName())) {

            } else if (cl->InheritsFrom("TTree")) {
                TTree *T = (TTree*)source->Get(key->GetName());
                TTree *newT = T->CloneTree();
                destination->WriteTObject(newT, key->GetName(), "WriteDelete");
                objectWritten = true;
            } else {
                std::unique_ptr<TObject> original_obj(key->ReadObj());
                std::unique_ptr<TObject> obj(original_obj->Clone());
                destination->WriteTObject(obj, key->GetName(), "WriteDelete");
                objectWritten = true;
            }

            if(objectWritten && isReference)
                std::cout << "Object '" << key->GetName() << "' taken from the reference file for '"
                          << destination->GetName()  << "'" << std::endl;
        }
        destination->SaveSelf(kTRUE);
    }

private:
    std::shared_ptr<TFile> originalFile, referenceFile, outputFile;
};
