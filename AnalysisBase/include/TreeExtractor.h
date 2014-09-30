/*!
 * \file TreeExtractor.h
 * \brief Definition of TreeExtractor class that extracts ntuple information.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-28 created
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

#include <limits>
#include <iostream>
#include <queue>
#include <fstream>

#include "EventDescriptor.h"

namespace analysis {

namespace detail {

template<typename Tree>
struct TreeVersionTag {
    static unsigned GetVersion() { return 1; }
};

template<>
struct TreeVersionTag<ntuple::GenEventTree> {
    static unsigned GetVersion() { return 2; }
};

const std::vector<std::string> treeNames = { "events", "electrons", "muons", "taus", "PFCand", "jets", "vertices",
                                             "genParticles", "triggers", "triggerObjects", "METs", "METsPF", "METsTC",
                                             "genMETs", "genEvents"
};

typedef std::tuple< std::shared_ptr<ntuple::EventTree>,
                    std::shared_ptr<ntuple::ElectronTree>,
                    std::shared_ptr<ntuple::MuonTree>,
                    std::shared_ptr<ntuple::TauTree>,
                    std::shared_ptr<ntuple::PFCandTree>,
                    std::shared_ptr<ntuple::JetTree>,
                    std::shared_ptr<ntuple::VertexTree>,
                    std::shared_ptr<ntuple::GenParticleTree>,
                    std::shared_ptr<ntuple::TriggerTree>,
                    std::shared_ptr<ntuple::TriggerObjectTree>,
                    std::shared_ptr<ntuple::METTree>,
                    std::shared_ptr<ntuple::METTree>,
                    std::shared_ptr<ntuple::METTree>,
                    std::shared_ptr<ntuple::GenMETTree>,
                    std::shared_ptr<ntuple::GenEventTree> > Forest;

template<typename Tree>
inline void CreateTree(std::shared_ptr<Tree>& tree, TFile& inputFile,
                       const std::string& treeName, bool extractMCtruth, unsigned maxVersion)
{
    const bool createTree = TreeVersionTag<Tree>::GetVersion() <= maxVersion && (!Tree::IsMCtruth() || extractMCtruth);
    tree = createTree ? std::shared_ptr<Tree>( new Tree(inputFile, treeName) ) : std::shared_ptr<Tree>();
}

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<Forest>::value >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth, unsigned maxVersion) {}

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<Forest>::value) >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth, unsigned maxVersion)
{
    CreateTree(std::get<N>(forest), inputFile, treeNames.at(N), extractMCtruth, maxVersion);
    CreateForest<N + 1>(forest, inputFile, extractMCtruth, maxVersion);
}

template<typename Tree, typename ObjectType>
void ReadTree(std::shared_ptr<Tree>& tree, ObjectType& container, Long64_t& current_entry, EventId& currentEventId)
{
    if (!tree || tree->GetEntries() <= 0 || current_entry + 1 >= tree->GetEntries()) return;

    if(currentEventId == EventId::Undef_event())
        ++current_entry;
    if(tree->GetEntry(current_entry) < 0)
        throw std::runtime_error("An I/O error while reading tree.");
    const EventId treeEventId(tree->data.run,tree->data.lumis,tree->data.EventId);
    if(currentEventId != EventId::Undef_event() && treeEventId != currentEventId)
        throw std::runtime_error("Inconsistent tree structure.");
    currentEventId = treeEventId;
    container = tree->data;
}

template<typename Tree, typename ObjectType>
void ReadTree(std::shared_ptr<Tree>& tree, std::vector<ObjectType>& container, Long64_t current_entry,
              EventId currentEventId)
{
    if (!tree || tree->GetEntries() <= 0) return;

    for(Long64_t n = tree->GetReadEntry(); n < tree->GetEntries();) {
        const EventId treeEventId(tree->RunId(),tree->LumiBlock(),tree->EventId());
        if (currentEventId != treeEventId) break;
        container.push_back(tree->data);
        if(tree->GetEntry(++n) < 0)
            throw std::runtime_error("An I/O error while reading tree.");
    }
}

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<Forest>::value, bool >::type
ReadForest(Forest& forest, EventTuple& data, Long64_t& current_entry,
           EventId currentEventId = EventId::Undef_event()) { return false; }

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<Forest>::value), bool >::type
ReadForest(Forest& forest, EventTuple& data, Long64_t& current_entry,
           EventId currentEventId = EventId::Undef_event())
{
    ReadTree(std::get<N>(forest), std::get<N>(data), current_entry, currentEventId);
    if(currentEventId == EventId::Undef_event()) return false;
    ReadForest<N + 1>(forest, data, current_entry, currentEventId);
    return true;
}

} // detail

class TreeExtractor{
public:
    TreeExtractor(const std::string& prefix, const std::string& input, bool _extractMCtruth, unsigned _maxTreeVersion)
        :  extractMCtruth(_extractMCtruth), maxTreeVersion(_maxTreeVersion)
    {
        if (input.find(".root") != std::string::npos)
            inputFileNames.push(input);
        else if (input.find(".txt") != std::string::npos){
            std::ifstream inputStream(input);
            while (inputStream.good()) {
                std::string inputFileName;
                std::getline(inputStream,inputFileName);
                if (inputFileName.size())
                    inputFileNames.push(prefix+inputFileName);
              }
        }
        else throw std::runtime_error("Unrecognized input");
        if (!OpenNextFile())
            throw std::runtime_error("No inputFile found");
    }

    bool ExtractNext(EventDescriptor& descriptor)
    {
        descriptor.Clear();
        do {
            if (detail::ReadForest(*forest, descriptor.data(), current_entry))
                return true;
        } while (OpenNextFile());
        return false;
    }

private:
    bool extractMCtruth;
    unsigned maxTreeVersion;
    std::shared_ptr<TFile> inputFile;
    std::queue<std::string> inputFileNames;
    std::shared_ptr<detail::Forest> forest;
    Long64_t current_entry;
    std::string prefix;



    bool OpenNextFile()
    {
        if (inputFileNames.empty()) return false;
        const std::string fileName = inputFileNames.front();
        inputFileNames.pop();
        forest = std::shared_ptr<detail::Forest>(new detail::Forest());
        inputFile = std::shared_ptr<TFile>(new TFile(fileName.c_str(), "READ"));
        if(inputFile->IsZombie()){
            std::ostringstream ss;
            ss << "Input file " << fileName << " not found." ;
            throw std::runtime_error(ss.str());
        }
        std::cout << "File " << fileName << " is opened." << std::endl;
        current_entry = -1;
        detail::CreateForest(*forest, *inputFile, extractMCtruth, maxTreeVersion);
        return true;
    }
};

} // analysis
