/*!
 * \file TreeExtractor.h
 * \brief Definition of TreeExtractor class that extracts ntuple information.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-28 created
 */

#pragma once

#include <limits>
#include <iostream>
#include <queue>
#include <fstream>
#include "EventDescriptor.h"

namespace analysis {

namespace detail {

const std::vector<std::string> treeNames = {"events", "electrons", "muons", "taus", "PFCand", "jets",  "vertices", "genParticles",
                                           "triggers", "triggerObjects", "METs", /*"METsMVAeTau", "METsMVAmuTau",
                                            "METsMVAtauTau",*/ "METsPF", "METsTC",
                                           "genMETs"};

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
//                    std::shared_ptr<ntuple::METTree>,
//                    std::shared_ptr<ntuple::METTree>,
//                    std::shared_ptr<ntuple::METTree>,
                    std::shared_ptr<ntuple::GenMETTree> > Forest;

template<typename Tree>
inline void CreateTree(std::shared_ptr<Tree>& tree, TFile& inputFile, const std::string& treeName ,bool extractMCtruth)
{
    tree = Tree::IsMCtruth() && !extractMCtruth ? std::shared_ptr<Tree>()
                                                : std::shared_ptr<Tree>( new Tree(inputFile, treeName) );
}

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<Forest>::value >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth) {}

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<Forest>::value) >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth)
{
    CreateTree(std::get<N>(forest), inputFile, treeNames.at(N), extractMCtruth);
    CreateForest<N + 1>(forest, inputFile, extractMCtruth);
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
    TreeExtractor(const std::string& prefix, const std::string& input, bool _extractMCtruth)
        :  extractMCtruth(_extractMCtruth)

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
    std::shared_ptr<TFile> inputFile;
    std::queue<std::string> inputFileNames;
    std::shared_ptr<detail::Forest> forest;
    Long64_t current_entry;
    std::string prefix;
    bool extractMCtruth;

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
        detail::CreateForest(*forest, *inputFile, extractMCtruth);
        return true;
    }
};

} // analysis
