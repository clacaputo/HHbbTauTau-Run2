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

#include "EventDescriptor.h"

namespace analysis {

namespace detail {

typedef std::tuple< std::shared_ptr<ntuple::EventTree>,
                    std::shared_ptr<ntuple::ElectronTree>,
                    std::shared_ptr<ntuple::MuonTree>,
                    std::shared_ptr<ntuple::TauTree>,
                    std::shared_ptr<ntuple::JetTree>,
                    std::shared_ptr<ntuple::VertexTree>,
                    std::shared_ptr<ntuple::GenParticleTree> > Forest;

template<typename Tree>
inline void CreateTree(std::shared_ptr<Tree>& tree, TFile& inputFile, bool extractMCtruth)
{
    tree = Tree::IsMCtruth() && !extractMCtruth ? std::shared_ptr<Tree>()
                                                : std::shared_ptr<Tree>( new Tree(inputFile) );
}

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<Forest>::value >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth) {}

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<Forest>::value) >::type
CreateForest(Forest& forest, TFile& inputFile, bool extractMCtruth)
{
    CreateTree(std::get<N>(forest), inputFile, extractMCtruth);
    CreateForest<N + 1>(forest, inputFile, extractMCtruth);
}

template<typename Tree, typename ObjectType>
void ReadTree(std::shared_ptr<Tree>& tree, ObjectType& container, Long64_t& current_entry, UInt_t& currentEventId)
{
    if (!tree || tree->GetEntries() <= 0 || current_entry + 1 >= tree->GetEntries()) return;

    if(currentEventId == std::numeric_limits<UInt_t>::max())
        ++current_entry;
    if(tree->GetEntry(current_entry) < 0)
        throw std::runtime_error("An I/O error while reading tree.");
    if(currentEventId != std::numeric_limits<UInt_t>::max() && tree->data.EventId != currentEventId)
        throw std::runtime_error("Inconsistent tree structure.");
    currentEventId = tree->data.EventId;
    container = tree->data;
}

template<typename Tree, typename ObjectType>
void ReadTree(std::shared_ptr<Tree>& tree, std::vector<ObjectType>& container, Long64_t current_entry,
              UInt_t currentEventId)
{
    if (!tree || tree->GetEntries() <= 0) return;

    for(Long64_t n = tree->GetReadEntry(); n < tree->GetEntries();) {
        if (currentEventId != tree->EventId()) break;
        container.push_back(tree->data);
        if(tree->GetEntry(++n) < 0)
            throw std::runtime_error("An I/O error while reading tree.");
    }
}

template<size_t N = 0>
inline typename std::enable_if< N == std::tuple_size<Forest>::value, bool >::type
ReadForest(Forest& forest, EventTuple& data, Long64_t& current_entry,
           UInt_t currentEventId = std::numeric_limits<UInt_t>::max()) { return false; }

template<size_t N = 0>
inline typename std::enable_if< (N < std::tuple_size<Forest>::value), bool >::type
ReadForest(Forest& forest, EventTuple& data, Long64_t& current_entry,
           UInt_t currentEventId = std::numeric_limits<UInt_t>::max())
{
    ReadTree(std::get<N>(forest), std::get<N>(data), current_entry, currentEventId);
    if(currentEventId == std::numeric_limits<UInt_t>::max()) return false;
    ReadForest<N + 1>(forest, data, current_entry, currentEventId);
    return true;
}

} // detail

class TreeExtractor{
public:
    TreeExtractor(const std::string& inputFileName, bool extractMCtruth)
        : inputFile(new TFile(inputFileName.c_str(), "READ")), forest(new detail::Forest()),
          current_entry(-1)
    {
        if(inputFile->IsZombie())
            throw std::runtime_error("Input file not found.");

        detail::CreateForest(*forest, *inputFile, extractMCtruth);
    }

    bool ExtractNext(EventDescriptor& descriptor)
    {
        descriptor.Clear();
        return detail::ReadForest(*forest, descriptor.data(), current_entry);
    }

private:
    std::shared_ptr<TFile> inputFile;
    std::shared_ptr<detail::Forest> forest;
    Long64_t current_entry;
};

} // analysis
