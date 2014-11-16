/*!
 * \file SmartTree.h
 * \brief Definition of SmartTree class.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 *
 * Copyright 2013, 2014 Konstantin Androsov <konstantin.androsov@gmail.com>
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

#include <stdexcept>
#include <sstream>
#include <memory>
#include <iostream>

#ifdef SMART_TREE_FOR_CMSSW
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#else
#include "TFile.h"
#endif

#include "TTree.h"

#define SIMPLE_TREE_BRANCH(type, name) \
private: type _##name; \
public:  type& name() { return _##name; }

#define VECTOR_TREE_BRANCH(type, name) \
private: std::vector< type > _##name; \
public:  std::vector< type >& name() { return _##name; }

#define SIMPLE_DATA_TREE_BRANCH(type, name) \
    type& name() { return data.name; }

#define VECTOR_DATA_TREE_BRANCH(type, name) \
    std::vector< type >& name() { return data.name; }

#define DECLARE_SIMPLE_BRANCH_VARIABLE(type, name) type name;
#define DECLARE_VECTOR_BRANCH_VARIABLE(type, name) std::vector< type > name;

#define ADD_SIMPLE_TREE_BRANCH(name) AddSimpleBranch(#name, _##name);
#define ADD_SIMPLE_DATA_TREE_BRANCH(name) AddSimpleBranch(#name, data.name);
#define ADD_VECTOR_TREE_BRANCH(name) AddVectorBranch(#name, _##name);
#define ADD_VECTOR_DATA_TREE_BRANCH(name) AddVectorBranch(#name, data.name);

#define DATA_CLASS(namespace_name, class_name, data_macro) \
    namespace namespace_name { \
        struct class_name : public root_ext::detail::BaseDataClass { data_macro() }; \
        typedef std::vector< class_name > class_name##Vector; \
    } \
    /**/

#define TREE_CLASS(namespace_name, tree_class_name, data_macro, data_class_name, tree_name, is_mc_truth) \
    namespace namespace_name { \
    class tree_class_name : public root_ext::SmartTree { \
    public: \
        static bool IsMCtruth() { return is_mc_truth; } \
        static const std::string& Name() { static const std::string name = tree_name; return name; } \
        tree_class_name() : SmartTree(Name(), "/") { Initialize(); } \
        tree_class_name(const std::string& name) : SmartTree(name, "/") { Initialize(); } \
        tree_class_name(TFile& file) : SmartTree(Name(), file) { Initialize(); } \
        tree_class_name(TFile& file, const std::string& name) : SmartTree(name, file) { Initialize(); } \
        data_class_name data; \
        data_macro() \
    private: \
        inline void Initialize(); \
    }; \
    } \
    /**/

#define TREE_CLASS_WITH_EVENT_ID(namespace_name, tree_class_name, data_macro, data_class_name, tree_name, is_mc_truth) \
    namespace namespace_name { \
    class tree_class_name : public root_ext::SmartTree { \
    public: \
        static bool IsMCtruth() { return is_mc_truth; } \
        static const std::string& Name() { static const std::string name = tree_name; return name; } \
        tree_class_name() : SmartTree(Name(), "/") { Initialize(); } \
        tree_class_name(const std::string& name) : SmartTree(name, "/") { Initialize(); } \
        tree_class_name(TFile& file) : SmartTree(Name(), file) { Initialize(); } \
        tree_class_name(TFile& file, const std::string& name) : SmartTree(name, file) { Initialize(); } \
        data_class_name data; \
        SIMPLE_TREE_BRANCH(UInt_t, RunId) \
        SIMPLE_TREE_BRANCH(UInt_t, LumiBlock) \
        SIMPLE_TREE_BRANCH(UInt_t, EventId) \
        data_macro() \
    private: \
        inline void Initialize(); \
    }; \
    } \
    /**/

#define TREE_CLASS_INITIALIZE(namespace_name, tree_class_name, data_macro) \
    namespace namespace_name { \
        inline void tree_class_name::Initialize() { \
            data_macro() \
            if (GetEntries() > 0) GetEntry(0); \
        } \
    } \
    /**/

#define TREE_CLASS_WITH_EVENT_ID_INITIALIZE(namespace_name, tree_class_name, data_macro) \
    namespace namespace_name { \
        inline void tree_class_name::Initialize() { \
            ADD_SIMPLE_TREE_BRANCH(RunId) \
            ADD_SIMPLE_TREE_BRANCH(LumiBlock) \
            ADD_SIMPLE_TREE_BRANCH(EventId) \
            data_macro() \
            if (GetEntries() > 0) GetEntry(0); \
        } \
    } \
    /**/

namespace root_ext {
namespace detail {
    struct BaseSmartTreeEntry {
        virtual ~BaseSmartTreeEntry() {}
        virtual void clear() = 0;
    };

    template<typename DataType>
    struct SmartTreeVectorPtrEntry : public BaseSmartTreeEntry {
        std::vector<DataType>* value;
        SmartTreeVectorPtrEntry(std::vector<DataType>& origin)
            : value(&origin) {}
        virtual void clear() { value->clear(); }
    };

    struct BaseDataClass {
        virtual ~BaseDataClass() {}
    };
} // detail

class SmartTree {
public:
    explicit SmartTree(const std::string& _name, const std::string& _directory = "", Long64_t maxVirtualSize = 10000000)
        : name(_name), file(nullptr), readMode(false), directory(_directory)
    {
#ifdef SMART_TREE_FOR_CMSSW
        edm::Service<TFileService> tfserv;
        if(directory.size()) {
            tfserv->file().cd(directory.c_str());
			tree = new TTree(name.c_str(), name.c_str());
		}
		else
			tree = tfserv->make<TTree>(name.c_str(), name.c_str());
#else
        tree = new TTree(name.c_str(), name.c_str());
#endif
        if(maxVirtualSize < 0)
            tree->SetDirectory(0); // detach tree from directory making it "in memory".
        else
            tree->SetMaxVirtualSize(maxVirtualSize);
    }
    SmartTree(const std::string& _name, TFile& _file)
        : name(_name), file(&_file), readMode(true)
    {
        tree = static_cast<TTree*>(file->Get(name.c_str()));
        if(!tree)
            throw std::runtime_error("Tree not found.");
        if(tree->GetNbranches())
            tree->SetBranchStatus("*", 0);
    }
    SmartTree(const SmartTree& other)
    {
        throw std::runtime_error("Can't copy a smart tree");
    }

    SmartTree(const SmartTree&& other)
    {
        entries = other.entries;
        readMode = other.readMode;
        tree = other.tree;
        directory = other.directory;
    }

    virtual ~SmartTree()
    {
        if(readMode) file->Delete(name.c_str());
        else delete tree;
    }

    void Fill()
    {
        tree->Fill();
        for(auto& entry : entries)
            entry.second->clear();
    }

    Long64_t GetEntries() const { return tree->GetEntries(); }
    Long64_t GetReadEntry() const { return tree->GetReadEntry(); }
    Int_t GetEntry(Long64_t entry) { return tree->GetEntry(entry); }
    void Write()
    {
#ifdef SMART_TREE_FOR_CMSSW
        edm::Service<TFileService> tfserv;
        if(directory.size())
            tfserv->file().cd(directory.c_str());
#endif
        tree->Write("", TObject::kWriteDelete);
    }

protected:
    template<typename DataType>
    void AddSimpleBranch(const std::string& name, DataType& value)
    {
        if(readMode) {
            try {
                EnableBranch(name);
                tree->SetBranchAddress(name.c_str(), &value);
                if(tree->GetReadEntry() >= 0)
                    tree->GetBranch(name.c_str())->GetEntry(tree->GetReadEntry());
            } catch(std::runtime_error& error) {
                std::cerr << "ERROR: " << error.what() << std::endl;
            }
        } else {
            TBranch* branch = tree->Branch(name.c_str(), &value);
            const Long64_t n_entries = tree->GetEntries();
            for(Long64_t n = 0; n < n_entries; ++n)
                branch->Fill();
        }
    }

    template<typename DataType>
    void AddVectorBranch(const std::string& name, std::vector<DataType>& value)
    {
        typedef detail::SmartTreeVectorPtrEntry<DataType> PtrEntry;
        auto entry = std::shared_ptr<PtrEntry>( new PtrEntry(value) );
        if(entries.count(name))
            throw std::runtime_error("Entry is already defined.");
        entries[name] = entry;
        if(readMode) {
            try {
                EnableBranch(name);
                tree->SetBranchAddress(name.c_str(), &entry->value);
                if(tree->GetReadEntry() >= 0)
                    tree->GetBranch(name.c_str())->GetEntry(tree->GetReadEntry());
            } catch(std::runtime_error& error) {
                std::cerr << "ERROR: " << error.what() << std::endl;
            }
        } else {
            TBranch* branch = tree->Branch(name.c_str(), entry->value);
            const Long64_t n_entries = tree->GetEntries();
            for(Long64_t n = 0; n < n_entries; ++n)
                branch->Fill();
        }
    }

    void EnableBranch(const std::string& name)
    {
        UInt_t n_found = 0;
        tree->SetBranchStatus(name.c_str(), 1, &n_found);
        if(n_found != 1) {
            std::ostringstream ss;
            ss << "Branch '" << name << "' is not found.";
            throw std::runtime_error(ss.str());
        }
    }

private:
    std::string name;
    TFile* file;
    std::map< std::string, std::shared_ptr<detail::BaseSmartTreeEntry> > entries;
    bool readMode;
    TTree* tree;
    std::string directory;
};

} // root_ext
