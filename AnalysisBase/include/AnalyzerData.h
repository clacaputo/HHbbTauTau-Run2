/*!
 * \class AnalyzerData AnalyzerData.h
 * \brief Base class for Analyzer data containers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2013-03-29 created
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

#include <vector>
#include <map>
#include <stdexcept>
#include <sstream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Compression.h>

#include "SmartHistogram.h"

#define ENTRY_1D(type, name, ...) \
    template<typename Key> \
    root_ext::SmartHistogram< type >& name(const Key& key) { \
        return Get((type*)nullptr, #name, key, ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< type >& name() { \
        static const size_t index = GetUniqueIndex(#name); \
        return GetFast((type*)nullptr, #name, index, ##__VA_ARGS__); \
    } \
    static std::string name##_Name() { return #name; } \
    /**/

#define ENTRY_2D(type, name, ...) \
    template<typename Key> \
    root_ext::SmartHistogram< root_ext::detail::Base2DHistogram<type>::Value >& name(const Key& key) { \
        return Get((root_ext::detail::Base2DHistogram<type>::Value*)nullptr, #name, key, ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< root_ext::detail::Base2DHistogram<type>::Value >& name() { return name(""); } \
    /**/

#define TH1D_ENTRY(name, nbinsx, xlow, xup) ENTRY_1D(TH1D, name, nbinsx, xlow, xup)
#define TH1D_ENTRY_FIX(name, binsizex, nbinsx, xlow) TH1D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx))
#define TH1D_ENTRY_CUSTOM(name, bins) ENTRY_1D(TH1D, name, bins)

#define TH1D_ENTRY_EX(name, nbinsx, xlow, xup, x_axis_title, y_axis_title, use_log_y, max_y_sf, store) \
    ENTRY_1D(TH1D, name, nbinsx, xlow, xup, x_axis_title, y_axis_title, use_log_y, max_y_sf, store)
#define TH1D_ENTRY_FIX_EX(name, binsizex, nbinsx, xlow, x_axis_title, y_axis_title, use_log_y, max_y_sf, store) \
    TH1D_ENTRY_EX(name, nbinsx, xlow, (xlow+binsizex*nbinsx), title, x_axis_title, y_axis_title, use_log_y, \
                  max_y_sf, store)
#define TH1D_ENTRY_CUSTOM_EX(name, bins, x_axis_title, y_axis_title, use_log_y, max_y_sf, store) \
    ENTRY_1D(TH1D, name, bins, x_axis_title, y_axis_title, use_log_y, max_y_sf, store)

#define TH2D_ENTRY(name, nbinsx, xlow, xup, nbinsy, ylow, yup) \
    ENTRY_2D(TH2D, name, nbinsx, xlow, xup, nbinsy, ylow, yup)
#define TH2D_ENTRY_FIX(name, binsizex, nbinsx, xlow, binsizey, nbinsy, ylow) \
    TH2D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx), nbinsy, ylow, (ylow+binsizey*nbinsy))

namespace root_ext {
class AnalyzerData {
private:
    static std::set<std::string>& HistogramNames()
    {
        static std::set<std::string> names;
        return names;
    }

    static std::set<std::string>& OriginalHistogramNames()
    {
        static std::set<std::string> names;
        return names;
    }

    static std::map<std::string, size_t>& IndexMap()
    {
        static std::map<std::string, size_t> index_map;
        return index_map;
    }

    static constexpr size_t MaxIndex = 1000;

public:
    static const std::set<std::string>& GetAllHistogramNames() { return HistogramNames(); }
    static const std::set<std::string>& GetOriginalHistogramNames() { return OriginalHistogramNames(); }

    static size_t GetUniqueIndex(const std::string& name)
    {
        const auto iter = IndexMap().find(name);
        if(iter != IndexMap().end())
            return iter->second;
        const size_t index = IndexMap().size();
        IndexMap()[name] = index;
        return index;
    }

    static TFile* CreateFile(const std::string& fileName)
    {
        TFile* file = new TFile(fileName.c_str(), "RECREATE", "", ROOT::kZLIB * 100 + 9);
        if(file->IsZombie()) {
            std::ostringstream ss;
            ss << "File '" << fileName << "' not created.";
            throw std::runtime_error(ss.str());
        }
        return file;
    }

public:
    AnalyzerData() : data_vector(MaxIndex), outputFile(nullptr), ownOutputFile(false), directory(nullptr)
    {
        data_vector.assign(MaxIndex, nullptr);
    }

    AnalyzerData(const std::string& outputFileName)
        : outputFile(CreateFile(outputFileName)), ownOutputFile(true)
    {
        data_vector.assign(MaxIndex, nullptr);
        directory = outputFile;
        cd();
    }

    AnalyzerData(TFile& _outputFile, const std::string& directoryName = "")
        : outputFile(&_outputFile), ownOutputFile(false), directory(nullptr)
    {
        data_vector.assign(MaxIndex, nullptr);
        if (directoryName.size()){
            outputFile->mkdir(directoryName.c_str());
            directory = outputFile->GetDirectory(directoryName.c_str());
            if(!directory)
                throw std::runtime_error("Unable to create analyzer data directory.");
        } else
            directory = outputFile;
        cd();
    }

    virtual ~AnalyzerData()
    {
        cd();
        for(const auto& iter : data) {
            if(directory)
                iter.second->WriteRootObject();
            delete iter.second;
        }
        if(ownOutputFile)
            delete outputFile;
    }

    TFile& getOutputFile() { return *outputFile; }

    void cd() const
    {
        if(directory && !directory->cd())
            throw std::runtime_error("Unable to cd to analyzer data directory.");
    }

    bool Contains(const std::string& name) const { return data.find(name) != data.end(); }

    void Erase(const std::string& name)
    {
        auto iter = data.find(name);
        if(iter != data.end()) {
            delete iter->second;
            data.erase(iter);
            auto index_iter = IndexMap().find(name);
            if(index_iter != IndexMap().end() && index_iter->second < MaxIndex)
                data_vector.at(index_iter->second) = nullptr;
        }
    }

    std::vector<std::string> KeysCollection() const
    {
        std::vector<std::string> keys;
        for(const auto& iter : data)
            keys.push_back(iter.first);
        return keys;
    }

    template<typename ValueType, typename KeySuffix, typename ...Args>
    SmartHistogram<ValueType>& Get(const ValueType* ptr, const std::string& name, const KeySuffix& suffix, Args... args)
    {

        std::ostringstream ss_suffix;
        ss_suffix << suffix;
        const std::string s_suffix = ss_suffix.str();
        const std::string full_name = s_suffix.size() ? name + "_" + s_suffix : name;
        return GetByFullName(ptr, name, full_name, args...);
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Get(const ValueType* null_value, const std::string& name)
    {
        return Get(null_value, name, "");
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Get(const std::string& name)
    {
        return Get((ValueType*)nullptr, name, "");
    }

    template<typename ValueType>
    SmartHistogram<ValueType>* GetPtr(const std::string& name)
    {
        if(!Contains(name)) return nullptr;
        return &Get<ValueType>(name);
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Clone(const SmartHistogram<ValueType>& original)
    {
        if(data.count(original.Name()))
            throw std::runtime_error("histogram already exists");
        cd();
        SmartHistogram<ValueType>* h = new SmartHistogram<ValueType>(original);
        data[h->Name()] = h;
        HistogramNames().insert(h->Name());
        h->SetOutputDirectory(directory);
        auto index_iter = IndexMap().find(h->Name());
        if(index_iter != IndexMap().end() && index_iter->second < MaxIndex)
            data_vector.at(index_iter->second) = h;
        return *h;
    }

protected:
    template<typename ValueType, typename ...Args>
    SmartHistogram<ValueType>& GetFast(const ValueType* ptr, const std::string& name, size_t index, Args... args)
    {
        if(index < MaxIndex && data_vector[index] != nullptr)
            return *static_cast< SmartHistogram<ValueType>* >(data_vector[index]);
        return GetByFullName(ptr, name, name, args...);
    }

private:
    template<typename ValueType, typename ...Args>
    SmartHistogram<ValueType>& GetByFullName(const ValueType*, const std::string& name, const std::string& full_name,
                                             Args... args)
    {
        auto iter = data.find(full_name);
        if(iter == data.end()) {
            cd();
            AbstractHistogram* h = HistogramFactory<ValueType>::Make(full_name, args...);
            data[full_name] = h;
            HistogramNames().insert(h->Name());
            OriginalHistogramNames().insert(name);
            h->SetOutputDirectory(directory);
            iter = data.find(full_name);
            auto index_iter = IndexMap().find(full_name);
            if(index_iter != IndexMap().end() && index_iter->second < MaxIndex)
                data_vector.at(index_iter->second) = h;
        }
        SmartHistogram<ValueType>* result = dynamic_cast< SmartHistogram<ValueType>* >(iter->second);
        if(!result) {
            std::ostringstream ss;
            ss << "Wrong type for histogram '" << full_name << "'.";
            throw std::runtime_error(ss.str());
        }
        return *result;
    }

private:
    typedef std::vector<AbstractHistogram*> DataVector;
    typedef std::map<std::string, AbstractHistogram*> DataMap;
    DataMap data;
    DataVector data_vector;


    TFile* outputFile;
    bool ownOutputFile;
    TDirectory* directory;
    std::string directoryName;
};
} // root_ext
