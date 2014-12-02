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

#include <map>
#include <stdexcept>
#include <sstream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "SmartHistogram.h"

#define ENTRY_1D(type, name, ...) \
    template<typename Key> \
    root_ext::SmartHistogram< type >& name(const Key& key) { \
        return Get((type*)nullptr, #name, key, ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< type >& name() { return name(""); } \
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

#define TH1D_ENTRY_EX(name, nbinsx, xlow, xup, x_axis_title, y_axis_title, use_log_y, max_y_sf) \
    ENTRY_1D(TH1D, name, nbinsx, xlow, xup, x_axis_title, y_axis_title, use_log_y, max_y_sf)
#define TH1D_ENTRY_FIX_EX(name, binsizex, nbinsx, xlow, x_axis_title, y_axis_title, use_log_y, max_y_sf) \
    TH1D_ENTRY_EX(name, nbinsx, xlow, (xlow+binsizex*nbinsx), title, x_axis_title, y_axis_title, use_log_y, max_y_sf)
#define TH1D_ENTRY_CUSTOM_EX(name, bins, x_axis_title, y_axis_title, use_log_y, max_y_sf) \
    ENTRY_1D(TH1D, name, bins, x_axis_title, y_axis_title, use_log_y, max_y_sf)

#define TH2D_ENTRY(name, nbinsx, xlow, xup, nbinsy, ylow, yup) \
    ENTRY_1D(TH2D, name, nbinsx, xlow, xup, nbinsy, ylow, yup)
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

public:
    static const std::set<std::string>& GetAllHistogramNames() { return HistogramNames(); }

public:
    AnalyzerData() : outputFile(nullptr), ownOutputFile(false) {}
    AnalyzerData(const std::string& outputFileName)
        : outputFile(new TFile(outputFileName.c_str(), "RECREATE")), ownOutputFile(true)
    {
        if(outputFile->IsZombie())
            throw std::runtime_error("Output file not created.");
        outputFile->cd();
    }
    AnalyzerData(TFile& _outputFile, const std::string& _directoryName = "")
        : outputFile(&_outputFile), ownOutputFile(false), directoryName(_directoryName)
    {
        outputFile->cd();
        if (directoryName.size()){
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
        }
    }

    virtual ~AnalyzerData()
    {
        if(outputFile) {
            cd();
            for(const auto& iter : data) {
                iter.second->WriteRootObject();
                delete iter.second;
            }
            if(ownOutputFile)
                delete outputFile;
        }
    }

    TFile& getOutputFile(){return *outputFile;}
    void cd() const
    {
        if(outputFile)
            outputFile->cd();
        if(directoryName.size())
            outputFile->cd(directoryName.c_str());
    }

    bool Contains(const std::string& name) const
    {
        return data.find(name) != data.end();
    }

    void Erase(const std::string& name)
    {
        auto iter = data.find(name);
        if(iter != data.end()) {
            delete iter->second;
            data.erase(iter);
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
    SmartHistogram<ValueType>& Get(const ValueType*, const std::string& name, const KeySuffix& suffix, Args... args)
    {
        std::ostringstream ss_suffix;
        ss_suffix << suffix;
        std::string full_name =  name;
        if(ss_suffix.str().size())
            full_name += "_" + ss_suffix.str();
        if(!data.count(full_name)) {
            cd();
            AbstractHistogram* h = HistogramFactory<ValueType>::Make(full_name, args...);
            data[full_name] = h;
            HistogramNames().insert(h->Name());
            if(!outputFile)
                h->DetachFromFile();
        }
        return *static_cast< SmartHistogram<ValueType>* >(data[full_name]);
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
        AbstractHistogram* h = new SmartHistogram<ValueType>(original);
        data[h->Name()] = h;
        HistogramNames().insert(h->Name());
        if(!outputFile)
            h->DetachFromFile();
        return *static_cast< SmartHistogram<ValueType>* >(h);
    }

private:
    typedef std::map< std::string, AbstractHistogram* > DataMap;
    DataMap data;

    TFile* outputFile;
    bool ownOutputFile;
    std::string directoryName;
};
} // root_ext
