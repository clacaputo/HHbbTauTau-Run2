/*!
 * \class AnalyzerData AnalyzerData.h
 * \brief Base class for Analyzer data containers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2013-03-29 created
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
        std::ostringstream ss_key; \
        ss_key << key; \
        type* t = nullptr; \
        return Get(t, #name, ss_key.str(), ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< type >& name() { \
        return name<std::string>(""); \
    } \
    /**/

#define ENTRY_2D(type, name, ...) \
    template<typename Key> \
    root_ext::SmartHistogram< std::pair<type, type> >& name(const Key& key) { \
        std::ostringstream ss_key; \
        ss_key << key; \
        std::pair<type, type>* t = nullptr; \
        return Get(t, #name, ss_key.str(), ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< std::pair<type, type> >& name() { \
        return name<std::string>(""); \
    } \
    /**/

#define TH1D_ENTRY(name, nbinsx, xlow, xup) ENTRY_1D(TH1D, name, nbinsx, xlow, xup)
#define TH1D_ENTRY_FIX(name, binsizex, nbinsx, xlow) TH1D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx))

#define TH2D_ENTRY(name, nbinsx, xlow, xup, nbinsy, ylow, yup) \
    ENTRY_1D(TH2D, name, nbinsx, xlow, xup, nbinsy, ylow, yup)
#define TH2D_ENTRY_FIX(name, binsizex, nbinsx, xlow, binsizey, nbinsy, ylow) \
    TH2D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx), nbinsy, ylow, (ylow+binsizey*nbinsy))

namespace root_ext {
class AnalyzerData {
public:
    AnalyzerData() : outputFile(nullptr), ownOutputFile(false) {}
    AnalyzerData(const std::string& outputFileName)
        : outputFile(new TFile(outputFileName.c_str(), "RECREATE")), ownOutputFile(true)
    {
        if(outputFile->IsZombie())
            throw std::runtime_error("Output file not created.");
        outputFile->cd();
    }
    AnalyzerData(TFile& _outputFile, const std::string& _directoryName)
        : outputFile(&_outputFile), ownOutputFile(false), directoryName(_directoryName)
    {
        outputFile->cd();
        outputFile->mkdir(directoryName.c_str());
        outputFile->cd(directoryName.c_str());
    }

    virtual ~AnalyzerData()
    {
        cd();
        for(const auto& iter : data) {
            iter.second->WriteRootHistogram();
            delete iter.second;
        }
        if(ownOutputFile)
            delete outputFile;
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

    template<typename ValueType, typename ...Args>
    SmartHistogram<ValueType>& Get(ValueType*, const std::string& name, const std::string& suffix, Args... args)
    {
        std::string full_name =  name;
        if(suffix.size())
            full_name += "_" + suffix;
        if(!data.count(full_name)) {
            AbstractHistogram* h = HistogramFactory<ValueType>::Make(full_name, args...);
            data[full_name] = h;
        }
        return *static_cast< SmartHistogram<ValueType>* >(data[full_name]);
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Get(const std::string& name)
    {
        ValueType *t = nullptr;
        return Get(t, name, "");
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Get(const std::string& name, const std::string& suffix)
    {
        ValueType *t = nullptr;
        return Get(t, name, suffix);
    }

private:
    typedef std::map< std::string, AbstractHistogram* > DataMap;
    DataMap data;

    TFile* outputFile;
    bool ownOutputFile;
    std::string directoryName;
};
} // root_ext
