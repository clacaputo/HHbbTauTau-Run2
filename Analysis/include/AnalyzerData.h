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
        return Get((type*)nullptr, #name, key, ##__VA_ARGS__); \
    } \
    root_ext::SmartHistogram< type >& name() { return name(""); } \
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
        cd();
        for(const auto& iter : data) {
            iter.second->WriteRootObject();
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

    template<typename ValueType, typename KeySuffix, typename ...Args>
    SmartHistogram<ValueType>& Get(ValueType*, const std::string& name, const KeySuffix& suffix, Args... args)
    {
        std::ostringstream ss_suffix;
        ss_suffix << suffix;
        std::string full_name =  name;
        if(ss_suffix.str().size())
            full_name += "_" + ss_suffix.str();
        if(!data.count(full_name)) {
            AbstractHistogram* h = HistogramFactory<ValueType>::Make(full_name, args...);
            data[full_name] = h;
        }
        return *static_cast< SmartHistogram<ValueType>* >(data[full_name]);
    }

    template<typename ValueType>
    SmartHistogram<ValueType>& Get(ValueType* null_value, const std::string& name)
    {
        return Get(null_value, name, "");
    }


    template<typename ValueType>
    SmartHistogram<ValueType>& Get(const std::string& name)
    {
        return Get((ValueType*)nullptr, name, "");
    }

private:
    typedef std::map< std::string, AbstractHistogram* > DataMap;
    DataMap data;

    TFile* outputFile;
    bool ownOutputFile;
    std::string directoryName;
};
} // root_ext
