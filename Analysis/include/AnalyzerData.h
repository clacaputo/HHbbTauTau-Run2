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

#define DATA_ENTRY(type, name, ...) \
    type&  name() { \
        if(!Contains(#name)) \
            data[#name] = new type(__VA_ARGS__); \
        return Get< type >(#name); } \
    /**/

#define MAP_DATA_ENTRY(type, name, ...) \
    template<typename Key> \
    type& name(const Key& key) { \
        std::ostringstream ss; \
        ss << #name << "_" << key; \
        const std::string& full_name = ss.str(); \
        if(!Contains(full_name)) { \
            type* object = new type(__VA_ARGS__); \
            TObject* clone = object->Clone(full_name.c_str()); \
            delete object; \
            data[full_name] = clone; \
        } \
        return Get< type >(full_name); \
    } \
    /**/

#define TH1D_ENTRY(name, nbinsx, xlow, xup) DATA_ENTRY(TH1D, name, #name, #name, nbinsx, xlow, xup)
#define TH1D_ENTRY_FIX(name, binsizex, nbinsx, xlow) TH1D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx))
#define TH1D_MAP_ENTRY(name, nbinsx, xlow, xup) MAP_DATA_ENTRY(TH1D,name, #name, #name, nbinsx, xlow, xup)
#define TH1D_MAP_ENTRY_FIX(name, binsizex, nbinsx, xlow) TH1D_MAP_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx))

#define TH2D_ENTRY(name, nbinsx, xlow, xup, nbinsy, ylow, yup) DATA_ENTRY(TH2D, name, #name, #name, nbinsx, xlow, xup,\
        nbinsy, ylow, yup)
#define TH2D_ENTRY_FIX(name, binsizex, nbinsx, xlow, binsizey, nbinsy, ylow) \
    TH2D_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx), nbinsy, ylow, (ylow+binsizey*nbinsy))
#define TH2D_MAP_ENTRY(name, nbinsx, xlow, xup, nbinsy, ylow, yup) \
    MAP_DATA_ENTRY(TH2D, name, #name, #name, nbinsx, xlow, xup, nbinsy, ylow, yup)
#define TH2D_MAP_ENTRY_FIX(name, binsizex, nbinsx, xlow, binsizey, nbinsy, ylow) \
    TH2D_MAP_ENTRY(name, nbinsx, xlow, (xlow+binsizex*nbinsx), nbinsy, ylow, (ylow+binsizey*nbinsy))


namespace root_ext {
class AnalyzerData {
protected:
    typedef std::map<std::string, TObject*> DataMap;
    typedef std::map< std::string, AbstractHistogram* > SmartDataMap;
    DataMap data;
    SmartDataMap smartData;

    bool Contains(const std::string& name) const
    {
        return data.find(name) != data.end();
    }

    template<typename T>
    T& Get(const std::string& name)
    {
        return *static_cast<T*>(data[name]);
    }

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
        for(DataMap::iterator iter = data.begin(); iter != data.end(); ++iter) {
            iter->second->Write();
            delete iter->second;
        }
        for(const auto& iter : smartData) {
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

    template<typename ValueType>
    SmartHistogram<ValueType, TH1D>& getSmart(const std::string& name)
    {
        if(!smartData.count(name)) {
            AbstractHistogram* h = HistogramFactory<ValueType>::Make();
            h->setName(name);
            smartData[name] = h;
        }

        return *static_cast< SmartHistogram<ValueType, TH1D>* >(smartData[name]);
    }

private:
    TFile* outputFile;
    bool ownOutputFile;
    std::string directoryName;
};
} // root_ext
