/*!
 * \file SmartHistogram.h
 * \brief Definition of class SmartHistogram that allows to create ROOT-compatible histograms.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-20 created
 */

#pragma once

#include <deque>
#include <string>
#include <limits>
#include <stdexcept>

#include <TH1D.h>

namespace root_ext {

template<typename ValueType, typename RootHistogram>
class BaseSmartHistogram {
public:
    typedef std::deque<ValueType>::const_iterator const_iterator;

    BaseSmartHistogram()
        : limits(std::numeric_limits<ValueType>::max(), std::numeric_limits<ValueType>::lowest()),
          minMax(std::numeric_limits<ValueType>::max(), std::numeric_limits<ValueType>::lowest()),
          numberOfBins(1), binSize(1), lowLimitIsFixed(false), highLimitIsFixed(false), numberOfBinsIsFixed(false),
          binSizeIsFixed(false), useBinSize(false)
    {}

    const std::string& Name() const { return name; }
    void setName(const std::string& _name) { name = _name; }
    const std::deque<ValueType>& Data() const { return data; }
    const size_t size() const { return data.size(); }
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }

    void setLowLimit(const ValueType& value) { limits.first = value; lowLimitIsFixed = true; }
    const ValueType& LowLimit() const { return limits.first; }
    void setHighLimit(const ValueType& value) { limits.second = value; highLimitIsFixed = true; }
    const ValueType& HighLimit() const { return limits.second; }
    const ValueType& Min() const { return minMax.first; }
    const ValueType& Max() const { return minMax.second; }

    void setNumberOfBins(size_t n)
    {
        if(!n)
            throw std::runtime_error("Number of bins can't be zero.");
        numberOfBins = n;
        numberOfBinsIsFixed = true;
        binSizeIsFixed = false;
        useBinSize = false;
    }
    size_t NumberOfBins() const { return numberOfBins; }
    void setBinSize(const ValueType& value)
    {
        if(!value)
            throw std::runtime_error("Bin size can't be zero.");
        binSize = value;
        binSizeIsFixed = true;
        numberOfBinsIsFixed = false;
        useBinSize = true;
    }
    const ValueType& BinSize() const { return binSize; }

    bool LowLimitIsFixed() const { return lowLimitIsFixed; }
    bool HighLimitIsFixed() const { return highLimitIsFixed; }
    bool NumberOfBinsIsFixed() const { return numberOfBinsIsFixed; }
    bool BinSizeIsFixed() const { return binSizeIsFixed; }
    bool UseBinSize() const { return useBinSize; }

    void ResetFixedParameters()
    {
        lowLimitIsFixed = highLimitIsFixed = numberOfBinsIsFixed = binSizeIsFixed = useBinSize = false;
    }

    void Fill(int value)
    {
        data.push_back(value);
        minMax.first = std::min(minMax.first, value);
        minMax.second = std::max(minMax.second, value);
    }

    virtual void Adjust()
    {
        if(!lowLimitIsFixed)
            limits.first = minMax.first;
        if(!highLimitIsFixed)
            limits.second = minMax.second;

        if(useBinSize)
            numberOfBins = static_cast<size_t>((limits.second - limits.first) / binSize);
        else
            binSize = (limits.second - limits.first) / numberOfBins;
    }

    RootHistogram* ProduceRootHistogram()
    {
        Adjust();
        RootHistogram* rootHistogram = CreateRootHistogram();
        for(const ValueType& value : data)
            rootHistogram->Fill(value);
        return rootHistogram;
    }

private:
    virtual RootHistogram* CreateRootHistogram() const = 0;

protected:
    std::string name;
    std::deque<ValueType> data;
    std::pair<ValueType, ValueType> limits, minMax;
    size_t numberOfBins;
    ValueType binSize;
    bool lowLimitIsFixed, highLimitIsFixed, numberOfBinsIsFixed, binSizeIsFixed, useBinSize;
};

template<typename ValueType, typename RootHistogram>
class SmartHistogram {};

template<>
class SmartHistogram<double, TH1D> : public BaseSmartHistogram<double, TH1D> {
private:

    virtual void Adjust()
    {
        static const size_t defaultNumberOfBins = 100;
        if(!lowLimitIsFixed)
            limits.first = minMax.first;
        if(!highLimitIsFixed)
            limits.second = minMax.second;

        const double interval = limits.second - limits.first;
        if(useBinSize) {
            if(!binSizeIsFixed)
                binSize = interval / defaultNumberOfBins;
            numberOfBins = static_cast<size_t>(interval / binSize);
        }
        else {
            if(!numberOfBinsIsFixed)
                numberOfBins = defaultNumberOfBins;
            binSize = interval / numberOfBins;
        }
    }

    virtual RootHistogram* CreateRootHistogram() const
    {
        return new TH1D(Name().c_str(), Name().c_str(), NumberOfBins(), LowLimit(), HighLimit());
    }
};

} // namespace root_ext
