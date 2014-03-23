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
#include <TH2D.h>

namespace root_ext {

class AbstractHistogram {
public:
    AbstractHistogram(const std::string& _name)
        : name(_name) {}

    virtual ~AbstractHistogram() {}

    virtual void WriteRootHistogram() = 0;

    const std::string& Name() const { return name; }

private:
    std::string name;
};

template<typename ValueType>
class HistogramProperties {
public:
    HistogramProperties()
        : limits(std::numeric_limits<ValueType>::max(), std::numeric_limits<ValueType>::lowest()),
          minMax(std::numeric_limits<ValueType>::max(), std::numeric_limits<ValueType>::lowest()),
          numberOfBins(1), binSize(1), lowLimitIsFixed(false), highLimitIsFixed(false), numberOfBinsIsFixed(false),
          binSizeIsFixed(false), useBinSize(false) {}
    virtual ~HistogramProperties() {}

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

    void ResetFixedParameters()
    {
        lowLimitIsFixed = highLimitIsFixed = numberOfBinsIsFixed = binSizeIsFixed = useBinSize = false;
    }

    bool LowLimitIsFixed() const { return lowLimitIsFixed; }
    bool HighLimitIsFixed() const { return highLimitIsFixed; }
    bool NumberOfBinsIsFixed() const { return numberOfBinsIsFixed; }
    bool BinSizeIsFixed() const { return binSizeIsFixed; }
    bool UseBinSize() const { return useBinSize; }

private:
    virtual void UpdateAfterFill(const ValueType& value)
    {
        minMax.first = std::min(minMax.first, value);
        minMax.second = std::max(minMax.second, value);
    }

protected:
    std::pair<ValueType, ValueType> limits, minMax;
    size_t numberOfBins;
    ValueType binSize;
    bool lowLimitIsFixed, highLimitIsFixed, numberOfBinsIsFixed, binSizeIsFixed, useBinSize;
};

namespace detail {

template<typename ValueType>
class Base1DHistogram : public AbstractHistogram, public HistogramProperties<ValueType> {
public:
    typedef typename std::deque<ValueType>::const_iterator const_iterator;

    Base1DHistogram(const std::string& name) : AbstractHistogram(name) {}
    Base1DHistogram(const std::string& name, size_t nbins, const ValueType& low, const ValueType& high)
        : AbstractHistogram(name)
    {
        this->setNumberOfBins(nbins);
        this->setLowLimit(low);
        this->setHighLimit(high);
    }

    const std::deque<ValueType>& Data() const { return data; }
    const size_t size() const { return data.size(); }
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }

    void Fill(const ValueType& value)
    {
        data.push_back(value);
        this->minMax.first = std::min(this->minMax.first, value);
        this->minMax.second = std::max(this->minMax.second, value);
    }

    TH1D* ProduceRootHistogram()
    {
        Adjust();
        TH1D* rootHistogram = CreateRootHistogram();
        for(const ValueType& value : data)
            rootHistogram->Fill(value);
        return rootHistogram;
    }

    void WriteRootHistogram()
    {
        TH1D* h = ProduceRootHistogram();
        h->Write();
        delete h;
    }

    virtual void Adjust() = 0;

private:
    virtual TH1D* CreateRootHistogram() const = 0;

private:
    std::deque<ValueType> data;
};

template<typename ValueType>
class Base1DFloatingPointHistogram : public Base1DHistogram<ValueType> {
public:
    Base1DFloatingPointHistogram(const std::string& name) : Base1DHistogram<ValueType>(name) {}
    Base1DFloatingPointHistogram(const std::string& name, size_t nbins, const ValueType& low, const ValueType& high) :
        Base1DHistogram<ValueType>(name, nbins, low, high) {}

    virtual void Adjust()
    {
        static const size_t defaultNumberOfBins = 100;
        if(!this->lowLimitIsFixed)
            this->limits.first = this->minMax.first;
        if(!this->highLimitIsFixed)
            this->limits.second = this->minMax.second;

        const ValueType interval = this->limits.second - this->limits.first;
        if(this->useBinSize) {
            if(!this->binSizeIsFixed)
                this->binSize = interval / defaultNumberOfBins;
            this->numberOfBins = static_cast<size_t>(interval / this->binSize);
        }
        else {
            if(!this->numberOfBinsIsFixed)
                this->numberOfBins = defaultNumberOfBins;
            this->binSize = interval / this->numberOfBins;
        }
    }

private:
    virtual TH1D* CreateRootHistogram() const
    {
        const char* _name = this->Name().c_str();
        return new TH1D(_name, _name, this->NumberOfBins() + 1, this->LowLimit(), this->HighLimit() + this->binSize);
    }
};

template<typename ValueType>
class Base1DIntegerHistogram : public Base1DHistogram<ValueType> {
public:
    Base1DIntegerHistogram(const std::string& name) : Base1DHistogram<ValueType>(name) {}
    Base1DIntegerHistogram(const std::string& name, size_t nbins, const ValueType& low, const ValueType& high) :
        Base1DHistogram<ValueType>(name, nbins, low, high) {}

    virtual void Adjust()
    {
        static const size_t defaultBinSize = 1;
        if(!this->lowLimitIsFixed)
            this->limits.first = this->minMax.first;
        if(!this->highLimitIsFixed)
            this->limits.second = this->minMax.second;

        const size_t interval = this->limits.second >= this->limits.first ?
                    static_cast<size_t>(this->limits.second - this->limits.first) + 1 : 1;
        if(this->useBinSize) {
            if(!this->binSizeIsFixed)
                this->binSize = defaultBinSize;
            this->numberOfBins = interval / this->binSize;
        }
        else {
            if(!this->numberOfBinsIsFixed)
                this->numberOfBins = interval;
            this->binSize = interval / this->numberOfBins;
        }
    }

private:
    virtual TH1D* CreateRootHistogram() const
    {
        const char* _name = this->Name().c_str();
        return new TH1D(_name, _name, this->NumberOfBins(), this->LowLimit() - 0.5, this->HighLimit() + 0.5);
    }
};

template<typename ValueType>
class Base2DHistogram : public AbstractHistogram, public HistogramProperties<ValueType> {
public:
    Base2DHistogram(const std::string& name) : AbstractHistogram(name) {}

    virtual ~Base2DHistogram()
    {
        delete x_histogram;
        delete y_histogram;
    }

    void Fill(const ValueType& x, const ValueType& y)
    {
        x_histogram->Fill(x);
        y_histogram->Fill(y);
    }

    TH2D* ProduceRootHistogram()
    {
        Adjust();
        TH2D* rootHistogram = CreateRootHistogram();
        auto x_iter = x_histogram->begin();
        auto y_iter = y_histogram->begin();
        for(; x_iter != x_histogram->end(); ++x_iter, ++y_iter)
            rootHistogram->Fill(*x_iter, *y_iter);
        return rootHistogram;
    }

    void WriteRootHistogram()
    {
        TH2D* h = ProduceRootHistogram();
        h->Write();
        delete h;
    }

    virtual void Adjust()
    {
        x_histogram->Adjust();
        y_histogram->Adjust();
    }

    HistogramProperties<ValueType>& XAxisProperties() { return *x_histogram; }
    const HistogramProperties<ValueType>& XAxisProperties() const { return *x_histogram; }
    HistogramProperties<ValueType>& YAxisProperties() { return *y_histogram; }
    const HistogramProperties<ValueType>& YAxisProperties() const { return *y_histogram; }

private:
    virtual TH2D* CreateRootHistogram() const = 0;

protected:
    Base1DHistogram<ValueType> *x_histogram, *y_histogram;
};

template<typename ValueType>
class Base2DFloatingPointHistogram : public Base2DHistogram<ValueType> {
public:
    Base2DFloatingPointHistogram(const std::string& name)
        : Base2DHistogram<ValueType>(name)
    {
        this->x_histogram = new Base1DFloatingPointHistogram<ValueType>(name);
        this->y_histogram = new Base1DFloatingPointHistogram<ValueType>(name);
    }

private:
    virtual TH2D* CreateRootHistogram() const
    {
        const char* _name = this->Name().c_str();
        const HistogramProperties<ValueType>& x = this->XAxisProperties();
        const HistogramProperties<ValueType>& y = this->YAxisProperties();
        return new TH2D(_name, _name,
                        x.NumberOfBins() + 1, x.LowLimit(), x.HighLimit() + x.BinSize(),
                        y.NumberOfBins() + 1, y.LowLimit(), y.HighLimit() + y.BinSize());
    }
};

template<typename ValueType>
class Base2DIntegerHistogram : public Base2DHistogram<ValueType> {
public:
    Base2DIntegerHistogram(const std::string& name)
        : Base2DHistogram<ValueType>(name)
    {
        this->x_histogram = new Base1DIntegerHistogram<ValueType>(name);
        this->y_histogram = new Base1DIntegerHistogram<ValueType>(name);
    }

private:
    virtual TH2D* CreateRootHistogram() const
    {
        const char* _name = this->Name().c_str();
        const HistogramProperties<ValueType>& x = this->XAxisProperties();
        const HistogramProperties<ValueType>& y = this->YAxisProperties();
        return new TH2D(_name, _name,
                        x.NumberOfBins(), x.LowLimit() - 0.5, x.HighLimit() + 0.5,
                        y.NumberOfBins(), y.LowLimit() - 0.5, y.HighLimit() + 0.5);
    }
};

} // namespace detail

template<typename ValueType>
class SmartHistogram {};

template<>
class SmartHistogram<double> : public detail::Base1DFloatingPointHistogram<double> {
public:
    SmartHistogram(const std::string& name) : Base1DFloatingPointHistogram<double>(name) {}
    SmartHistogram(const std::string& name, size_t nbins, double low, double high) :
        Base1DFloatingPointHistogram<double>(name, nbins, low, high) {}
};

template<>
class SmartHistogram<float> : public detail::Base1DFloatingPointHistogram<float> {
public:
    SmartHistogram(const std::string& name) : Base1DFloatingPointHistogram<float>(name) {}
    SmartHistogram(const std::string& name, size_t nbins, float low, float high) :
        Base1DFloatingPointHistogram<float>(name, nbins, low, high) {}
};

template<>
class SmartHistogram<int> : public detail::Base1DIntegerHistogram<int> {
public:
    SmartHistogram(const std::string& name) : Base1DIntegerHistogram<int>(name) {}
};

template<>
class SmartHistogram<bool> : public detail::Base1DIntegerHistogram<bool> {
public:
    SmartHistogram(const std::string& name) : Base1DIntegerHistogram<bool>(name) {}
};

template<>
class SmartHistogram< std::pair<double, double> > : public detail::Base2DFloatingPointHistogram<double> {
public:
    SmartHistogram(const std::string& name) : Base2DFloatingPointHistogram<double>(name) {}
};

template<>
class SmartHistogram< std::pair<float, float> > : public detail::Base2DFloatingPointHistogram<float> {
public:
    SmartHistogram(const std::string& name) : Base2DFloatingPointHistogram<float>(name) {}
};

template<>
class SmartHistogram< std::pair<int, int> > : public detail::Base2DIntegerHistogram<int> {
public:
    SmartHistogram(const std::string& name) : Base2DIntegerHistogram<int>(name) {}
};

template<>
class SmartHistogram< std::pair<bool, bool> > : public detail::Base2DIntegerHistogram<bool> {
public:
    SmartHistogram(const std::string& name) : Base2DIntegerHistogram<bool>(name) {}
};

template<>
class SmartHistogram<TH1D> : public TH1D, public AbstractHistogram {
public:
    SmartHistogram(const std::string& name, size_t nbins, double low, double high)
        : TH1D(name.c_str(), name.c_str(), nbins, low, high), AbstractHistogram(name) {}

    virtual void WriteRootHistogram()
    {
        this->Write();
    }
};

template<>
class SmartHistogram<TH2D> : public TH2D, public AbstractHistogram {
public:
    SmartHistogram(const std::string& name,
                   size_t nbinsx, double xlow, double xup,
                   size_t nbinsy, double ylow, double yup)
        : TH2D(name.c_str(), name.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup), AbstractHistogram(name) {}

    virtual void WriteRootHistogram()
    {
        this->Write();
    }
};

template<typename ValueType>
struct HistogramFactory {
    template<typename ...Args>
    static SmartHistogram<ValueType>* Make(const std::string& name, Args... args)
    {
        return new SmartHistogram<ValueType>(name, args...);
    }
};

} // root_ext
