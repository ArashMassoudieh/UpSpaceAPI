#ifndef MAPASTIMESERIESSET_H
#define MAPASTIMESERIESSET_H

#include "BTCSet.h"


class MapAsTimeSeriesSet: public CTimeSeriesSet
{
    public:
        MapAsTimeSeriesSet();
        virtual ~MapAsTimeSeriesSet();
        MapAsTimeSeriesSet(const MapAsTimeSeriesSet& other);
        MapAsTimeSeriesSet& operator=(const MapAsTimeSeriesSet& other);
        void append(CTimeSeries& _BTC, double x);
        bool writetofile(const string &filename);
        vector<double> x;
    protected:

    private:
};

#endif // MAPASTIMESERIESSET_H
