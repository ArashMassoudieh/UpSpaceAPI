#include "MapAsTimeSeriesSet.h"
#include "StringOP.h"

MapAsTimeSeriesSet::MapAsTimeSeriesSet():CTimeSeriesSet::CTimeSeriesSet()
{
    //ctor
}

MapAsTimeSeriesSet::~MapAsTimeSeriesSet()
{
    //dtor
}

MapAsTimeSeriesSet::MapAsTimeSeriesSet(const MapAsTimeSeriesSet& other):CTimeSeriesSet::CTimeSeriesSet(other)
{
    x = other.x;
}

MapAsTimeSeriesSet& MapAsTimeSeriesSet::operator=(const MapAsTimeSeriesSet& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    CTimeSeriesSet::operator=(rhs);
    x = rhs.x;
    //assignment operator

    return *this;
}

void MapAsTimeSeriesSet::append(CTimeSeries& _BTC, double xx)
{
    CTimeSeriesSet::append(_BTC, numbertostring(xx));
    x.push_back(xx);
}

bool MapAsTimeSeriesSet::writetofile(const string &filename)
{
    CTimeSeriesSet::writetofile(filename,1,true);
    return true;
}
