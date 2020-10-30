#include "ThreeDGrid.h"

ThreeDGrid::ThreeDGrid()
{
    //ctor
}

ThreeDGrid::~ThreeDGrid()
{
    //dtor
}

ThreeDGrid::ThreeDGrid(const ThreeDGrid& other)
{
    nx = other.nx;
    ny = other.ny;
    nz = other.nz;
    values = other.values;
}

ThreeDGrid& ThreeDGrid::operator=(const ThreeDGrid& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    nx = rhs.nx;
    ny = rhs.ny;
    nz = rhs.nz;
    values = rhs.values;
    return *this;
}

bool ThreeDGrid::ReadFromModFlowFile(const string &filename)
{

}

void ThreeDGrid::setdimension(int nx, int ny, int nz)
{

}
