#ifndef THREEDGRID_H
#define THREEDGRID_H

#include<vector>

class ThreeDGrid
{
    public:
        ThreeDGrid();
        virtual ~ThreeDGrid();
        ThreeDGrid(const ThreeDGrid& other);
        ThreeDGrid& operator=(const ThreeDGrid& other);
        bool ReadFromModFlowFile(const string &filename);
        void setdimension(int nx, int ny, int nz);
    protected:

    private:
    vector<vector<vector<double>>> values;
    int nx;
    int ny;
    int nz;
};

#endif // THREEDGRID_H
