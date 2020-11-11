#ifndef THREEDGRID_H
#define THREEDGRID_H

#include <vector>
#include <string>

using namespace std;


class ThreeDGrid
{
    public:
        ThreeDGrid();
        virtual ~ThreeDGrid();
        ThreeDGrid(const ThreeDGrid& other);
        ThreeDGrid& operator=(const ThreeDGrid& other);
        bool ReadFromModFlowFile(const string &filename);
        bool ReadFromModFlowFile(const string &filename, double time);
        void setdimension(int nx, int ny, int nz);
        void setincrements(double _dx, double _dy, double _dz) {dz=_dz; dx=_dx; dy=_dy;}
        bool create_vts_test(const string &filename);
        bool create_vts(const string &filename, int x_limit=0, int y_limit=0, int z_limit=0, int x_inteval=0, int y_interval=0, int z_interval=0);

    protected:

    private:
    vector<vector<vector<double>>> values;
    int nx;
    int ny;
    int nz;
    double dx, dy, dz;
};

#endif // THREEDGRID_H
