#ifndef ONEDGRID_H
#define ONEDGRID_H

#include "Vector.h"
#include <string>
#include "BTCSet.h"

using namespace std;

class OneDGrid
{
    public:
        OneDGrid();
        virtual ~OneDGrid();
        OneDGrid(const OneDGrid& other);
        OneDGrid& operator=(const OneDGrid& other);
        double D;
        double delta_x;
        bool Generate_omega_field(int n, double delta_x, double correlation);
        void write_omega_field(const string &filename);
        CTimeSeries getDist(int nbins);
        void SetInitialCondition_Basedon_Omega(double omega_start, double omega_end);
        void SetInitialCondition_Basedon_U(double omega_start, double omega_end);
    protected:

    private:
        CVector omega;
        CVector v;
        CVector C;
        CVector u;
};

#endif // ONEDGRID_H
