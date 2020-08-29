#ifndef ONEDGRID_H
#define ONEDGRID_H

#include "Vector.h"
#include <string>
#include "BTCSet.h"
#include "Matrix_arma.h"

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
        bool Generate_u_field(int n);
        void write_omega_field(const string &filename);
        void write_u_field(const string &filename);
        CTimeSeries getDist(int nbins);
        CTimeSeries getDistU(int nbins);
        void SetInitialCondition_Basedon_Omega(double omega_start, double omega_end);
        void SetInitialCondition_Basedon_U(double omega_start, double omega_end);
        void AssignConcentration(int i, double val);
        CTimeSeries GetConcentrationDistributionOverOmega(int nbins);
        CTimeSeries GetConcentrationDistributionOverU(int nbins);
        void AssignConcentration(double val);
        void OneStepSolve(const double &dt, const double &D);
        CMatrix_arma getTransportMatrix(const double &dt, const double &D);
        void Solve(const double &dt, const double &t_end, const double &D, int n_w_bins=50);
        double w = 1;
        CBTCSet ANSCX;
        CBTCSet ANSCW;
        CBTCSet ANSCU;
        CBTC u_dev;
        void Solve_U(const double &dt, const double &t_end, const double &D, const double &correlation);
        CVector create_ou_exchange(const double &D, const double &correlation);
        void Initialize_U(int i, double value);
        void Initialize_U(const CVector &values);
        void AssignConcentration_based_on_u(double u_min, double u_max, double val);

    protected:

    private:
        CVector omega;
        CVector v;
        CVector C;
        CVector u;
        CMatrix_arma M;
        CBTC steady_state_u_dist;
        double du=0;

};

#endif // ONEDGRID_H
