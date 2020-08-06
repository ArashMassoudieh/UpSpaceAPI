#include "OneDGrid.h"

#include "NormalDist.h"
#include "BTC.h"

OneDGrid::OneDGrid()
{
    //ctor
}

OneDGrid::~OneDGrid()
{
    //dtor
}

OneDGrid::OneDGrid(const OneDGrid& other)
{
    omega = other.omega;
    v = other.v;
    C = other.C;
    delta_x = other.delta_x;
}

OneDGrid& OneDGrid::operator=(const OneDGrid& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    omega = rhs.omega;
    v = rhs.v;
    C = rhs.C;
    delta_x = rhs.delta_x;
    return *this;
}

bool OneDGrid::Generate_omega_field(int n, double deltax, double correlation_ls)
{
    delta_x = deltax;
    v = CVector(n);
    C = CVector(n);
    omega = CVector(n);
    v = CVector(n);
    u = CVector(n);
    C = CVector(n);
    omega[0] = getnormalrand(0,1);
    double gamma = exp(-pow(delta_x,2.0)/pow(correlation_ls,2.0));
    for (int i=1; i<n; i++)
    {
        omega[i] = gamma*omega[i-1] + sqrt(1.0-pow(gamma,2))*getnormalrand(0,1);
    }

    return true;
}

void OneDGrid::write_omega_field(const string &filename)
{
    CBTC towrite;
    for (int i=0; i<omega.num; i++)
    {
        towrite.append((i+0.5)*delta_x,omega[i]);
    }
    towrite.writefile(filename);
}

CTimeSeries OneDGrid::getDist(int nbins)
{
    CBTC towrite;
    for (int i=0; i<omega.num; i++)
    {
        towrite.append((i+0.5)*delta_x,omega[i]);
    }
    return towrite.distribution(nbins);
}

void OneDGrid::AssignConcentration(int i, double val)
{
    C[i] = val;
}

void OneDGrid::AssignConcentration(double val)
{
    C = val;
}


CTimeSeries OneDGrid::GetConcentrationDistributionOverOmega(int nbins)
{
    CBTC omegadist = getDist(nbins);
    CBTC cdist = omegadist;
    cdist = 0;
    double dp = omegadist.t[1]-omegadist.t[0];
    double p_start = omegadist.t[0] + dp/2;
    for (int i=0; i<C.num; i++)
        cdist.C[int((C[i]-p_start)/dp)+1] += C[i]/omegadist.n/dp;

    return cdist;
}
