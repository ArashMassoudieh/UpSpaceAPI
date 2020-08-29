#include "OneDGrid.h"

#include "NormalDist.h"
#include "BTC.h"
#include "Matrix_arma.h"
#include "Vector_arma.h"
#include "StringOP.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "Distribution.h"
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
    u = CVector(n);

    omega[0] = getnormalrand(0,1);
    double gamma = exp(-pow(delta_x,2.0)/pow(correlation_ls,2.0));
    for (int i=1; i<n; i++)
    {
        omega[i] = gamma*omega[i-1] + sqrt(1.0-pow(gamma,2))*getnormalrand(0,1);

    }
    omega = NormalizetoGaussian(omega);

    for (int i=1; i<n; i++)
        u[i] = gsl_cdf_gaussian_P(omega[i],1);

    return true;
}

bool OneDGrid::Generate_u_field(int n)
{
    u = CVector(n);

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

void OneDGrid::write_u_field(const string &filename)
{
    CBTC towrite;
    for (int i=0; i<omega.num; i++)
    {
        towrite.append((i+0.5)*delta_x,u[i]);
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


CTimeSeries OneDGrid::getDistU(int nbins)
{
    CBTC towrite;
    for (int i=0; i<u.num; i++)
    {
        towrite.append((i+0.5)*delta_x,u[i]);
    }
    steady_state_u_dist = towrite.distribution(nbins);
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

void OneDGrid::AssignConcentration_based_on_u(double u_min, double u_max, double val)
{
    for (int i=0; i<omega.num; i++)
        if (u[i]<=u_max && u[i]>=u_min)
            C[i] = val;
}


CTimeSeries OneDGrid::GetConcentrationDistributionOverOmega(int nbins)
{
    CBTC omegadist = getDist(nbins);
    CBTC cdist = omegadist;
    cdist = 0;
    CBTC wdist = cdist;
    double dp = omegadist.t[1]-omegadist.t[0];
    double p_start = omegadist.t[0] + dp/2;
    for (int i=0; i<C.num; i++)
    {   cdist.C[int((omega[i]-p_start)/dp)+1] += C[i];
        wdist.C[int((omega[i]-p_start)/dp)+1] += 1;
    }

    return cdist/wdist.integrate();
}

CTimeSeries OneDGrid::GetConcentrationDistributionOverU(int nbins)
{
    CBTC udist = getDistU(nbins);
    CBTC cdist = udist;
    cdist = 0;
    CBTC wdist = cdist;
    double dp = udist.t[1]-udist.t[0];
    double p_start = udist.t[0] + dp/2;
    for (int i=0; i<C.num; i++)
    {   cdist.C[int((u[i]-p_start)/dp)+1] += C[i];
        wdist.C[int((u[i]-p_start)/dp)+1] += 1;
    }

    return cdist/wdist.integrate();
}


void OneDGrid::Solve(const double &dt, const double &t_end, const double &D, int n_w_bins)
{
    w = 1;
    set_progress_value(0);
    M = getTransportMatrix(dt, D);
    cout<<"Matrix created!"<<endl;
    ANSCX = CBTCSet(omega.num);
    for (int i=0; i<omega.num; i++)
    {
        ANSCX.names[i] = ("x=" + numbertostring((i - 0.5)*delta_x));
    }

    ANSCX.append(0, C.vec);
    ANSCW.append(GetConcentrationDistributionOverOmega(n_w_bins),"t = " + numbertostring(0));
    ANSCU.append(GetConcentrationDistributionOverU(n_w_bins),"t = " + numbertostring(0));
    u_dev.append(0, diff(GetConcentrationDistributionOverU(n_w_bins),steady_state_u_dist));
    for (double t = dt; t<t_end; t+=dt)
    {
        set_progress_value(t);
        OneStepSolve(dt, D);
        ANSCX.append(t, C.vec);
        ANSCW.append(GetConcentrationDistributionOverOmega(n_w_bins),"t = " + numbertostring(t));
        ANSCU.append(GetConcentrationDistributionOverU(n_w_bins),"t = " + numbertostring(t));
        u_dev.append(t, diff(GetConcentrationDistributionOverU(n_w_bins),steady_state_u_dist));
    }
    u_dev.append(t_end, diff(GetConcentrationDistributionOverU(n_w_bins),steady_state_u_dist));

}

void OneDGrid::OneStepSolve(const double &dt, const double &D)
{

    CVector_arma RHS(omega.num);

    for (int i=1; i<omega.num-1; i++)
    {
        RHS[i] += (1/dt - 2*(1-w)*D/pow(delta_x,2))*C[i];
        RHS[i] += (1-w)*D/pow(delta_x,2)*C[i-1];
        RHS[i] += (1-w)*D/pow(delta_x,2)*C[i+1];
    }

    RHS[0]=0;
    RHS[omega.num-1]=0;

    CVector_arma ANS = RHS/M;
    C = ANS;

}

CMatrix_arma OneDGrid::getTransportMatrix(const double &dt, const double &D)
{
    CMatrix_arma M(omega.num);
    for (int i=1; i<omega.num-1; i++)
    {
        M(i,i) = 1/dt + 2*w*D/pow(delta_x,2);
        M(i,i-1) = -w*D/pow(delta_x,2);
        M(i,i+1) = -w*D/pow(delta_x,2);
    }

    M(0,1) = 1;
    M(0,0)= -1;

    M(omega.num-1,omega.num-2) = 1;
    M(omega.num-1,omega.num-1) = -1;

    return M;
}



void OneDGrid::Initialize_U(int i, double value)
{
    u[i]=value;
}

void OneDGrid::Initialize_U(const CVector &values)
{
    u=values;
}


void OneDGrid::Solve_U(const double &dt, const double &t_end, const double &D, const double &correlation)
{

    ANSCU = CBTCSet(u.num);
    du = 1.0/double(u.num);
    for (int i=0; i<u.num; i++)
    {
        ANSCU.names[i] = ("x=" + numbertostring((i - 0.5)*du));
    }


    CVector exchange = create_ou_exchange(D, correlation);
    CMatrix_arma M(u.num);
    for (int i=1; i<u.num-1; i++)
    {
        M(i,i) = 1/dt + w*(exchange[i-1]+exchange[i])/pow(du,2);
        M(i,i-1) = -w*exchange[i-1]/pow(du,2);
        M(i,i+1) = -w*exchange[i]/pow(du,2);
    }

    M(0,1) = -w*(exchange[0])/pow(du,2);
    M(0,0)= 1/dt + w*(exchange[0])/pow(du,2);

    M(u.num-1,u.num-2) = -w*(exchange[u.num-2])/pow(du,2);
    M(u.num-1,u.num-1) = 1/dt + w*(exchange[u.num-2])/pow(du,2);
    ANSCU.append(0,u.vec);
    u_dev.append(0,(u-1).norm2());
    CVector_arma RHS(u.num);
    for (double t=dt; t<t_end; t+=dt)
    {
        for (int i=1; i<u.num-1; i++)
        {
            RHS[i] = 1/dt*u[i] - (1-w)*(exchange[i-1]+exchange[i])/pow(du,2)*u[i];
            RHS[i] += (1-w)*exchange[i-1]/pow(du,2)*u[i-1];
            RHS[i] += (1-w)*exchange[i]/pow(du,2)*u[i+1];
        }

        RHS[0] = (1-w)*(exchange[0])/pow(du,2)*u[1];
        RHS[0] += 1/dt*u[0] - (1-w)*(exchange[0])/pow(du,2)*u[0];

        RHS[u.num-1] = (1-w)*(exchange[u.num-2])/pow(du,2)*u[u.num-2];
        RHS[u.num-1] += 1/dt*u[u.num-1] - (1-w)*(exchange[u.num-2])/pow(du,2)*u[u.num-1];

        CVector_arma ANS = (RHS/M);
        u = ANS;
        ANSCU.append(t,u.vec);
        u_dev.append(t,(u-1).norm2());

    }

}


CVector OneDGrid::create_ou_exchange(const double &D, const double &correlation)
{
	CVector exchanges(u.num-1);
	for (int j = 0; j < u.num-1; j++)
	{
        exchanges[j] = (2.0*D/pow(correlation,2))*pow(std_normal_phi_inv((double(j+1))*du), 2);
	}
	return exchanges;
}
