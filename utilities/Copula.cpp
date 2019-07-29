#include "Copula.h"
#include "NormalDist.h"


CCopula::CCopula()
{
}


CCopula::~CCopula()
{
}

double CCopula::evaluate11(double u1, double u2)
{
	double a1;
	if (copula == "gaussian")
	{
		double x = exp(-parameters[0]);
		CMatrix M_inv(2);
		M_inv[0][0] = x*x / (1 - x*x);
		M_inv[1][1] = x*x / (1 - x*x);
		M_inv[0][1] = -x / (1 - x*x);
		M_inv[1][0] = -x / (1 - x*x);
		CVector Y0(2);
		Y0[0] = stdnormal_inv(u1);
		Y0[1] = stdnormal_inv(u2);
		a1 = 1/sqrt(1 - pow(x,2))*exp(-0.5*dotproduct((M_inv*Y0), Y0));
	}
	return a1;

}
