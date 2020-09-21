#include "Copula.h"
#include "NormalDist.h"
#include "gsl/gsl_randist.h"
#include "StringOP.h"

CCopula::CCopula()
{
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 0);
}


CCopula::~CCopula()
{
}

double CCopula::evaluate11(double u1, double u2)
{
	double a1;
	if (tolower(copula) == "gaussian")
	{
		CVector Y0(2);
		Y0[0] = stdnormal_inv(u1);
		Y0[1] = stdnormal_inv(u2);
		a1 = 1/sqrt(1 - pow(correlation,2))*exp(-0.5*dotproduct((M_inv*Y0), Y0));
	}
	if (tolower(copula)== "frank")
	{
        a1 = evaluate_frank_copula_density(u1,u2);

	}
	return a1;

}


void CCopula::SetCorrelation(const double &r)
{
    correlation = exp(-r);
    M_inv = CMatrix(2);
    M_inv[0][0] = correlation*correlation / (1 - correlation*correlation);
	M_inv[1][1] = correlation*correlation / (1 - correlation*correlation);
	M_inv[0][1] = -correlation / (1 - correlation*correlation);
	M_inv[1][0] = -correlation / (1 - correlation*correlation);
}

double CCopula::get_random_at(const double &u1)
{
    double w1 = gsl_cdf_gaussian_Pinv(u1, 1);
    double w2 = gsl_ran_gaussian(rng_ptr,sqrt(1.0-correlation*correlation))+correlation*w1;
    w = w2;
    return gsl_cdf_gaussian_P(w2,1);
}

double CCopula::evaluate_frank_copula_density(const double &u1, const double &u2)
{
    double oneminusalpha = (-double(1.0)+exp(-Frank_copula_alpha*(u1-1))+exp(-Frank_copula_alpha*(u2-1))-exp(-Frank_copula_alpha*(u1+u2-1)))/(1-exp(-Frank_copula_alpha));

    if (oneminusalpha==0)
        cout<<"stop";
    double lnout = log(Frank_copula_alpha)-log(1-exp(-Frank_copula_alpha))-Frank_copula_alpha*(u1+u2-2)-2.0*log(oneminusalpha);
    return exp(lnout);
}

