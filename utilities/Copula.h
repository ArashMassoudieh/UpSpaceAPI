#pragma once

#include <string>
#include <vector>
#include "Vector.h"
#include "Matrix.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_rng.h"
#include "2DMap.h"

class CCopula
{
public:
	string copula;
	vector<double> parameters;
	double evaluate11(double u1, double u2);
	double evaluate(const double &u1, const double &u2);
	double evaluate_frank_copula_density(const double &u1, const double &u2);
	CCopula();
	~CCopula();
	void SetCorrelation(const double &r);
	void SetDiffusionParams(double D, double cls, double D_cls)
	{diffusion_coeff = D;
	if (D_cls!=0)
        diffusion_correlation_ls = D_cls;
	else
        diffusion_correlation_ls = 1000;
	correlation_ls = cls;}
	double get_random_at(const double &u1);
	double w;
	gsl_rng *RngPtr() {return rng_ptr;}
	double correlation_ls;
    double diffusion_correlation_ls;
    double diffusion_coeff;
    double Frank_copula_alpha=1;
    void SetCopulaMap(TDMap &cmap) {copula = "experimental"; copulamap = cmap; }
private:
    gsl_rng *rng_ptr;
    CMatrix M_inv;
    double correlation;
    TDMap copulamap;


};

