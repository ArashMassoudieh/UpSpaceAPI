#pragma once

#include <string>
#include <vector>
#include "Vector.h"
#include "Matrix.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_rng.h"


class CCopula
{
public:
	string copula;
	vector<double> parameters;
	double evaluate11(double u1, double u2);
	CCopula();
	~CCopula();
	void SetCorrelation(const double &r);
	double get_random_at(const double &u1);
	double w;
	gsl_rng *RngPtr() {return rng_ptr;}
private:
    gsl_rng *rng_ptr;
    CMatrix M_inv;
    double correlation;

};

