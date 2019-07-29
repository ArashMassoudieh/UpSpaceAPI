#pragma once

#include <string>
#include <vector>
#include "Vector.h"
#include "Matrix.h"


class CCopula
{
public:
	string copula;
	vector<double> parameters;
	double evaluate11(double u1, double u2);
	CCopula();
	~CCopula();
};

