#include "Pathway.h"
#include <iostream>
#include <fstream>
#include "gsl/gsl_cdf.h"
#include "Copula.h"



CPathway::CPathway()
{
}

CPathway::CPathway(const CPathway &P)
{
	maxx = P.maxx;
	minx = P.minx;
	positions = P.positions;
	uniform = P.uniform;
    weight = P.weight;
}

void CPathway::create_range_x(double x_min, double x_max, double dx)
{


}

CPathway& CPathway::operator=(const CPathway &P)
{
	maxx = P.maxx;
	minx = P.minx;
	positions = P.positions;
    weight = P.weight;
    uniform = P.uniform;
	return *this;
}



CPathway::~CPathway()
{
}

void CPathway::append(CPosition pos)
{
	positions.push_back(pos);
}

CVector CPathway::get_velocity_at_x(double x)
{
	for (int i = 0; i < int(positions.size())-1; i++)
		if (positions[i].x<x && positions[i + 1].x>=x)
			return (x - positions[i].x) / (positions[i + 1].x - positions[i].x)*(positions[i + 1].v - positions[i].v) + positions[i].v;

	return CVector();
}

CPosition CPathway::get_position_at_x(double x)
{
	CPosition p;
	for (int i = 0; i < int(positions.size()) - 1; i++)
		if (positions[i].x<x && positions[i + 1].x>=x)
		{
			p = (x - positions[i].x) / (positions[i + 1].x - positions[i].x)*(positions[i + 1] - positions[i]) + positions[i];
			return p;
		}

	return CPosition();
}

CPosition CPathway::get_position_at_t(double t)
{
	CPosition p;
	for (int i = 0; i < int(positions.size()) - 1; i++)
		if (positions[i].t<t && positions[i + 1].t>t)
		{
			p = (t - positions[i].t) / (positions[i + 1].t - positions[i].t)*(positions[i + 1] - positions[i]) + positions[i];
			return p;
		}

	return CPosition();
}

CVector CPathway::get_velocity_at_t(double t)
{
	for (int i = 0; i < int(positions.size()) - 1; i++)
		if (positions[i].t<t && positions[i + 1].t>t)
			return (t - positions[i].t) / (positions[i + 1].t - positions[i].t)*(positions[i + 1].t - positions[i].t) + positions[i].t;
	return CVector();
}

CPathway CPathway::make_uniform_x(double dx)
{
	CPathway pathout;
	double backward = 1;
	if (positions[positions.size()-1].x>positions[0].x)
        backward = 1;
    else
        backward = -1;

	pathout.uniform = true;
    pathout.weight = weight;
	double x;
	if (positions.size()>0)
    {
        if (backward == 1)
        {
            x = positions[0].x;
            pathout.append(positions[0]);
        }
        else
        {
            x = positions[positions.size()-1].x;
            pathout.append(positions[positions.size()-1]);
        }
    }
    else
        return pathout;

	x += dx;
	if (backward == 1)
    {
        for (int i = 1; i < int(positions.size()); i++)
        {

            while (positions[i].x > x)
            {
                CPosition p = positions[i - 1] + (positions[i] - positions[i - 1]) / (positions[i].x - positions[i - 1].x)*(x - positions[i - 1].x);
                pathout.append(p);
                x += dx;
            }
        }
    }
    else
    {
        for (int i = int(positions.size()-1); i >= 0 ; i--)
        {
            while (positions[i].x > x)
            {
                CPosition p = positions[i - 1] + (positions[i] - positions[i - 1]) / (positions[i].x - positions[i - 1].x)*(x - positions[i - 1].x);
                pathout.append(p);
                x += dx;
            }
        }
    }
	return pathout;
}

CPathway CPathway::make_uniform_t(double dt)
{
	CPathway pathout;
	double t = positions[0].t;
        pathout.weight = weight;
	pathout.positions[0] = positions[0];
	t += dt;
	for (int i = 1; i < int(positions.size()); i++)
	{

		while (positions[i].t > t)
		{
			pathout.append(positions[i - 1] + (positions[i] - positions[i - 1]) / (positions[i].t - positions[i - 1].t)*(t - positions[i - 1].t));
			t += dt;
		}
	}

	return pathout;
}

double CPathway::max_x()
{
	if (maxx != minx) return maxx;
	maxx = -1e12;
	for (int i = 0; i < int(positions.size()); i++)
		maxx = max(positions[i].x, maxx);
	return maxx;
}

double CPathway::min_x()
{
	if (maxx != minx) return minx;
	minx = 1e12;
	for (int i = 0; i < int(positions.size()); i++)
		minx = min(positions[i].x, minx);
	return minx;
}

double CPathway::max_t()
{
	if (maxt != mint) return maxt;
	maxt = 1e12;
	for (int i = 0; i < int(positions.size()); i++)
		maxx = max(positions[i].t, maxt);
	return maxt;
}

double CPathway::min_t()
{
	if (maxt != mint) return mint;
	mint = 1e12;
	for (int i = 0; i < int(positions.size()); i++)
		minx = min(positions[i].t, mint);
	return mint;
}

void CPathway::write(string filename)
{
	ofstream file;
	file.open(filename.c_str());
	file << "t,x,y,vx,vy,u,z" << endl;
	for (int i = 0; i < int(positions.size()); i++)
	{
		file <<  positions[i].t << "," << positions[i].x << "," << positions[i].y << "," << positions[i].v[0] << "," << positions[i].v[1] << "," << positions[i].u << "," << positions[i].z << endl;
	}
	file.close();

}

void CPathway::create_ou(CDistribution *dist, double x_min, double x_max, double kappa, double dx, double _weight)
{
	CPosition p;
    weight = _weight;
	p.u = unitrandom();
	p.z = gsl_cdf_gaussian_Pinv(p.u,1);
	p.v[0] = dist->inverseCDF(p.u);
	p.x = x_min;
	p.y = 0;
	p.t = 0;
	append(p);
	while (p.x < x_max)
	{
		p.z += -dx*kappa*p.z + sqrt(2*kappa*dx)*getnormalrand(0, 1);
		p.u = gsl_cdf_gaussian_P(p.z, 1);
		p.x += dx;
		p.y = 0;
		p.v[0] = dist->inverseCDF(p.u);
		p.t += dx / p.v[0];
		append(p);
	}
}

void CPathway::create_copula(CDistribution *dist, double x_min, double x_max, double epsilon, double r, double dx, double _weight)
{
    CCopula C;
    C.SetCorrelation(r);
    CPosition p;
    weight = _weight;
	p.u = unitrandom();
	p.z = gsl_cdf_gaussian_Pinv(p.u, 1);

	p.v[0] = dist->inverseCDF(p.u);
	p.x = x_min;
	p.y = p.u;
	p.t = 0;
	append(p);
	while (p.x < x_max)
	{
		if (unitrandom()<1-exp(-dx/epsilon))
        {
            p.u = C.get_random_at(p.u);
            p.z = C.w;
            p.v[0] = dist->inverseCDF(p.u);
        }

		p.x += dx;
		p.y = p.u;
		p.t += dx / p.v[0];
		append(p);
	}
}

void CPathway::create_copula(CDistribution *dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double _weight)
{
    CPosition p;
    weight = _weight;
	p.u = unitrandom();
	p.z = gsl_cdf_gaussian_Pinv(p.u, 1);
	//p.v[0] = dist->inverseCDF(p.u,true);
	//*** after change ***///
	p.v[0] = dist->inverseCDF(p.u,false);
	p.u = dist->evaluate_CDF(p.v[0],true);
	double v1 = p.v[0];
	double v2 = v1;
	p.x = x_min;
	p.y = p.u;
	p.t = 0;
	append(p);
	while (p.x < x_max)
	{
		if (unitrandom()<1-exp(-dx/epsilon))
        {
            p.u = copula->get_random_at(p.u);
            p.z = copula->w;
            p.v[0] = dist->inverseCDF(p.u,true);
            v2 = p.v[0];
            p.y = p.u;

        }

		p.x += dx;
		p.t += (v1+v2)*dx / (2*v1*v2) ;
		append(p);
	}
}

vtkSmartPointer<vtkPolyData> CPathway::pathway_vtk_pdt_vtp(double z_factor, double offset)
{
	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> values_vx =
		vtkSmartPointer<vtkFloatArray>::New();
	values_vx->SetNumberOfComponents(1);
	values_vx->SetName("vx");

	vtkSmartPointer<vtkFloatArray> values_u =
		vtkSmartPointer<vtkFloatArray>::New();
	values_u->SetNumberOfComponents(1);
	values_u->SetName("u");

	vtkSmartPointer<vtkFloatArray> values_z =
		vtkSmartPointer<vtkFloatArray>::New();
	values_z->SetNumberOfComponents(1);
	values_z->SetName("z");

	vtkSmartPointer<vtkFloatArray> values_v_eff =
		vtkSmartPointer<vtkFloatArray>::New();
	values_v_eff->SetNumberOfComponents(1);
	values_v_eff->SetName("v_eff");

	vtkSmartPointer<vtkFloatArray> values_t_eff =
		vtkSmartPointer<vtkFloatArray>::New();
	values_t_eff->SetNumberOfComponents(1);
	values_t_eff->SetName("t_eff");

	vtkSmartPointer<vtkFloatArray> values_t =
		vtkSmartPointer<vtkFloatArray>::New();
	values_t->SetNumberOfComponents(1);
	values_t->SetName("arrival_time");


	for (int i = 0; i<int(positions.size()); i++)
	{
		double vx[1] = { positions[i].v[0] };
		double u[1] = { positions[i].u };
		double z[1] = { positions[i].z };
		double t[1] = { positions[i].t };
		double v_eff[1] = { positions[i].getvar("v_eff") };
		double t_eff[1] = { positions[i].getvar("t_eff") };
		double p[3] = { positions[i].x , positions[i].y, positions[i].v[0] * z_factor + offset };
		points->InsertNextPoint(p);
		values_vx->InsertNextTuple(vx);
		values_z->InsertNextTuple(z);
		values_u->InsertNextTuple(u);
		values_t->InsertNextTuple(t);
		values_v_eff->InsertNextTuple(v_eff);
		values_t_eff->InsertNextTuple(t_eff);
	}
	vtkSmartPointer<vtkPolyLine> polyLine =
		vtkSmartPointer<vtkPolyLine>::New();

	polyLine->GetPointIds()->SetNumberOfIds(positions.size());
	for (unsigned int i = 0; i < positions.size(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	polyData->SetPoints(points);

	// Add the lines to the dataset
	polyData->SetLines(cells);

	polyData->GetPointData()->SetScalars(values_vx);
	polyData->GetPointData()->AddArray(values_u);
	polyData->GetPointData()->AddArray(values_z);
	polyData->GetPointData()->AddArray(values_t);
	polyData->GetPointData()->AddArray(values_v_eff);
	polyData->GetPointData()->AddArray(values_t_eff);

	// Visualization

	return polyData;
}

CBTCSet CPathway::get_distribution(bool _log, int n_bins)
{
	CBTCSet B;
	//CBTC tBTC=get_distribution("t", "vx", _log, n_bins);
	CBTC tBTC=get_distribution("t", _log, n_bins);
	CBTC xBTC=get_distribution("x", _log, n_bins);
	CBTC yBTC=get_distribution("y", _log, n_bins);
	CBTC vxBTC=get_distribution("vx", _log, n_bins);
	CBTC vyBTC=get_distribution("vy", _log, n_bins);
	CBTC uBTC=get_distribution("u", _log, n_bins);
	CBTC zBTC=get_distribution("z", _log, n_bins);
	CBTC v_effBTC=get_distribution("v_eff", _log, n_bins);
	CBTC t_effBTC=get_distribution("t_eff", _log, n_bins);


	B.append(tBTC);
	B.append(xBTC);
	B.append(yBTC);
	B.append(vxBTC);
	B.append(vyBTC);
	B.append(uBTC);
	B.append(zBTC);
	B.append(v_effBTC);
	B.append(t_effBTC);

	B.setname(0, "t");
	B.setname(1, "x");
	B.setname(2, "y");
	B.setname(3, "vx");
	B.setname(4, "vy");
	B.setname(5, "u");
	B.setname(6, "z");
	B.setname(7, "v_eff");
	B.setname(8, "t_eff");

	return B;

}

CBTC CPathway::get_distribution(string var, bool _log, int n_bins)
{
	CBTC out(n_bins + 2);

	if (!_log)
	{
		vector<double> min_max = minmax(var);
		double p_start = min_max[0];
		double p_end = min_max[1] * 1.001;
		double dp = abs(p_end - p_start) / n_bins;
		if (dp == 0){
                    CBTC X = CBTC();
                    return X;
                }
		out.t[0] = p_start - dp / 2;

		for (int i = 0; i < n_bins + 1; i++)
			out.t[i + 1] = out.t[i] + dp;


		for (int i = 0; i < int(positions.size()); i++)
			out.C[int((positions[i].getvar(var) - p_start) / dp) + 1] += 1.0 / positions.size() / dp;
	}
	else
	{
		vector<double> min_max = minmax(var);
		double p_start = min_max[0];
		if (p_start <= 0){
                    CBTC X = CBTC();
                    return X;
                }
		if (min_max[0] == min_max[1]) return out;
		double p_end = min_max[1] * 1.001;
		double dp = (log(p_end) - log(p_start)) / n_bins;
		if (dp == 0)
                {
                    CBTC X = CBTC();
                    return X;
                }
		out.t[0] = exp(log(p_start) - dp / 2.0);

		for (int i = 0; i < n_bins + 1; i++)
			out.t[i + 1] = out.t[i] * exp(dp);

		for (int i = 0; i < positions.size(); i++)
			out.C[int((log(positions[i].getvar(var)) - log(p_start)) / dp) + 1] += 1.0 / positions.size() / (out.t[int((log(positions[i].getvar(var)) - log(p_start)) / dp) + 2]- out.t[int((log(positions[i].getvar(var)) - log(p_start)) / dp)])*2;
	}
	return out;
}


CBTC CPathway::get_distribution(string var, string weight_var, bool _log, int n_bins)
{
	CBTC out(n_bins + 2);

	if (!_log)
	{
		vector<double> min_max = minmax(var);
		double p_start = min_max[0];
		double p_end = min_max[1] * 1.001;
		double dp = abs(p_end - p_start) / n_bins;
		if (dp == 0){
                    CBTC X = CBTC();
                    return X;
                }
		out.t[0] = p_start - dp / 2;

		for (int i = 0; i < n_bins + 1; i++)
			out.t[i + 1] = out.t[i] + dp;

        double sum_weights = 0;
        for (int i = 0; i < int(positions.size()); i++)
		{
		    sum_weights += positions[i].getvar(weight_var,true);
		}

		for (int i = 0; i < int(positions.size()); i++)
		{
		    out.C[int((positions[i].getvar(var) - p_start) / dp) + 1] += positions[i].getvar(weight_var,true) / sum_weights / dp;
		}
	}
	else
	{
		vector<double> min_max = minmax(var);
		double p_start = min_max[0];
		if (p_start <= 0){
                    CBTC X = CBTC();
                    return X;
                }
		if (min_max[0] == min_max[1]) return out;
		double p_end = min_max[1] * 1.001;
		double dp = (log(p_end) - log(p_start)) / n_bins;
		if (dp == 0)
                {
                    CBTC X = CBTC();
                    return X;
                }
		out.t[0] = exp(log(p_start) - dp / 2.0);

        double sum_weights = 0;
        for (int i = 0; i < int(positions.size()); i++)
		{
		    sum_weights += positions[i].getvar(weight_var,true);
		}

		for (int i = 0; i < n_bins + 1; i++)
			out.t[i + 1] = out.t[i] * exp(dp);

		for (int i = 0; i < positions.size(); i++)
			out.C[int((log(positions[i].getvar(var)) - log(p_start)) / dp) + 1] += positions[i].getvar(weight_var,true) / sum_weights / (out.t[int((log(positions[i].getvar(var)) - log(p_start)) / dp) + 2]- out.t[int((log(positions[i].getvar(var)) - log(p_start)) / dp)])*2;
	}
	return out;
}

vector<double> CPathway::minmax(string var)
{
	vector<double> min_max(2);
	min_max[0] = 1e23;
	min_max[1] = -1e23;
	for (int i = 0; i < positions.size(); i++)
	{
		min_max[0] = min(min_max[0], positions[i].getvar(var));
		min_max[1] = max(min_max[1], positions[i].getvar(var));
	}

	return min_max;
}

double CPathway::get_cross_time(double xx)
{
    if (uniform)
    {
        double dx = positions[1].x - positions[0].x;
        int i = int(xx / dx);
        if (i < positions.size() - 1)
            return positions[i].t + (positions[i + 1].t - positions[i].t) / dx*(xx - positions[i].x);
        else
            return positions[positions.size()-2].t + (positions[positions.size() - 1].t - positions[positions.size() - 2].t) / dx*(xx - positions[positions.size() - 2].x);
}
    else
    {
        for (int i = 0; i < positions.size()-1; i++)
        {
            if (xx < positions[i + 1].x && xx >= positions[i].x)
            {
                return positions[i].t + (positions[i + 1].t - positions[i].t) / (positions[i+1].x-positions[i].x)*(xx - positions[i].x);
            }
        }
        if (xx > positions[positions.size() - 1].x)
        {
            return positions[positions.size() - 2].t + (positions[positions.size() - 1].t - positions[positions.size() - 2].t) / (positions[positions.size() - 1].x - positions[positions.size() - 2].x)*(xx - positions[positions.size() - 2].x);
        }
    }
}


vector<double> CPathway::get_cross_time_vx(double xx)
{
    vector<double> out(2);
    if (uniform)
    {
        double dx = positions[1].x - positions[0].x;
        int i = int(xx / dx);
        if (i < positions.size() - 1)
        {
            out[0] = positions[i].t + (positions[i + 1].t - positions[i].t) / dx*(xx - positions[i].x);
            out[1] = positions[i].v[0] + (positions[i + 1].v[0] - positions[i].v[0]) / dx*(xx - positions[i].x);
            return out;
        }
        else
        {
            out[0] = positions[positions.size()-2].t + (positions[positions.size() - 1].t - positions[positions.size() - 2].t) / dx*(xx - positions[positions.size() - 2].x);
            out[1] = positions[positions.size()-2].v[0] + (positions[positions.size() - 1].v[0] - positions[positions.size() - 2].v[0]) / dx*(xx - positions[positions.size() - 2].x);
            return out;
        }
}
    else
    {
        for (int i = 0; i < positions.size()-1; i++)
        {
            if (xx < positions[i + 1].x && xx >= positions[i].x)
            {
                out[0] = positions[i].t + (positions[i + 1].t - positions[i].t) / (positions[i+1].x-positions[i].x)*(xx - positions[i].x);
                out[1] = positions[i].v[0] + (positions[i + 1].v[0] - positions[i].v[0]) / (positions[i+1].x-positions[i].x)*(xx - positions[i].x);
                return out;
            }
        }
        if (xx > positions[positions.size() - 1].x)
        {
            out[0] = positions[positions.size() - 2].t + (positions[positions.size() - 1].t - positions[positions.size() - 2].t) / (positions[positions.size() - 1].x - positions[positions.size() - 2].x)*(xx - positions[positions.size() - 2].x);
            out[1] = positions[positions.size() - 2].v[0] + (positions[positions.size() - 1].v[0] - positions[positions.size() - 2].v[0]) / (positions[positions.size() - 1].x - positions[positions.size() - 2].x)*(xx - positions[positions.size() - 2].x);
            return out;
        }
    }
}
