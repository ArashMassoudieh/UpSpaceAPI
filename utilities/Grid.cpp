#include "Grid.h"
#include "StringOP.h"
#include "NormalDist.h"
#ifdef QT_version
#include <qtextbrowser.h>
#include "heteval.h"
#include "qdebug.h"
#include <QCursor>
#endif // Qt_version
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include <math.h>
#include "MapAsTimeSeriesSet.h"
#include "omp.h"


#include <sys/resource.h>

#ifdef PowerEdge820
#define DEFAULT_FILE_PATH "/mnt/E/Projects/Upscaling/"
#endif // Arash

#ifdef PowerEdge720
#define DEFAULT_FILE_PATH "/mnt/3rd900/Projects/UpscalingFiles/"
#endif // Arash


#ifdef Arash
#define DEFAULT_FILE_PATH "/home/arash/Projects/UpscalingInputfiles/"
#endif // Arash

#ifdef Abdullelah
#define DEFAULT_FILE_PATH "/home/Abdulelah/UpscalingInputfiles/"
#endif // Abdullelah

#define symetrical
//#define Qt_version

CGrid::CGrid()
{
}

CGrid::~CGrid()
{
}

vector<ijval> CGrid::get_closest_K_dets(int i, int j, int n)
{
	vector<ijval> out;
	double max_dist = 0;
	for (int k = 1; k < max(GP.nx, GP.ny); k++)
	{
		double min_dist = 1e12;
		int jj = j - k;
		if (jj>=0)
			for (int ii = max(i - k,0); ii <= min(i + k,GP.nx-1); ii++)
			{
				double dist2 = ((i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x*(i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x + (j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y*(j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		jj = j + k;
		if (jj < GP.ny)
			for (int ii = max(i - k, 0); ii <= min(i + k, GP.nx - 1); ii++)
			{
				double dist2 = ((i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x*(i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x + (j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y*(j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		int ii = i - k;
		if (ii >= 0)
			for (int jj = max(j - k+1, 0); jj <= min(j + k-1, GP.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x*(i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x + (j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y*(j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		ii = i + k;
		if (ii < GP.nx)
			for (int jj = max(j - k + 1, 0); jj <= min(j + k - 1, GP.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x*(i - ii)*GP.dx/field_gen.k_correlation_lenght_scale_x + (j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y*(j - jj)*GP.dy/field_gen.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
                                            ijval pp;
                                            pp.i = ii;
                                            pp.j = jj;
                                            pp.val = sqrt(dist2);
                                            out.push_back(pp);
                                            max_dist = max(dist2, max_dist);
					}
				}

			}
		if ((int(out.size()) >= n) && min_dist > max_dist)
		{
                    k = max(GP.nx, GP.ny);
                    ijval pp;
                    pp.i = i;
                    pp.j = j;
                    pp.val = 0;
                    out.push_back(pp);
                    return out;
		}
	}
	ijval pp;
	pp.i = i;
	pp.j = j;
	pp.val = 0;
	out.push_back(pp);
	return out;
}

correl_mat_vec CGrid::get_correll_matrix_vec(int i, int j)
{
	correl_mat_vec M;
	vector<ijval> ij = get_top_n(get_closest_K_dets(i, j, min(n_k_dets,field_gen.max_correl_n)+1), min(n_k_dets, field_gen.max_correl_n)+1);
#ifdef  arma

	M.M_22 = CMatrix_arma(ij.size() - 1);
	M.V_21 = CVector_arma(ij.size() - 1);
	M.V_RHS = CVector_arma(ij.size() - 1);
#else
	M.M_22 = CMatrix(ij.size() - 1);
	M.V_21 = CVector(ij.size() - 1);
	M.V_RHS = CVector(ij.size() - 1);
#endif //  arma
	for (int ii = 1; ii < int(ij.size()); ii++)
	{
		M.V_21[ii-1] = exp(-sqrt((i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y) );
		M.V_RHS[ii - 1] = p[ij[ii].i][ij[ii].j].K_gauss[0];
		for (int jj = 1; jj < int(ij.size()); jj++)
		{
#ifdef arma
		M.M_22(ii-1,jj-1) = exp(-sqrt((ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y));
#else
		M.M_22[ii - 1][jj - 1] = exp(-sqrt((ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y));
#endif // arma
		}
	}
	return M;
}

void CGrid::assign_K_gauss(int i, int j)
{
	CNormalDist ND;
	correl_mat_vec M = get_correll_matrix_vec(i, j);
	double mu;
	double sigma;
	if (M.V_RHS.num == 0)
	{
		mu = 0;
		sigma = 1;
	}
	else
	{
            CMatrix_arma M_inv = inv(M.M_22);
            mu = dotproduct(M_inv*M.V_21, M.V_RHS);
            sigma = 1.0 - dotproduct(M_inv*M.V_21, M.V_21);
	}

	double K_gauss = getnormalrand(mu, sigma);
	p[i][j].k_det = true;
	n_k_dets++;
	p[i][j].K_gauss[0] = K_gauss;
	p[i][j].K[0] = map_to_KCDF(getnormalcdf(K_gauss));
}


void CGrid::assign_K_gauss()
{
	int n_filled = 0;
	srand(time(NULL));
	while (n_filled<GP.nx*GP.ny)
	{
        if (n_filled%100==0)
		set_progress_value(double(n_filled) / double(GP.nx * GP.ny));
            int i = unitrandom()*(GP.nx-1) + 0.5;
            int j = unitrandom()*(GP.ny-1) + 0.5;
            if (!p[i][j].k_det)
            {
                assign_K_gauss(i, j);
                n_filled++;
            }
    }
    cout<<endl;

}

bool CGrid::getfromfile(string filename, int _nx, int _ny)
{

    bool file_not_found;
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
			p[i][j].V = CVector(2);
	}

	ifstream file(filename);
	vector<string> s;
	if (file.good() == false)
	{
		file_not_found = true;
		return false;
	}

	double x_0, y_0;
	if (file.good())
		for (int j = 0; j < GP.ny; j++)
			for (int i = 0; i < GP.nx; i++)
			{

				s = getline(file);
				if ((i == 0) && (j == 0))
				{
					x_0 = atof(s[0].c_str());
					y_0 = atof(s[1].c_str());
				}
				if (i == 1)
					GP.dx = atof(s[0].c_str()) - x_0;

				if (j == 1)
					GP.dy = atof(s[1].c_str()) - y_0;

				p[i][j].V[0] = atof(s[2].c_str());
				p[i][j].V[1] = atof(s[3].c_str());
			}
	file.close();
	return true;

}

bool CGrid::getKfromfile(string filename, int _nx, int _ny, int _nx_act, int _ny_act)
{
	bool file_not_found;
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
		}
	}

	ifstream file(filename);
	vector<string> s;
	if (file.good() == false)
	{
		file_not_found = true;
		return false;
	}

	double x_0, y_0;
	if (file.good())
		for (int j = 0; j < _ny_act; j++)
			for (int i = 0; i < _nx_act; i++)
			{

				s = getline(file);
				if ((i == 0) && (j == 0))
				{
					x_0 = atof(s[0].c_str());
					y_0 = atof(s[1].c_str());
				}
				if (i == 1)
					GP.dx = atof(s[0].c_str()) - x_0;

				if (j == 1)
					GP.dy = atof(s[1].c_str()) - y_0;

				//p[i][j].K[0] = atof(s[2].c_str());
				if ((i<GP.nx) && (j<GP.ny)) p[i][j].K[0] = atof(s[3].c_str());
			}
	file.close();
	return true;

}

bool CGrid::createuniform(double K, int _nx, int _ny)
{
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K[0] = K;
			p[i][j].K[1] = K;
		}

	}

	return false;
}




CGrid::CGrid(string filename)
{
	ifstream file(filename);
	if (!file.good())
        {
#if QT_version
            show_in_window("File " + filename + "was not found!");
#endif
            filename = DEFAULT_FILE_PATH + filename;
            //filename = "/home/Arash/Projects/UpscalingInputfiles/" + filename;
            file.open(filename);
            if (!file.good())
            {
                cout << "File " + filename + "was not found!" << endl;
                return;
            }

        }
        else
        {
#if QT_version
            show_in_window("Reading " + filename);
#endif
            cout << "Reading " + filename + "..." << endl;

        }
	vector<string> s;

	int n_for = 1;
	int last_for_location = 0;
	while (!file.eof())
	{
		s = getline(file);
		if (s.size())
		{

			if (tolower(s[0]) == "pathout")
            {
                pathout = s[1].c_str();
                pathout.erase( std::remove(pathout.begin(), pathout.end(), '\r'), pathout.end() );
                cout<<"Output path is: " + pathout << endl;
            }
			if (tolower(s[0]) == "for")
			{
				cout << "For loop encountered"<<endl;
				n_for = atoi(s[1].c_str());
				last_for_location = commands.size();
			}
			if (tolower(s[0]) == "endfor")
			{
                _command command;
                command.command = "clear_all";
                commands.push_back(command);
                int end_for_location = commands.size();
                for (int i = 0; i < n_for-1; i++)
                {
                    _command write_commend;
                    write_commend.command = "write";
                    write_commend.parameters["content"] = "starting realization " + numbertostring(i + 1);
                    commands.push_back(write_commend);
                    for (int j = last_for_location; j < end_for_location; j++)
                    {
                        _command new_command = commands[j];
                        if (new_command.parameters["filename"].size() > 0)
                                new_command.parameters["filename"] = insert_counter_in_file_name(new_command.parameters["filename"], i+1);
                        if (new_command.parameters["filename_x"].size() > 0)
                                new_command.parameters["filename_x"] = insert_counter_in_file_name(new_command.parameters["filename_x"], i+1);
                        if (new_command.parameters["filename_y"].size() > 0)
                                new_command.parameters["filename_y"] = insert_counter_in_file_name(new_command.parameters["filename_y"], i+1);
                        if (new_command.parameters["filename_mag"].size() > 0)
                                new_command.parameters["filename_mag"] = insert_counter_in_file_name(new_command.parameters["filename_mag"], i+1);
                        if (new_command.parameters["dist_filename"].size() > 0)
                                new_command.parameters["dist_filename"] = insert_counter_in_file_name(new_command.parameters["dist_filename"], i+1);
                        if (new_command.parameters["ranks_filename"].size() > 0)
                                new_command.parameters["ranks_filename"] = insert_counter_in_file_name(new_command.parameters["ranks_filename"], i+1);
                        if (new_command.parameters["ranks_filename"].size() > 0)
                                new_command.parameters["ranks_filename"] = insert_counter_in_file_name(new_command.parameters["normal_filename"], i+1);
                        if (new_command.parameters["OU_parameters_filename"].size() > 0)
                                new_command.parameters["OU_parameters_filename"] = insert_counter_in_file_name(new_command.parameters["OU_parameters_filename"], i+1);

                        commands.push_back(new_command);

                    }
                    _command command;
                    command.command = "clear_all";
                    commands.push_back(command);

                }
                n_for = 1;

            }
            if (tolower(s[0]) == "k_dist")
            {

                marginal_K_dist_type = s[1];
                for (int j = 2; j < s.size(); j++)
                    marginal_K_dist_params.push_back(atof(s[j].c_str()));
                cout << "Hydraulic conductivity marginal distribution initialized ..." << endl;
            }
            if (tolower(s[0]) == "command")
            {
                _command command;
                command.command = tolower(s[1]);
                for (int j = 2; j < s.size(); j++)
                {   if (split(s[j], '=').size()>1)
                        command.parameters[split(s[j], '=')[0]] = split(s[j], '=')[1];

                }
                cout << "command: " << s[1] << endl;
                commands.push_back(command);
            }
            if (tolower(s[0]) == "x0_trajs") trajectory_params.x0_trajs = atof(s[1].c_str());
		}
	}


}

void CGrid::creategrid(int _nx, int _ny, double _dx, double _dy)
{
	omp_set_num_threads(12);
	cout<< "max_threads = " << omp_get_max_threads() <<endl;
	cout<< "num_threads = " << omp_get_num_threads() <<endl;

	#pragma omp parallel
      {
        #pragma omp single
        printf("num_threads = %d\n", omp_get_num_threads());
      }

	GP.nx = _nx;
	GP.ny = _ny;
	GP.dx = _dx;
	GP.dy = _dy;
	p.resize(GP.nx);
	for (int i = 0; i < GP.nx; i++)
	{
		p[i].resize(GP.ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K_gauss = CVector(2);
		}
	}

}

void CGrid::writeasmatrix(string filename, int component)
{
	ofstream file(filename);
	for (int i = 0; i < GP.nx; i++)
	{
		for (int j = 0; j < GP.ny; j++)
			file << p[i][j].V[component]<<",";
		file << endl;
	}


}

void CGrid::writeasmatrixK(string filename, int component)
{
	ofstream file(filename);
	for (int i = 0; i < GP.nx; i++)
	{
		for (int j = 0; j < GP.ny; j++)
			file << p[i][j].K[component] << ",";
		file << endl;
	}
}


CTimeSeries CGrid::getvelocity_gradient_samples()
{
    CTimeSeries out;
    for (int i=1; i<GP.nx; i++)
        for (int j=1; j<GP.ny; j++)
        {
            point p;
            p.x = i*GP.dx;
            p.y = j*GP.dy;
            out.append(i*GP.ny+j,getvelocity_gradient(p));
        }
    return out;
}

double CGrid::getvelocity_gradient(const point &pp, const string &direction)
{
    CVector out(2);
    int i_floar_x = int(pp.x / GP.dx);
	int j_floar_x = int(pp.y / GP.dy+0.5);
	int i_floar_y = int(pp.x / GP.dx+0.5);
	int j_floar_y = int(pp.y / GP.dy);
	if (pp.x<=0 || pp.x>=(GP.nx-1)*GP.dx || pp.y<=0 || pp.y>=GP.dy*(GP.ny-1))
        {   //qDebug() << "Empty CVector returned";
            return -999;

        }

    double vx1 = vx[i_floar_x][max(j_floar_x-1,0)] + 1.0/GP.dx*(pp.x - GP.dx*i_floar_x)*(vx[min(i_floar_x+1,GP.nx-1)][max(j_floar_x-1,0)]-vx[i_floar_x][max(j_floar_x-1,0)]);
    double vx2 = vx[i_floar_x][min(j_floar_x,GP.ny-2)] + 1.0/GP.dx*(pp.x - GP.dx*i_floar_x)*(vx[min(i_floar_x+1,GP.nx-1)][min(j_floar_x,GP.ny-2)]-vx[i_floar_x][min(j_floar_x,GP.ny-2)]);
    double vx_grad = (vx2-vx1)/GP.dy;

    return vx_grad;
}


CVector CGrid::getvelocity(point pp)
{
	int i_floar_x = int(pp.x / GP.dx);
	int j_floar_x = int(pp.y / GP.dy+0.5);
	int i_floar_y = int(pp.x / GP.dx+0.5);
	int j_floar_y = int(pp.y / GP.dy);
	if (pp.x<=0 || pp.x>=(GP.nx-1)*GP.dx || pp.y<=0 || pp.y>=GP.dy*(GP.ny-1))
        {   //qDebug() << "Empty CVector returned";
            return CVector();

        }

    double vx1 = vx[i_floar_x][max(j_floar_x-1,0)] + 1.0/GP.dx*(pp.x - GP.dx*i_floar_x)*(vx[min(i_floar_x+1,GP.nx-1)][max(j_floar_x-1,0)]-vx[i_floar_x][max(j_floar_x-1,0)]);
    double vx2 = vx[i_floar_x][min(j_floar_x,GP.ny-2)] + 1.0/GP.dx*(pp.x - GP.dx*i_floar_x)*(vx[min(i_floar_x+1,GP.nx-1)][min(j_floar_x,GP.ny-2)]-vx[i_floar_x][min(j_floar_x,GP.ny-2)]);
    double vx_interp = vx1 + (vx2-vx1)/GP.dy*(pp.y-GP.dy*(double(j_floar_x)-0.5));

    double vy1 = vy[max(i_floar_y-1,0)][j_floar_y] + 1.0/GP.dy*(pp.y - GP.dy*j_floar_y)*(vy[max(i_floar_y-1,0)][min(j_floar_y+1,GP.nx-1)]-vy[max(i_floar_y-1,0)][j_floar_y]);
    double vy2 = vy[min(i_floar_y,GP.nx-2)][j_floar_y] + 1.0/GP.dy*(pp.y - GP.dy*j_floar_y)*(vy[min(i_floar_y,GP.nx-2)][min(j_floar_y+1,GP.ny-1)]-vy[min(i_floar_y,GP.nx-2)][j_floar_y]);
    double vy_interp = vy1 + (vy2-vy1)/GP.dx*(pp.x-GP.dx*(double(i_floar_y)-0.5));

	CVector V(2);
	V[0] = vx_interp;
	V[1] = vy_interp;
	return V;
}


CVector CGrid::getvelocity_exact(point pp)
{
	int i_floar = int(pp.x / GP.dx);
	int j_floar = int(pp.y / GP.dy);
	if ((i_floar > GP.nx - 2) || (j_floar>GP.ny - 2) || (i_floar<0) || (j_floar<0))
        {   //qDebug() << "Empty CVector returned";
            return CVector();

        }
	CVector V1 = p[i_floar][j_floar].V + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar].V - p[i_floar][j_floar].V);
	CVector V2 = p[i_floar][j_floar+1].V + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar+1].V - p[i_floar][j_floar+1].V);
	CVector V = V1 + (1.0 / GP.dy*(pp.y - GP.dy*j_floar))*(V2 - V1);
	return V;
}

double CGrid::interpolate_K(double x, double y)
{
	point pp;
	pp.x = x;
	pp.y = y;
	int i_floar = int(pp.x / GP.dx);
	int j_floar = int(pp.y / GP.dy);
	if ((i_floar > GP.nx - 2) || (j_floar>GP.ny - 2) || (i_floar<0) || (j_floar<0)) return 0;
	CVector K1 = p[i_floar][j_floar].K + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar].K - p[i_floar][j_floar].K);
	CVector K2 = p[i_floar][j_floar + 1].K + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar + 1].K - p[i_floar][j_floar + 1].K);
	CVector K = K1 + (1.0 / GP.dx*(pp.y - GP.dy*j_floar))*(K2 - K1);
	return K[0];
}

vector<double> CGrid::interpolate_V(double x, double y)
{
	point pp;
	pp.x = x;
	pp.y = y;
	return getvelocity(pp).vec;
}

CPathway CGrid::gettrajectory(CPosition pp, double dt, double t_end)
{
	CPosition pt = pp;
	CPathway Trajectory;
        Trajectory.weight = getvelocity(pt)[0];
	double t = 0;

	while (t < t_end)
	{
            CVector V = getvelocity(pt);
            if (V.getsize() == 0) return Trajectory;
            CPosition ps;
            ps.x = pt.x;
            ps.y = pt.y;
            ps.t = t;
            ps.v = V;
            Trajectory.append(ps);
            pt.x += V[0] * dt;
            pt.y += V[1] * dt;
            t += dt;
	}

	return Trajectory;
}

CPathway CGrid::gettrajectory_vdt(CPosition pp, double dt0, double x_end, double tol, double diffusion)
{
	CPosition pt = pp;
	CPathway Trajectory;
        Trajectory.weight = getvelocity(pt)[0];
	double t = 0;
    double dt = dt0;

	while (pt.x < x_end)
	{
        CPosition ps;
        double err = tol;
        while (err>=tol)
        {
            dt /= (err/tol);
            CVector V = getvelocity(pt);
            if (V.getsize() == 0) return Trajectory;
            ps.x = pt.x;
            ps.y = pt.y;
            ps.t = t;
            ps.v = V;
            pt.x += V[0] * dt;
            pt.y += V[1] * dt;
            CVector V2 = getvelocity(pt);
            if (V2.num == 0)
            {
                Trajectory.append(ps);
                return Trajectory;
            }
            err = norm(V-V2);
            pt.x += -(V-V2)[0]*dt/2.0;
            pt.y += -(V-V2)[1]*dt/2.0;
        }
        t += dt;
        pt.x += gsl_ran_gaussian(Copula.RngPtr(),sqrt(D*dt));
        pt.y += gsl_ran_gaussian(Copula.RngPtr(),sqrt(D*dt));
        dt *= pow(1.05,1.0-err/tol);
        Trajectory.append(ps);
	}

	return Trajectory;
}


CPathway CGrid::gettrajectory_fix_dx(CPosition pp, double dx, double x_end)
{
    CPosition pt = pp;
    CPathway Trajectory;
    if (weighted)
        Trajectory.weight = getvelocity(pt)[0];
    else
        Trajectory.weight = 1;
    double t = 0;

    bool ex = false;
    while (pt.x < x_end && ex==false)
    {
        CVector V = getvelocity(pt);
        if (V.num == 2)
        {
            double dt = fabs(dx / sqrt(pow(V[0],2)+pow(V[1],2)));
            if (V.getsize() == 0) return Trajectory;
            CPosition ps;
            CVector V = getvelocity(pt);
            ps.x = pt.x;
            ps.y = pt.y;
            ps.t = t;
            ps.v = V;
            Trajectory.append(ps);
            pt.x += V[0] * dt;
            pt.y += V[1] * dt;
            t += dt;
        }
        else ex = true;
    }

    return Trajectory;
}


CBTCSet CGrid::get_correlation_based_on_random_samples(int nsamples, double dx0, double x_inc, bool magnitude)
{
    CBTCSet output(2);
    for (int i=0; i<nsamples; i++)
    {
        CPosition pt = getrandompoint();
        CVector V = v_correlation_single_point(pt,dx0,x_inc, magnitude);
        if (V.num==2)
        {
            output.BTC[0].append(i,V[0]);
            output.BTC[1].append(i,V[1]);
        }
    }
    return output;
}

CBTCSet CGrid::get_correlation_based_on_random_samples_dt(int nsamples, double dt0, double t_inc)
{
    CBTCSet output(2);
    for (int i=0; i<nsamples; i++)
    {
        CPosition pt = getrandompoint();
        CVector V = v_correlation_single_point_dt(pt,dt0,t_inc);
        if (V.num==2)
        {
            output.BTC[0].append(i,V[0]);
            output.BTC[1].append(i,V[1]);
        }
    }
    return output;
}

CBTCSet CGrid::get_correlation_based_on_random_samples_diffusion(int nsamples, double dt0, double diffusion_coeff, double increment, bool fixed_increment, const string &direction)
{
    CBTCSet output(2);
    for (int i=0; i<nsamples; i++)
    {
        CPosition pt = getrandompoint();
        CVector V = v_correlation_single_point_diffusion(pt,dt0,diffusion_coeff, increment, fixed_increment, direction);
        if (V.num==2)
        {
            output.BTC[0].append(i,V[0]);
            output.BTC[1].append(i,V[1]);
        }
    }
    return output;
}


CPosition CGrid::getrandompoint()
{
    CPosition out;
    out.x = unitrandom()*GP.dx*(GP.nx-1);
    out.y = unitrandom()*GP.dy*(GP.ny-1);
    return out;
}

CVector CGrid::v_correlation_single_point(const CPosition &pp, double dx0, double x_inc, bool magnitude)
{
    CPosition pt = pp;
    CPosition p_new;
    CVector Vout;
    double x_end = pt.x + dx0;
    bool ex = false;
    while (pt.x < x_end && ex==false)
    {
        CVector V = getvelocity(pt);
        if (V[0]<0)
        {
            cout<<"Negative Velocity! Quiting! "<<endl;
            return Vout;
        }

        if (V.getsize() == 0)
            {
                return Vout;
            }
        if (V[0]<0) return Vout;
        double dx = min(x_inc/(sqrt(pow(V[0],2)+pow(V[1],2)))*V[0],x_end-pt.x);

        bool changed_sign = true;
        while (changed_sign)
        {   if (V.num == 2)
            {
                if (V[0]==0)
                {
                    cout<<"Vx = 0!"<<endl;
                }
                double dt = dx/V[0];
                p_new.x = pt.x + dt*V[0];
                p_new.y = pt.y + dt*V[1];
                CVector V_new = getvelocity(p_new);

                if (V_new.getsize() == 0)
                {
                    return Vout;
                }
                if (V_new[0]<=0) return CVector();
                if (V_new[0]*V[0]<0)
                {
                    cout<<"V changed sign!"<<endl;
                    dx/=2;
                }
                else
                {   changed_sign = false;
                    if (V[0]==V_new[0])
                        dt = dx/V[0];
                    else
                        dt = dx*(log(V_new[0])-log(V[0]))/(V_new[0]-V[0]);
                    if (dt<0)
                    {
                        cout<<"Negative dt!"<<endl;
                    }
                    p_new.y = pt.y + 0.5*(V[1]+V_new[1])*dt;
                    p_new.x = pt.x + dx;
                    V_new = getvelocity(p_new);
                    if (V_new.num!=2)
                        return CVector();
                    p_new.v = V_new;
                    p_new.t += dt;
                    pt = p_new;

                }
            }
            else
            {
                return Vout;
            }
        }
    }
    Vout = CVector(2);
    if (!magnitude)
    {
        Vout[0] = getvelocity(pp)[0];
        Vout[1] = getvelocity(p_new)[0];
    }
    else
    {
        Vout[0] = sqrt(pow(getvelocity(pp)[0],2)+pow(getvelocity(pp)[1],2));
        Vout[1] = sqrt(pow(getvelocity(p_new)[0],2)+pow(getvelocity(p_new)[1],2));
    }

    return Vout;

}

CVector CGrid::v_correlation_single_point_dt(const CPosition &pp, double dt0, double t_inc)
{
    if (t_inc>dt0) t_inc = dt0/10.0;
    CPosition pt = pp;
    CPosition p_new;
    CVector Vout;
    bool ex = false;
    while (pt.t < dt0 && ex==false)
    {
        CVector V = getvelocity(pt);
        if (V.getsize() == 0)
            {
                return Vout;
            }
        double dx = V[0]*t_inc;
        double dy = V[1]*t_inc;

        if (V.num == 2)
        {
            p_new.x = pt.x + dx;
            p_new.y = pt.y + dy;
            CVector V_new = getvelocity(p_new);

            if (V_new.num!=2)
                return CVector();
            if (V_new[0]<=0) return CVector();
            p_new.y = pt.y + 0.5*(V[1]+V_new[1])*t_inc;
            p_new.x = pt.x + 0.5*(V[0]+V_new[0])*t_inc;
            V_new = getvelocity(p_new);
            if (V_new.num!=2)
                return CVector();
            p_new.t += t_inc;
            pt = p_new;
         }
    }
    Vout = CVector(2);
    Vout[0] = getvelocity(pp)[0];
    Vout[1] = getvelocity(p_new)[0];

    return Vout;

}

CVector CGrid::v_correlation_single_point_diffusion(const CPosition &pp, double dt0, double diffusion_coefficient, double increment, bool fixed_increment, const string &direction)
{
    if (increment==0) increment = dt0/10.0;
    CPosition pt = pp;
    CVector V1=getvelocity(pp);
    CPosition p_new;
    CVector Vout;
    if (V1[0]<=0) return Vout;
    if (V1.num!=2) return Vout;
    bool ex = false;
    p_new = pp;
    if (!fixed_increment)
    for (double dtt=0; dtt<dt0; dtt+=increment)
    {
        double zx = gsl_cdf_gaussian_Pinv(unitrandom(),1);
        double zy = gsl_cdf_gaussian_Pinv(unitrandom(),1);
        p_new.x += sqrt(2*increment*diffusion_coefficient)*zx;
        p_new.y += sqrt(2*increment*diffusion_coefficient)*zy;
    }
    else
    {
        double u;
        if (direction == "y")
            u = (int(unitrandom()-0.5) + 0.5)*3.141521;
        else if (direction == "x")
            u = int(unitrandom()-0.5)*3.141521;
        else
            u = unitrandom()*2*3.141521;

        double zx;
        double zy;
        if (direction=="g")
        {
            zx = dt0*gsl_cdf_gaussian_Pinv(unitrandom(),1);
            zy = dt0*gsl_cdf_gaussian_Pinv(unitrandom(),1);

        }
        else
        {
            zx = dt0*cos(u)*gsl_cdf_gaussian_Pinv(unitrandom(),1);
            zy = dt0*sin(u)*gsl_cdf_gaussian_Pinv(unitrandom(),1);
        }

        p_new.x += zx;
        p_new.y += zy;
    }

    CVector V_new = getvelocity(p_new);

    if (V_new.num!=2) return Vout;
    if (V_new[0]<=0) return Vout;
    Vout = CVector(2);
    Vout[0] = V1[0];
    Vout[1] = V_new[0];

    return Vout;

}



CPathway CGrid::gettrajectory_fix_dx_2nd_order(CPosition pp, double dx0, double x_end, double D, double tol, double weight)
{
    double x0 = pp.x;
    double backward = 1;
    if (x_end<x0)
        backward = -1;

    CPosition pt = pp;
    CPathway Trajectory;
    if (weighted)
        Trajectory.weight = getvelocity(pt)[0];
    else
        Trajectory.weight = 1;
    double t = 0;

    bool ex = false;
    int counter = 0;
    while ( backward*x0<=backward*pt.x && backward*pt.x<=backward*x_end && ex==false)
    {
        counter ++;
        CVector V = getvelocity(pt);
        if (V.getsize() == 0)
            {
                return Trajectory;
            }
        if (V[0]==0)
        {
            cout<< "Vx = 0!" <<endl;
        }

        double dx = fabs(dx0/(sqrt(pow(V[0],2)+pow(V[1],2)))*V[0]);

        bool changed_sign = true;
        while (changed_sign)
        {   if (V.num == 2)
            {
                if (V[0]==0)
                {
                    cout<<"Vx = 0!"<<endl;
                }

                double dt0 = fabs(dx/V[0]);
                if (dt0<0)
                {
                    cout<<"negative dt0!"<<endl;
                }

                CPosition p_new;
                p_new.x = pt.x + backward*dt0*V[0];//+diffusion[0];
                p_new.y = pt.y + backward*dt0*V[1];//+diffusion[1];
                CVector V_new = getvelocity(p_new);
                if (V_new.getsize() == 0)
                {
                    return Trajectory;
                }
                if (V_new[0]*V[0]<0)
                {
                    cout<<"V changed sign!"<<endl;
                    dx/=2;
                }
                else
                {   changed_sign = false;
                    double err = norm(V-V_new);
                    double dt1;
                    while (err>=tol)
                    {

                        if (V[0]==V_new[0])
                            dt1 = fabs(dx/V[0]);
                        else
                            dt1 = fabs(dx*(log(fabs(V_new[0]))-log(fabs(V[0])))/(fabs(V_new[0])-fabs(V[0])));

                        p_new.y = pt.y + 0.5*(V[1]+V_new[1])*dt1*backward;// + diffusion[1]*sqrt(dt1/dt0);
                        p_new.x = pt.x + dx*backward; // + diffusion[0]*sqrt(dt1/dt0);
                        CVector V_star = getvelocity(p_new);
                        if (V_star.getsize()!=2)
                            return Trajectory;
                        err = norm(V_star - V_new);
                        V_new = V_star;
                    }
                    p_new.v = V_new;
                    p_new.t = t+dt1;
                    CVector diffusion(2);

                    diffusion[0] = gsl_ran_gaussian(Copula.RngPtr(),sqrt(D*dt1));
                    diffusion[1] = gsl_ran_gaussian(Copula.RngPtr(),sqrt(D*dt1));
                    p_new.x = p_new.x + diffusion[0];
                    p_new.y = p_new.y + diffusion[1];

                    if (isnan(p_new.x))
                    {
                        cout<< "nan!" << endl;
                    }
                    Trajectory.append(p_new);
                    pt = p_new;

                    t += dt1;
                }
            }
            else ex = true;
        }
    }

    return Trajectory;
}



CPathwaySet CGrid::gettrajectories(double dt, double t_end)
{
    CPathwaySet X;
    for (int i = 0; i < int(pts.size()); i++)
    {
        cout << i << endl;
        CPathway X1 = gettrajectory(pts[i], dt, t_end);
        if (weighted)
        {   X.weighted = true;
            X.append(X1);
        }
        else
        {   X.weighted = false;
            X.append(X1);
        }
        set_progress_value(double(i) / double(pts.size()));
    }
    cout<<endl;
    return X;
}

CPathwaySet CGrid::gettrajectories_vdt(double dt, double t_end, double tol, double diffusion)
{

    CPathwaySet X;
    for (int i = 0; i < int(pts.size()); i++)
    {
        CPathway X1 = gettrajectory_vdt(pts[i], dt, t_end, tol, diffusion);
        if (weighted)
        {   X.weighted = true;
            X.append(X1);
        }
        else
        {   X.weighted = false;
            X.append(X1);
        }
        set_progress_value(double(i) / double(pts.size()));
    }
    pts.clear();
    for (int i=0; i<X.paths.size(); i++)
        pts.push_back(X.paths[i].positions[X.paths[i].positions.size()-1]);

    cout<<endl;
    return X;
}

CPathwaySet CGrid::gettrajectories_fixed_dx(double dx, double x_end, double tol, double diffusion)
{
    //qDebug() << "Simulating trajectories"<<endl;

    CPathwaySet X(pts.size());
    unsigned int counter = 0;
    #pragma omp parallel for
    for (int i = 0; i < int(pts.size()); i++)
    {
//	qDebug() << i << endl;

        CPathway X1 = gettrajectory_fix_dx_2nd_order(pts[i], dx, x_end, diffusion, tol);
        //cout << "\r" << "Trajectory #"<<  i << " Weight: " << X1.weight<< std::flush;
        if (weighted)
            {   X.weighted = true;
                X.paths[i] = X1;
            }
            else
            {   X.weighted = false;
                X.paths[i] = X1;
            }
        #pragma omp critical
        {
            counter++;
            set_progress_value(double(counter) / double(pts.size()));
        }
    }
    pts.clear();

    for (int i=0; i<X.paths.size(); i++)
        pts.push_back(X.paths[i].positions[X.paths[i].positions.size()-1]);
    cout<<endl;
    return X;
}


CBTC CGrid::initialize(int numpoints, double x_0, double smoothing_factor, bool flow_weighted, string boundary_dist_filename, bool _weighted)
{
    pts.clear();
    weighted = _weighted;
    cout << "weighted = " << weighted<<endl;
    CBTC vels;
    if (!_weighted)
    {
        cout << "weighted = " << weighted<<endl;
        int burnout = 0;
        CBTC boundary_v_dist = get_v_btc(x_0,0);
        double v_max = boundary_v_dist.maxC();
        if (boundary_dist_filename != "")
        {
            boundary_v_dist.distribution(40).writefile(boundary_dist_filename);
        }
        double y_0 = unitrandom()*GP.dy*(GP.ny-1);
        point pt_0; pt_0.x = x_0; pt_0.y = y_0;
        double v_x = getvelocity(pt_0)[0];
        pts.push_back(pt_0);
        for (int i = 1; i < numpoints; i++)
        {
            bool accepted = false;
            while (!accepted)
            {
                y_0 = unitrandom()*GP.dy*(GP.ny-1);
                pt_0.x = x_0; pt_0.y = y_0;
                v_x = getvelocity(pt_0)[0];
                double u = unitrandom();
                if (flow_weighted)
                    if (u < (v_x / v_max/5)) accepted = true;
                else
                    accepted = true;

            }
            pts.push_back(pt_0);
            vels.append(i,v_x);

            set_progress_value(double(i) / double(numpoints));
        }

    }
    else
    {   double y_0 = unitrandom()*GP.dy*(GP.ny-1);
        point pt_0; pt_0.x = x_0; pt_0.y = y_0;
        double v_x = getvelocity(pt_0)[0];
        pts.push_back(pt_0);
        for (int i = 1; i < numpoints; i++)
        {
            y_0 = unitrandom()*GP.dy*(GP.ny-1);
            pt_0.x = x_0; pt_0.y = y_0;
            double v_xp = v_x;
            v_x = getvelocity(pt_0)[0];
            double u = unitrandom();
            pts.push_back(pt_0);
            vels.append(i,v_x,v_x);
            set_progress_value(double(i) / double(numpoints));
        }
        cout<<endl;
    }
    return vels.distribution(40,(vels.maxC()-vels.minC())*smoothing_factor);
}

CMatrix_arma_sp CGrid::create_stiffness_matrix_arma()
{
	string averaging = "arithmetic";
//	qDebug() << "Creating stiffness matrix" << endl;
        CMatrix_arma_sp K((GP.nx+1)*(GP.ny+1), (GP.nx + 1)*(GP.ny + 1));
//	qDebug() << "Stiffness matrix created" << endl;
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) / (GP.dx*GP.dx);
			K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -(avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) + avg(p[i][j-1].K[0], p[i][j].K[0], averaging)) / (GP.dx*GP.dx);
			K.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = avg(p[i][j-1].K[0], p[i][j].K[0], averaging) / (GP.dx*GP.dx);

			K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) / (GP.dy*GP.dy);
			K.matr(get_cell_no(i, j), get_cell_no(i, j)) += -(avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) + avg(p[i-1][j].K[0], p[i][j].K[0], averaging)) / (GP.dy*GP.dy);
			K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = avg(p[i-1][j].K[0], p[i][j].K[0], averaging) / (GP.dy*GP.dy);
		}
		// top boundary
		int j = GP.ny;
		K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i+1, j)) = 1;

	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = 1;
		}

	return K;

}

CMatrix CGrid::create_stiffness_matrix()
{
	string averaging = "arithmetic";
	CMatrix K((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	CVector V((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			K[get_cell_no(i, j)][get_cell_no(i - 1, j)] = avg(p[i - 1][j - 1].K[0], p[i - 1][j].K[0], averaging) / (GP.dx*GP.dx);
			K[get_cell_no(i, j)][get_cell_no(i, j)] = -(avg(p[i - 1][j - 1].K[0], p[i - 1][j].K[0], averaging) + avg(p[i][j - 1].K[0], p[i][j].K[0], averaging)) / (GP.dx*GP.dx);
			K[get_cell_no(i, j)][get_cell_no(i + 1, j)] = avg(p[i][j - 1].K[0], p[i][j].K[0], averaging) / (GP.dx*GP.dx);

			K[get_cell_no(i, j)][get_cell_no(i, j - 1)] = avg(p[i - 1][j - 1].K[0], p[i][j - 1].K[0], averaging) / (GP.dy*GP.dy);
			K[get_cell_no(i, j)][get_cell_no(i, j)] += -(avg(p[i - 1][j - 1].K[0], p[i][j - 1].K[0], averaging) + avg(p[i - 1][j].K[0], p[i][j].K[0], averaging)) / (GP.dy*GP.dy);
			K[get_cell_no(i, j)][get_cell_no(i, j + 1)] = avg(p[i - 1][j].K[0], p[i][j].K[0], averaging) / (GP.dy*GP.dy);
		}
		// top boundary
		int j = GP.ny;
		K[get_cell_no(i, j)][get_cell_no(i, j - 1)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i, j)] = -1;

		// bottom boundary
		j = 0;
		K[get_cell_no(i, j)][get_cell_no(i, j + 1)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i, j)] = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K[get_cell_no(i, j)][get_cell_no(i, j)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i + 1, j)] = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K[get_cell_no(i, j)][get_cell_no(i, j)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i - 1, j)] = 1;
	}

	return K;

}


CVector_arma CGrid::create_RHS_arma()
{
	CVector_arma V((GP.nx + 1)*(GP.ny + 1));

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2*leftboundary_h;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2*rightboundary_h;


	return V;
}

CVector CGrid::create_RHS()
{
	CVector V((GP.nx + 1)*(GP.ny + 1));

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2 * leftboundary_h;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2 * rightboundary_h;


	return V;

}


int CGrid::get_cell_no(int i, int j)
{
	return i*(GP.ny+1) + j;
}

int CGrid::get_cell_no_OU(int i, int j)
{
	return i*(GP.ny) + j;
}

CMatrix CGrid::solve()
{
    set_progress_value(0);
    CMatrix_arma_sp K = create_stiffness_matrix_arma();
    CVector_arma V = create_RHS_arma();
 //   qDebug() << "Solving the system of equations"<< endl;
    CVector_arma S = solve_ar(K, V);
//    qDebug() << "Solved"<< endl;
    H = CMatrix(GP.nx+1,GP.ny+1);
    vx = CMatrix(GP.nx, GP.ny-1);
    vy = CMatrix(GP.nx - 1, GP.ny);
    for (int i = 0; i < GP.nx+1; i++)
        for (int j = 0; j < GP.ny+1; j++)
            H[i][j] = S[get_cell_no(i, j)];

    for (int i = 0; i < GP.nx; i++)
	for (int j = 0; j < GP.ny; j++)
	{
            double Kx1 = 0.5*(p[i][max(j - 1,0)].K[0] + p[i][j].K[0]);
            double Kx2 = 0.5*(p[i][j].K[0] + p[i][min(j+1,GP.ny-1)].K[0]);
            double Ky1 = 0.5*(p[max(i - 1, 0)][j].K[0] + p[i][j].K[0]);
            double Ky2 = 0.5*(p[i][j].K[0] + p[min(i+1,GP.nx-1)][j].K[0]);
            p[i][j].V[0] = -(0.5*Kx1 * (H[i + 1][j] - H[i][j]) / GP.dx + 0.5*Kx2 * (H[i + 1][j+1] - H[i][j+1]) / GP.dx);
            p[i][j].V[1] = -(0.5*Ky1 * (H[i][j+1] - H[i][j]) / GP.dy + 0.5*Ky2 * (H[i + 1][j + 1] - H[i+1][j]) / GP.dy);
            p[i][j].Vbx = -Kx1*(H[i+1][j] - H[i][j]) / GP.dx;
            p[i][j].Vtx = -Kx2*(H[i+1][j+1] - H[i][j+1]) / GP.dx;
            p[i][j].Vby = -Ky1*(H[i][j+1] - H[i][j]) / GP.dy;
            p[i][j].Vtx = -Ky2*(H[i+1][j+1] - H[i+1][j]) / GP.dy;
	}

	for (int i=0; i<GP.nx; i++)
		for (int j = 0; j < GP.ny - 1; j++)
		{
			double Kx = 0.5*(p[i][j].K[0] + p[i][j + 1].K[0]);
			vx[i][j] = -Kx * (H[i + 1][j + 1] - H[i][j + 1]) / GP.dx;
		}

        set_progress_value(0.5);
	for (int i = 0; i<GP.nx-1; i++)
		for (int j = 0; j < GP.ny; j++)
		{
			double Ky = 0.5*(p[i][j].K[0] + p[i+1][j].K[0]);
			vy[i][j] = -Ky * (H[i + 1][j+1] - H[i+1][j]) / GP.dy;
		}

        set_progress_value(1);

	max_v_x = max_vx();
	min_v_x = min_vx();
	return H;
    cout<<endl;
}

CTimeSeriesSet CGrid::get_Eulerian_vdist()
{
    CTimeSeriesSet out(2);
    for (int i=0; i<GP.nx; i++)
		for (int j = 0; j < GP.ny - 1; j++)
            out.BTC[0].append(i*10000+j,vx[i][j]);

        set_progress_value(0.5);
	for (int i = 0; i<GP.nx-1; i++)
		for (int j = 0; j < GP.ny; j++)
            out.BTC[1].append(i*10000+j,vy[i][j]);

        set_progress_value(1);


}

void CGrid::Assign_Linear_Velocity_Field(double V0, double V_slope)
{
    H = CMatrix(GP.nx+1,GP.ny+1);
    vx = CMatrix(GP.nx, GP.ny-1);
    vy = CMatrix(GP.nx - 1, GP.ny);
    for (int i=0; i<GP.nx; i++)
		for (int j = 0; j < GP.ny - 1; j++)
		{
			vx[i][j] = V0 + V_slope*(double(j)+0.5)*GP.dy;
		}


	for (int i = 0; i<GP.nx-1; i++)
		for (int j = 0; j < GP.ny; j++)
		{
			vy[i][j] = 0;
		}

	for (int i = 0; i < GP.nx; i++)
	for (int j = 0; j < GP.ny; j++)
	{
            p[i][j].K[0] = 1;
            p[i][j].V[0] = V0 + V_slope*(double(j))*GP.dy;
            p[i][j].V[1] = 0;
            p[i][j].Vbx = V0 + V_slope*(double(j)-0.5)*GP.dy;
            p[i][j].Vtx = V0 + V_slope*(double(j)+0.5)*GP.dy;
            p[i][j].Vby = 0;
            p[i][j].Vtx = 0;
	}

}

vector<int> CGrid::get_ij(int k)
{
	vector<int> out(2);
	out[0] = k / (GP.ny + 1);
	out[1] = k % (GP.ny + 1);
	return out;

}

CBTC CGrid::get_v_btc(int k)
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].V[k]);

	return out;
}

CBTC CGrid::get_v_dist_MODFlow(const string &filename)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    CBTC out;
    int rownum=0;
    getline(file);
    while (!file.eof())
    {
        vector<double> s1 = ATOF(getline(file,' '));
        if (s1.size()>=5)
            out.append(rownum,s1[3]);
        rownum++;
    }
    return out;


}

CBTC CGrid::get_v_dist_frac(const string &filename)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    CBTC out;
    int rownum=0;
    getline(file);
    while (!file.eof())
    {
        vector<double> s1 = ATOF(getline(file,' '));
        if (s1.size()>=5)
            out.append(rownum,s1[6]);
        rownum++;
    }
    return out;
}

CBTCSet CGrid::get_BTC_frac(const string &filename, const double &x_min, const double &x_max)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    CBTCSet out(1);
    out.setname(0,"time");
    //out.setname(1,"v");
    int rownum=0;
    getline(file);
    while (!file.eof())
    {
        vector<double> s1 = ATOF(getline(file,' '));
        if (s1.size()==2)
            set_progress_value("time= " + numbertostring(s1[1]));
        if (s1.size()>=5)
        {
            if (s1[2]<x_max && s1[2]>=x_min)
            {      out.BTC[0].append(s1[1],s1[1]);
                   //out.BTC[1].append(s1[1],s1[1]);
            }
        }
        rownum++;
    }
    return out;
}

CBTCSet CGrid::get_BTC_mf(const string &filename, const double &x_min, const double &x_max)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    CBTCSet out(1);
    out.setname(0,"time");
    //out.setname(1,"v");
    int rownum=0;
    getline(file);
    while (!file.eof())
    {
        vector<double> s1 = ATOF(getline(file,' '));
        if (s1.size()==2)
            set_progress_value("time= " + numbertostring(s1[1]));
        if (s1.size()>=5)
        {
            if (s1[1]<x_max && s1[1]>=x_min)
            {      out.BTC[0].append(s1[4],s1[4]);
                   //out.BTC[1].append(s1[1],s1[1]);
            }
        }
        rownum++;
    }
    return out;
}


CBTCSet CGrid::get_BTC_frac(CPathwaySet &pthwayset, const double &x_min, const double &x_max, bool consider_reaction)
{

    show_in_window("Getting Breakthrough curves from trajectories already loaded...");
    CBTCSet out(3);
    out.setname(0,"time");
    out.setname(1,"v");
    out.setname(2,"x");

    double mean_x = (x_min+x_max)*0.5;
    for (int i=0; i<pthwayset.paths.size(); i++)
    {
        for (int j=0; j<pthwayset.paths[i].size(); j++)
        {
            if (pthwayset.paths[i].positions[j].x<x_max && pthwayset.paths[i].positions[j].x>=x_min && (!consider_reaction || !pthwayset.paths[i].positions[j].reacted))
            {
                double t;
                if (pthwayset.paths[i].positions[j].v[0]!=0)
                    t = pthwayset.paths[i].positions[j].t - (pthwayset.paths[i].positions[j].x - mean_x)/pthwayset.paths[i].positions[j].v[0];
                else
                    t = pthwayset.paths[i].positions[j].t;
                out.BTC[0].append(t,t);
                out.BTC[1].append(t,pthwayset.paths[i].positions[j].v[0]);
                out.BTC[2].append(t,pthwayset.paths[i].positions[j].x);
            }
        }

    }
    return out;
}

CBTC CGrid::get_v_mag_btc()
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, sqrt(pow(p[i][j].V[0],2)+pow(p[i][j].V[1],2)));

	return out;
}

CBTC CGrid::get_kg_btc(int k)
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].K_gauss[k]);

	return out;
}

void CGrid::remap_K(int k)
{
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			p[i][j].K[k] = map_to_KCDF(getnormalcdf(p[i][j].K_gauss[k]));

}

CBTC CGrid::get_v_btc(double x, int k)
{
	CBTC out;
	int i = int(x / GP.dx);
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].V[k]);
	return out;

}

CBTC CGrid::get_v_dist(double x, int k, int nbins, double smoothing_factor)
{
    CBTC v = get_v_btc(x,k);
    return v.distribution(nbins,(v.maxC()-v.minC())*smoothing_factor);
}

void CGrid::set_inv_K_dist(int ninc)
{
	inv_K_dist.clear();
	double epsilon = 1e-8;

	for (int i = 0; i < ninc+1; i++)
	{
		double u = double(i) / double(ninc)*(1 - 2 * epsilon) + epsilon;
		if (i==0)
			inv_K_dist.append(u, inv_function_s(u));
		else
			inv_K_dist.append(u, inv_function_s(u));
	}

}

double CGrid::map_to_KCDF(double u)
{
	return inv_K_dist.interpol(u);
}

double CGrid::K_CDF(double x)
{
	double out = 0;
	if (marginal_K_dist_type == "lognormal")
	{
		int n = marginal_K_dist_params.size() / 3.0;

		for (int i = 0; i < n; i++)
			out += marginal_K_dist_params[3 * i] * 0.5*(1.0 + erf((log(x) - log(marginal_K_dist_params[3 * i + 1])) / (sqrt(2.0)*marginal_K_dist_params[3 * i + 2])));
	}
	return out;

}


double CGrid::inv_function(double u, double guess)
{
	double x = guess;
	double dx = 1e-8;
	double tol = 1e-6;
	double err = K_CDF(x) - u;
	double err_p;
	double lambda = 1;
	while (fabs(err) > tol)
	{
		double dfdx = (K_CDF(x + dx) - K_CDF(x)) / dx;
		x -= lambda*err / dfdx;
		err_p = err;
		err = K_CDF(x) - u;
		if (x < 0 || !(x == x) || fabs(err)>fabs(err_p))
		{
			x += lambda*err_p / dfdx;
			lambda /= 2;
			err = err_p;
		}
		else lambda = min(lambda*1.3, 1.0);

	}
	return x;
}

double CGrid::inv_function_s(double u)
{
	double x2 = 1000;
	double x1 = 0.001;
	double tol = 1e-8;
	double err1 = K_CDF(x1) - u;
	double err2 = K_CDF(x2) - u;
	while (err1 > 0)
	{
		x1 /= 2;
		err1 = K_CDF(x1) - u;
	}

	while (err2 < 0)
	{
		x2 *= 2;
		err2 = K_CDF(x2) - u;
	}

	while (min(fabs(err1),fabs(err2)) > tol && fabs(x1-x2)>tol)
	{
		double slope = (err2 - err1) / (log(x2) - log(x1));
		double x_p = exp(log(x1) - err1 / slope);
		double err_p = K_CDF(x_p) - u;
		if (err_p > 0)
		{
			x2 = x_p;
			err2 = K_CDF(x2) - u;
		}
		else
		{
			x1 = x_p;
			err1 = K_CDF(x1) - u;
		}

	}
	if (fabs(err1) > fabs(err2))
		return x2;
	else
		return x1;
}

CBTC CGrid::get_K_CDF(double x0, double x1, double log_inc)
{
	CBTC out;
	for (double logx = log(x0); logx <= log(x1); logx += log_inc )
	{
		double x = exp(logx);
		out.append(x, K_CDF(x));
	}

	return out;
}

CBTC CGrid::get_V_PDF(double x0, double x1, double log_inc, bool _log)
{
	CBTC out;
	if (_log)
    {   for (double logx = log(x0); logx <= log(x1); logx += log_inc )
        {
            double x = exp(logx);
            out.append(x, dist.evaluate(x));
        }
    }
    else
    {
        for (double x = x0; x <= x1; x += log_inc )
        {
            out.append(x, dist.evaluate(x));
        }
    }

	return out;
}

CBTC CGrid::get_margina_traj_v_dist(double vmin, double vmax, double nbins, string val)
{
	CBTC out(nbins + 1);
	int count = 0;
	for (int i = 0; i < int(Traj.paths.size()); i++)
	{
		for (int j = 0; j < int(Traj.paths[i].positions.size()); j++)
			out.C[int((Traj.paths[i].positions[j].getvar(val) - vmin) / nbins)] += 1;
		count += 1;
	}
	out = out/count;
	return out;
}


double avg(double x, double y, string type)
{
	if (type == "arithmetic")
            return 0.5*(x + y);
	if (type == "geometric")
            return sqrt(x*y);
	if (type == "harmonic")
            return (2.0*x*y / (x + y));
        else
            return 0.5*(x + y);

}

vector<ijval> get_top_n(vector<ijval> vec, int n)
{
	vector<ijval> out;
	vector<bool> extracted(vec.size());
	for (int i = 0; i < int(vec.size()); i++) extracted[i] = false;
	int smallest_dist = -1;

	for (int i = 0; i < n; i++)
	{
		double min_dist = 1e12;
		for (int j = 0; j < int(vec.size()); j++)
			if ((vec[j].val < min_dist) && (extracted[j]==false))
			{
				smallest_dist = j;
				min_dist = vec[j].val;
			}
		out.push_back(vec[smallest_dist]);
		extracted[smallest_dist] = true;
	}

	return out;
}





double CGrid::max_K()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].K[0]);

	return max_k;

}
double CGrid::min_K()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].K[0]);
	return min_k;
}
double CGrid::max_vx()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].V[0]);
	return max_k;
}

double CGrid::min_vx()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].V[0]);
	return min_k;
}
double CGrid::max_vy()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].V[1]);
	return max_k;

}
double CGrid::min_vy()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].V[1]);
	return min_k;
}

void CGrid::set_K_transport(double dt, double D, double weight)
{
	string averaging = "arithmetic";
	Kv = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	KD = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	Kt = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));

	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			Kv.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*pos(vx[i-1][j-1]) / (GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(neg(vx[i-1][j-1])/GP.dx + pos(vx[i][j-1])/GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*neg(vx[i][j - 1]) / (GP.dx);

			Kv.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*pos(vy[i - 1][j - 1]) / GP.dy;
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i-1][j]) / GP.dy);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -weight*neg(vy[i-1][j]) / GP.dy;

			KD.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*D / (GP.dx*GP.dx);

			KD.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -weight*D / (GP.dy*GP.dy);

			Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) += 1.0 / dt;

		}
		// top boundary
		int j = GP.ny;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = -1;
	}


}

CVector_arma CGrid::create_RHS_transport(int species_counter, double dt, double weight,double D, double decay_coefficient, double decay_order)
{
	CVector_arma RHS((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			double rhs = 0;
			rhs += (1-weight)*pos(vx[i - 1][j - 1]) / (GP.dx)*C[species_counter][i-1][j];
			rhs += -(1-weight)*(neg(vx[i - 1][j - 1]) / GP.dx + pos(vx[i][j-1]) / GP.dx)*C[species_counter][i][j];
			rhs += (1-weight)*neg(vx[i][j - 1]) / (GP.dx)*C[species_counter][i+1][j];

			rhs += (1-weight)*pos(vy[i - 1][j - 1]) / GP.dy*C[species_counter][i][j-1];
			rhs += -(1 - weight)*(neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i - 1][j]) / GP.dy)*C[species_counter][i][j];
			rhs += (1-weight)*neg(vy[i - 1][j]) / GP.dy*C[species_counter][i][j+1];

			rhs += (1-weight)*D / (GP.dx*GP.dx)*C[species_counter][i-1][j];
			rhs += -2 * (1-weight)*D / (GP.dx*GP.dx)*C[species_counter][i][j];
			rhs += (1-weight)*D / (GP.dx*GP.dx)*C[species_counter][i + 1][j];

			rhs += (1-weight)*D / (GP.dy*GP.dy)*C[species_counter][i][j-1];
			rhs += -2 * (1-weight)*D / (GP.dy*GP.dy)*C[species_counter][i][j];
			rhs += (1-weight)*D / (GP.dy*GP.dy)*C[species_counter][i][j+1];

			if (numberofspecies==1)
                rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*pow(C[species_counter][i][j], decay_order);
            else if (numberofspecies==2)
                rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*C[0][i][j]*C[1][i][j];
            else if (numberofspecies==3)
            {
                if (species_counter<2)
                {
                    rhs += 1.0 / dt*C[species_counter][i][j] - decay_coefficient*C[0][i][j]*C[1][i][j];
                }
                else
                    rhs += 1.0 / dt*C[species_counter][i][j] + decay_coefficient*C[0][i][j]*C[1][i][j];
            }
			RHS[get_cell_no(i, j)] = rhs;
		}
		// top boundary
		int j = GP.ny;
		RHS[get_cell_no(i, j)] = 0;


		// bottom boundary
		j = 0;
		RHS[get_cell_no(i, j)] = 0;

	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 2*leftboundary_C[species_counter];


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 0;

	return RHS;

}

CVector_arma CGrid::create_RHS_transport_laplace(int species_counter, double weight, double D, double s)
{
	CVector_arma RHS((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		// top boundary
		int j = GP.ny;
		RHS[get_cell_no(i, j)] = 0;


		// bottom boundary
		j = 0;
		RHS[get_cell_no(i, j)] = 0;

	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 2 * leftboundary_C[species_counter];


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 0;

	return RHS;

}

void CGrid::solve_transport(double t_end, vector<double> decay_coeff, vector<double> decay_order)
{
	set_K_transport(dt, D, time_weight);
	C.resize(numberofspecies);
	for (int species_counter=0; species_counter<numberofspecies; species_counter++)
	{
        C[species_counter] = CMatrix(GP.nx + 1, GP.ny + 1);// = leftboundary_C;
    }
    for (int i = 0; i < GP.nx; i++)
        for (int j = 0; j < GP.ny; j++)
            p[i][j].C.resize(int(t_end/dt), numberofspecies);
	CMatrix_arma_sp K = KD + Kt + Kv;
	//Kt.writetofile(pathout + "Kt_matrix.txt");
	//KD.writetofile(pathout + "KD_matrix.txt");
	//Kv.writetofile(pathout + "Kv_matrix.txt");
	//K.writetofile(pathout + "transport_matrix.txt");
	set_progress_value(0);
        int counter=0;
        for (double t = 0; t < t_end; t += dt)
        {
            for (int species_counter=0; species_counter<numberofspecies; species_counter++)
            {   CVector_arma RHS = create_RHS_transport(species_counter, dt, time_weight, D, decay_coeff[species_counter], decay_order[species_counter]);
                CVector_arma S = solve_ar(K, RHS);

                for (int i=0; i<GP.nx+1; i++)
                    for (int j=0; j<GP.ny+1; j++)
                        C[species_counter][i][j] = S[get_cell_no(i, j)];
            }
            for (int i = 0; i < GP.nx; i++)
                for (int j = 0; j < GP.ny; j++)
                {
                    vector<double> cc;
                    for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
                        p[i][j].C.setvalue(counter,species_counter,C[species_counter][i+1][j]);
                }
            counter++;
            set_progress_value(t / t_end);

        }
        cout<<endl;
}

void CGrid::solve_transport_laplace(double s)
{
	set_K_transport_laplace(D, s);
	C.resize(numberofspecies);
	for (int species_counter=0; species_counter<numberofspecies; species_counter++)
	{
        C[species_counter] = CMatrix(GP.nx + 1, GP.ny + 1);// = leftboundary_C;
    }

	CMatrix_arma_sp K = KD + Kt + Kv;
    for (int species_counter=0; species_counter<numberofspecies; species_counter++)
    {
        CVector_arma RHS = create_RHS_transport_laplace(species_counter, dt, time_weight, D);
        CVector_arma S = solve_ar(K, RHS);

        for (int i = 0; i<GP.nx + 1; i++)
            for (int j = 0; j<GP.ny + 1; j++)
                C[species_counter][i][j] = S[get_cell_no(i, j)];

        for (int i = 0; i < GP.nx; i++)
            for (int j = 0; j < GP.ny; j++)
                p[i][j].C.setvalue(i,j,0.25*(C[species_counter][i][j] + C[species_counter][i + 1][j] + C[species_counter][i][j + 1] + C[species_counter][i + 1][j + 1]));

    }


}


void CGrid::set_K_transport_laplace(double D, double s)
{
	string averaging = "arithmetic";
	Kv = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	KD = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	Kt = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));

	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			Kv.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -pos(vx[i - 1][j - 1]) / (GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += (neg(vx[i - 1][j - 1]) / GP.dx + pos(vx[i][j - 1]) / GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -neg(vx[i][j - 1]) / (GP.dx);

			Kv.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -pos(vy[i - 1][j - 1]) / GP.dy;
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += (neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i - 1][j]) / GP.dy);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -neg(vy[i - 1][j]) / GP.dy;

			KD.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2 * D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -D / (GP.dx*GP.dx);

			KD.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2 * D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -D / (GP.dy*GP.dy);

			Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) += s;

		}
		// top boundary
		int j = GP.ny;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = -1;
	}


}


void CGrid::create_inverse_K_OU(double dt)
{
	CMatrix_arma M(GP.ny*(GP.nx+2), GP.ny*(GP.nx+2));

	for (int i = 1; i < GP.nx+1; i++)
	{
		for (int j = 0; j < GP.ny; j++)
		{
			// Advection
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1.0 / dt;
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.FinvU[j] / GP.dx;
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*OU.FinvU[j] / GP.dx;

			/*if (i < GP.nx - 1)
			{
				M(get_cell_no(i, j), get_cell_no(i, j)) = 2 * time_weight*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
				M(get_cell_no(i, j), get_cell_no(i-1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
				M(get_cell_no(i, j), get_cell_no(i+1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			}*/

			//Diffusion
			if (i < GP.nx+ 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += 2 * time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i + 1, j)) += -time_weight*D / pow(GP.dx, 2);
			}
			else
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*D / pow(GP.dx, 2);
			}
			//Exchange
			if (j > 0 && j < GP.ny - 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*(OU.Exchanges[j-1]+OU.Exchanges[j])/pow(GP.dy,2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j - 1)) -= time_weight*OU.Exchanges[j-1] / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j + 1)) -= time_weight*OU.Exchanges[j] / pow(GP.dy, 2);
			}
			else if (j == 0)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.Exchanges[j] / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j + 1)) -= time_weight*OU.Exchanges[j] / pow(GP.dy, 2);
			}
			else if (j == GP.ny - 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.Exchanges[j-1] / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j - 1)) -= time_weight*OU.Exchanges[j-1] / pow(GP.dy, 2);

			}
		}

	}

	int i = 0;
	for (int j = 0; j < GP.ny; j++)
	{
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1;
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i+1, j)) = 1;
	}

	i = GP.nx+1;
	for (int j = 0; j < GP.ny; j++)
	{
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1;
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) = -1;
	}

	OU.Inv_M = inv(M);
	CMatrix MM_inv(OU.Inv_M);

}

void CGrid::create_inv_K_Copula(double dt, double Diffusion_coefficient)
{
    CMatrix_arma_sp M(GP.ny*(GP.nx+2), GP.ny*(GP.nx+2));

    for (int i = 1; i < GP.nx+1; i++)
    {
        for (int j = 0; j < GP.ny; j++)
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1.0 / dt;
            if (OU.FinvU[j]>=0)
            {   M.matr(j + GP.ny*i,j + GP.ny*i) += time_weight*OU.FinvU[j] / GP.dx;
                M.matr(j + GP.ny*i,j + GP.ny*(i - 1)) += -time_weight*OU.FinvU[j] / GP.dx;
            }
            else
            {
                M.matr(j + GP.ny*i,j + GP.ny*i) += -time_weight*OU.FinvU[j] / GP.dx;
                M.matr(j + GP.ny*i,j + GP.ny*(i + 1)) += time_weight*OU.FinvU[j] / GP.dx;
            }
            M.matr(j + GP.ny*i,j + GP.ny*i) += 2*time_weight*Diffusion_coefficient/pow(GP.dx,2);
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) -= time_weight*Diffusion_coefficient/pow(GP.dx,2);
            M.matr(j + GP.ny*i,j + GP.ny*min(i+1,GP.nx)) -= time_weight*Diffusion_coefficient/pow(GP.dx,2);

            for (int k = 0; k < GP.ny; k++)
                M.matr(j + GP.ny*i,k + GP.ny*i) += -time_weight*copula_params.K[j][k] / copula_params.epsilon*GP.dy;

        }
    }
    int i = 0;
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>=0)
        {   M.matr(j + GP.ny*i,j + GP.ny*i) = 0.5;
            M.matr(j + GP.ny*i,j + GP.ny*(i+1)) = 0.5;
        }
        else
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i+1)) = -1;
        }

    }

    i = GP.nx + 1;
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>=0)
        {   M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) = -1;
        }
        else
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) = -1;
        }
    }


    //Changed from inverse
    copula_params.Inv_M = M;


}

void CGrid::create_inv_K_Copula_diffusion(double dt, double Diffusion_coefficient)
{
    CMatrix_arma_sp M(GP.ny*(GP.nx+2), GP.ny*(GP.nx+2));

    for (int i = 1; i < GP.nx+1; i++)
    {
        for (int j = 0; j < GP.ny; j++)
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1.0 / dt;
            if (OU.FinvU[j]>=0)
            {   M.matr(j + GP.ny*i,j + GP.ny*i) += time_weight*OU.FinvU[j] / GP.dx;
                M.matr(j + GP.ny*i,j + GP.ny*(i - 1)) += -time_weight*OU.FinvU[j] / GP.dx;
            }
            else
            {
                M.matr(j + GP.ny*i,j + GP.ny*i) += -time_weight*OU.FinvU[j] / GP.dx;
                M.matr(j + GP.ny*i,j + GP.ny*(i + 1)) += time_weight*OU.FinvU[j] / GP.dx;
            }
            M.matr(j + GP.ny*i,j + GP.ny*i) += 2*time_weight*Diffusion_coefficient/pow(GP.dx,2);
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) -= time_weight*Diffusion_coefficient/pow(GP.dx,2);
            M.matr(j + GP.ny*i,j + GP.ny*min(i+1,GP.nx)) -= time_weight*Diffusion_coefficient/pow(GP.dx,2);

            for (int k = 0; k < GP.ny; k++)
                M.matr(j + GP.ny*i,k + GP.ny*i) += -time_weight*(copula_params.K_disp[j][k] / copula_params.epsilon+2*Diffusion_coefficient*copula_params.K_diff[j][k] / pow(copula_params.tau,2))*GP.dy;

        }
    }
    int i = 0;
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>=0)
        {   M.matr(j + GP.ny*i,j + GP.ny*i) = 0.5;
            M.matr(j + GP.ny*i,j + GP.ny*(i+1)) = 0.5;
        }
        else
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i+1)) = -1;
        }

    }

    i = GP.nx + 1;
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>=0)
        {   M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) = -1;
        }
        else
        {
            M.matr(j + GP.ny*i,j + GP.ny*i) = 1;
            M.matr(j + GP.ny*i,j + GP.ny*(i-1)) = -1;
        }
    }

    //Changed from inverse
    copula_params.Inv_M = M;


}



void CGrid::create_k_mat_copula()
{
    copula_params.K = CMatrix(GP.ny);
    #ifdef symetric
    for (int i = 0; i < GP.ny; i++)
        for (int j = 0; j < GP.ny; j++)
            if (i != j)
            {
                double u1 = double(i)*GP.dy + GP.dy / 2;
                double u2 = double(j)*GP.dy + GP.dy / 2;
                copula_params.K[i][j] = Copula.evaluate11(u1, u2)*mean(dist.inverseCDF(u1),dist.inverseCDF(u2));
                copula_params.K[i][i] -= copula_params.K[i][j];

            }
#else
    for (int i = 0; i < GP.ny; i++)
        for (int j = 0; j < GP.ny; j++)
            if (i != j)
            {
                double u1 = double(i)*GP.dy + GP.dy / 2;
                double u2 = double(j)*GP.dy + GP.dy / 2;
                copula_params.K[i][j] = Copula.evaluate11(u1, u2)*(fabs(mean(dist.inverseCDF(u2),dist.inverseCDF(u1))) + 2.0*Copula.correlation_ls*Copula.diffusion_coeff/pow(Copula.diffusion_correlation_ls,2));
                copula_params.K[i][i] -= Copula.evaluate11(u1, u2)*(fabs(mean(dist.inverseCDF(u1),dist.inverseCDF(u2))) + 2.0*Copula.correlation_ls*Copula.diffusion_coeff/pow(Copula.diffusion_correlation_ls,2));

            }
#endif
}


void CGrid::create_k_mat_copula_only_dispersion()
{
    copula_params.K_disp = CMatrix(GP.ny);
    for (int i = 0; i < GP.ny; i++)
    for (int j = 0; j < GP.ny; j++)
        if (i != j)
        {
            double u1 = double(i)*GP.dy + GP.dy / 2;
            double u2 = double(j)*GP.dy + GP.dy / 2;
            copula_params.K_disp[i][j] = Copula.evaluate11(u1, u2)*(fabs(mean(dist.inverseCDF(u2),dist.inverseCDF(u1))));
            copula_params.K_disp[i][i] -= Copula.evaluate11(u1, u2)*(fabs(mean(dist.inverseCDF(u1),dist.inverseCDF(u2))));

        }

}


void CGrid::create_k_mat_copula_only_diffusion()
{
    copula_params.K_diff = CMatrix(GP.ny);
    for (int i = 0; i < GP.ny; i++)
    for (int j = 0; j < GP.ny; j++)
        if (i != j)
        {
            double u1 = double(i)*GP.dy + GP.dy / 2;
            double u2 = double(j)*GP.dy + GP.dy / 2;
            copula_params.K_diff[i][j] = Copula_diffusion.evaluate11(u1, u2);
            copula_params.K_diff[i][i] -= Copula_diffusion.evaluate11(u1, u2);

        }

}


void CGrid::create_f_inv_u()
{
	OU.FinvU = CVector(GP.ny);
	for (int j = 0; j < GP.ny; j++)
	{
		double u1 = double(j)*GP.dy + GP.dy / 2;
		OU.FinvU[j] = dist.inverseCDF(u1);
	}
}

void CGrid::create_ou_exchange()
{
	OU.Exchanges = CVector(GP.ny-1);
	for (int j = 0; j < GP.ny-1; j++)
	{
        double v = (OU.FinvU[j] + OU.FinvU[j+1])/2.0;
        OU.Exchanges[j] = (v/OU.lc + 2.0*OU.diffusion/pow(OU.ld,2))*pow(std_normal_phi_inv((double(j+1))*GP.dy), 2)  ;
	}
}


CVector_arma CGrid::create_RHS_Copula(int species_counter, double dt, double diffusion_coeff, double decay_coeff, double decay_order)
{
    CVector_arma RHS(GP.ny*(GP.nx+2));

    for (int i = 1; i < GP.nx+1; i++)
    {
        for (int j = 0; j < GP.ny; j++)
        {
            if (numberofspecies==1)
                RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*pow(C[species_counter][i][j], decay_order);
            else if (numberofspecies==2)
                RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*C[0][i][j]*C[1][i][j];
            else if (numberofspecies==3)
            {
                if (species_counter<2)
                {
                    RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*C[0][i][j]*C[1][i][j];
                }
                else
                    RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] + decay_coeff*C[0][i][j]*C[1][i][j];
            }

            if (OU.FinvU[j]>0)
            {
                RHS[j + GP.ny*i] -= (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i][j];
                RHS[j + GP.ny*i] += (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i - 1][j];
            }
            else
            {
                RHS[j + GP.ny*i] += (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i][j];
                RHS[j + GP.ny*i] -= (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i + 1][j];
            }
            RHS[j + GP.ny*i] -= 2*(1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][i][j];
            RHS[j + GP.ny*i] += (1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][i - 1][j];
            RHS[j + GP.ny*i] += (1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][min(i + 1,GP.nx)][j];

            for (int k = 0; k < GP.ny; k++)
                RHS[j + GP.ny*i] += (1 - time_weight)*copula_params.K[j][k] / copula_params.epsilon*GP.dy*C[species_counter][i][k];

        }
    }
    int i = 0;
    #ifdef symetrical
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>0)
            RHS[j + GP.ny*i] = leftboundary_C[species_counter];
        else
            RHS[j + GP.ny*i] = 0;

    }
    #else
    double sum = 1;
    for (int j = 0; j < GP.ny; j++)
        sum += OU.FinvU[j]*GP.dy;
    for (int j = 0; j < GP.ny; j++)
        RHS[j + GP.ny*i] = 1.0/OU.FinvU[j]*sum;
    #endif

    i = GP.nx + 1;
    for (int j = 0; j < GP.ny; j++)
    {   if (OU.FinvU[j]>0)
            RHS[j + GP.ny*i] = 0;
        else
            RHS[j + GP.ny*i] = 0;
    }

    return RHS;
}


CVector_arma CGrid::create_RHS_Copula_diffusion(int species_counter, double dt, double diffusion_coeff, double decay_coeff, double decay_order)
{
    CVector_arma RHS(GP.ny*(GP.nx+2));

    for (int i = 1; i < GP.nx+1; i++)
    {
        for (int j = 0; j < GP.ny; j++)
        {
            if (numberofspecies==1)
                RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*pow(C[species_counter][i][j], decay_order);
            else if (numberofspecies==2)
                RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*C[0][i][j]*C[1][i][j];
            else if (numberofspecies==3)
            {
                if (species_counter<2)
                {
                    RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] - decay_coeff*C[0][i][j]*C[1][i][j];
                }
                else
                    RHS[j + GP.ny*i] += 1.0 / dt*C[species_counter][i][j] + decay_coeff*C[0][i][j]*C[1][i][j];
            }

            if (OU.FinvU[j]>0)
            {
                RHS[j + GP.ny*i] -= (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i][j];
                RHS[j + GP.ny*i] += (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i - 1][j];
            }
            else
            {
                RHS[j + GP.ny*i] += (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i][j];
                RHS[j + GP.ny*i] -= (1 - time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i + 1][j];
            }
            RHS[j + GP.ny*i] -= 2*(1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][i][j];
            RHS[j + GP.ny*i] += (1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][i - 1][j];
            RHS[j + GP.ny*i] += (1 - time_weight)*diffusion_coeff / pow(GP.dx,2)*C[species_counter][min(i + 1,GP.nx)][j];

            for (int k = 0; k < GP.ny; k++)
                RHS[j + GP.ny*i] += (1 - time_weight)*(copula_params.K_disp[j][k] / copula_params.epsilon + copula_params.K_diff[j][k] / copula_params.tau)*GP.dy*C[species_counter][i][k];

        }
    }
    int i = 0;
    for (int j = 0; j < GP.ny; j++)
    {
        if (OU.FinvU[j]>0)
            RHS[j + GP.ny*i] = 1;
        else
            RHS[j + GP.ny*i] = 0;

    }

    i = GP.nx + 1;
    for (int j = 0; j < GP.ny; j++)
    {   if (OU.FinvU[j]>0)
            RHS[j + GP.ny*i] = 0;
        else
            RHS[j + GP.ny*i] = 1;
    }

    return RHS;
}



CVector_arma CGrid::create_RHS_OU(int species_counter, double dt, double decay_coeff, double decay_order)
{
	CVector_arma RHS(GP.ny*(GP.nx+2));

	for (int i = 1; i < GP.nx+1; i++)
	{
		for (int j = 0; j < GP.ny; j++)
		{
			// Advection
			RHS[get_cell_no_OU(i, j)] = 1.0 / dt*C[species_counter][i][j] - decay_coeff*pow(C[species_counter][i][j],decay_order);
			RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i][j];
			RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.FinvU[j] / GP.dx*C[species_counter][i-1][j];

			/*if (i < GP.nx - 1)
			{
			M(get_cell_no(i, j), get_cell_no(i, j)) = 2 * time_weight*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			M(get_cell_no(i, j), get_cell_no(i-1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			M(get_cell_no(i, j), get_cell_no(i+1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			}*/

			//Diffusion
			if (i < GP.nx - 1)
			{
				RHS[get_cell_no_OU(i, j)] -= 2 * (1-time_weight)*D / pow(GP.dx, 2)*C[species_counter][i][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[species_counter][i-1][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[species_counter][i+1][j];
			}
			else
			{
				RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*D / pow(GP.dx, 2)*C[species_counter][i][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[species_counter][i-1][j];
			}
			//Exchange
			if (j > 0 && j < GP.ny - 1)
			{
				RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*(OU.Exchanges[j-1]+OU.Exchanges[j]) *C[species_counter][i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.Exchanges[j-1]*C[species_counter][i][j-1] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.Exchanges[j]*C[species_counter][i][j+1] / pow(GP.dy, 2);
			}
			else if (j == 0)
			{
				RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*OU.Exchanges[j]*C[species_counter][i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.Exchanges[j]*C[species_counter][i][j + 1] / pow(GP.dy, 2);
			}
			else if (j == GP.ny - 1)
			{
				RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*OU.Exchanges[j-1]*C[species_counter][i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.Exchanges[j-1]*C[species_counter][i][j - 1] / pow(GP.dy, 2);
			}
		}

	}

	int i = 0;
	for (int j = 0; j < GP.ny; j++)
	{
		RHS[get_cell_no_OU(i, j)] = 2;
	}

	i = GP.nx+1;
	for (int j = 0; j < GP.ny; j++)
	{
		RHS[get_cell_no_OU(i, j)] = 0;
	}
	return RHS;

}

void CGrid::solve_transport_OU(double t_end, double decay_coeff, double decay_order)
{
	create_f_inv_u();
	create_ou_exchange();
	create_inverse_K_OU(dt);
	C.resize(numberofspecies);
	for (int species_counter=0; species_counter<numberofspecies; species_counter++)
        C[species_counter] = CMatrix(GP.nx+2, GP.ny);
	OU.BTCs = CBTCSet(GP.nx+2);
	OU.BTC_normal = CBTCSet(GP.nx + 2);
	OU.BTCs_fw = CBTCSet(GP.nx + 2);
	OU.BTC_normal_fw = CBTCSet(GP.nx + 2);
	for (int i = 0; i < GP.nx+2; i++) OU.BTCs.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	//K.writetofile(pathout + "transport_matrix.txt");
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs.BTC[i].append(0, 0);
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.BTC[i].append(0, 0);
	set_progress_value(0);
    for (double t = 0; t < t_end; t += dt)
    {
        for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
        {
            CVector_arma RHS = create_RHS_OU(species_counter, dt,decay_coeff, decay_order);
            CVector_arma S = OU.Inv_M*RHS;

            for (int i = 0; i < GP.nx+2; i++)
            {
                double sum = 0;
                double sum_fw = 0;

                for (int j = 0; j < GP.ny; j++)
                {
                    C[species_counter][i][j] = S[get_cell_no_OU(i, j)];
                    sum += C[species_counter][i][j] * GP.dy;
                    sum_fw += C[species_counter][i][j] * OU.FinvU[j]/OU.FinvU.sum();
                }
                OU.BTCs.BTC[i].append(t+dt, sum);
                OU.BTCs_fw.BTC[i].append(t + dt, sum_fw);
                OU.BTC_normal.BTC[i].append((t + dt)/((i-0.5)*GP.dx), sum);
                OU.BTC_normal_fw.BTC[i].append((t + dt) / ((i - 0.5)*GP.dx), sum_fw);
            }


        }

        for (int i = 0; i < GP.nx; i++)
            for (int j = 0; j < GP.ny; j++)
            {
                for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
                    p[i][j].C.setvalue(i,j,C[species_counter][i+1][j]);

            }
        #if QT_version
                set_progress_value(t / t_end);
        tbrowse->append("t = " + QString::number(t));
        #else
        cout<<"t = " <<t <<endl;
        #endif // QT_version

    }

}

void CGrid::solve_transport_Copula(double t_end, double Diffusion_coeff, vector<double> decay_coeff, vector<double> decay_order)
{
	create_f_inv_u();
	create_k_mat_copula();
	create_inv_K_Copula(dt,Diffusion_coeff);
	C.resize(numberofspecies);
	for (int species_counter=0; species_counter<numberofspecies; species_counter++)
        C[species_counter] = CMatrix(GP.nx+2, GP.ny);
    for (int i = 0; i < GP.nx; i++)
        for (int j = 0; j < GP.ny; j++)
            p[i][j].C.resize(int(t_end/dt), numberofspecies);
	OU.BTCs = CBTCSet(GP.nx+2);
	OU.BTC_normal = CBTCSet(GP.nx + 2);
	OU.BTCs_fw = CBTCSet(GP.nx + 2);
	OU.BTC_normal_fw = CBTCSet(GP.nx + 2);
	for (int i = 0; i < GP.nx+2; i++) OU.BTCs.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	//K.writetofile(pathout + "transport_matrix.txt");
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs.BTC[i].append(0, 0);
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.BTC[i].append(0, 0);
	set_progress_value(0);
	int counter = 0;
    for (double t = 0; t < t_end; t += dt)
    {

        for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
        {
            CVector_arma RHS = create_RHS_Copula(species_counter, dt, Diffusion_coeff, decay_coeff[species_counter], decay_order[species_counter]);
            CVector_arma S = RHS/copula_params.Inv_M;

            for (int i = 0; i < GP.nx+2; i++)
            {
                double sum = 0;
                double sum_fw = 0;

                for (int j = 0; j < GP.ny; j++)
                {
                    C[species_counter][i][j] = S[get_cell_no_OU(i, j)];
                    sum += C[species_counter][i][j] * GP.dy;
                    sum_fw += C[species_counter][i][j] * OU.FinvU[j]/OU.FinvU.sum();
                }
                OU.BTCs.BTC[i].append(t+dt, sum);
                OU.BTCs_fw.BTC[i].append(t + dt, sum_fw);
                OU.BTC_normal.BTC[i].append((t + dt)/((i-0.5)*GP.dx), sum);
                OU.BTC_normal_fw.BTC[i].append((t + dt) / ((i - 0.5)*GP.dx), sum_fw);
            }

            #if QT_version
                    set_progress_value(t / t_end);
            tbrowse->append("t = " + QString::number(t));
            #else
            set_progress_value(t);
            #endif // QT_version

        }
        for (int i = 0; i < GP.nx; i++)
            for (int j = 0; j < GP.ny; j++)
            {
                for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
                    p[i][j].C.setvalue(counter,species_counter,C[species_counter][i+1][j]);
            }
        counter++;
    }

}


void CGrid::solve_transport_Copula_diffusion(double t_end, double Diffusion_coeff, vector<double> decay_coeff, vector<double> decay_order)
{
	create_f_inv_u();
	create_k_mat_copula_only_dispersion();
	create_k_mat_copula_only_diffusion();
	create_inv_K_Copula_diffusion(dt,Diffusion_coeff);
	C.resize(numberofspecies);
	for (int i = 0; i < GP.nx; i++)
        for (int j = 0; j < GP.ny; j++)
            p[i][j].C.resize(int(t_end/dt),numberofspecies);
	for (int species_counter=0; species_counter<numberofspecies; species_counter++)
        C[species_counter] = CMatrix(GP.nx+2, GP.ny);
	OU.BTCs = CBTCSet(GP.nx+2);
	OU.BTC_normal = CBTCSet(GP.nx + 2);
	OU.BTCs_fw = CBTCSet(GP.nx + 2);
	OU.BTC_normal_fw = CBTCSet(GP.nx + 2);
	for (int i = 0; i < GP.nx+2; i++) OU.BTCs.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs.BTC[i].append(0, 0);
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.BTC[i].append(0, 0);
	set_progress_value(0);
    int counter=0;
    for (double t = 0; t < t_end; t += dt)
    {
        for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
        {
            CVector_arma RHS = create_RHS_Copula_diffusion(species_counter, dt, Diffusion_coeff, decay_coeff[species_counter], decay_order[species_counter]);
            CVector_arma S = RHS/copula_params.Inv_M;

            for (int i = 0; i < GP.nx+2; i++)
            {
                double sum = 0;
                double sum_fw = 0;

                for (int j = 0; j < GP.ny; j++)
                {
                    C[species_counter][i][j] = S[get_cell_no_OU(i, j)];
                    sum += C[species_counter][i][j] * GP.dy;
                    sum_fw += C[species_counter][i][j] * OU.FinvU[j]/OU.FinvU.sum();
                }
                OU.BTCs.BTC[i].append(t+dt, sum);
                OU.BTCs_fw.BTC[i].append(t + dt, sum_fw);
                OU.BTC_normal.BTC[i].append((t + dt)/((i-0.5)*GP.dx), sum);
                OU.BTC_normal_fw.BTC[i].append((t + dt) / ((i - 0.5)*GP.dx), sum_fw);


                #if QT_version
                        set_progress_value(t / t_end);
                tbrowse->append("t = " + QString::number(t));
                #else
                set_progress_value(t/t_end);
                #endif // QT_version
            }
            for (int i = 0; i < GP.nx; i++)
                for (int j = 0; j < GP.ny; j++)
                {

                    for (int species_counter = 0; species_counter<numberofspecies; species_counter++)
                        p[i][j].C.setvalue(counter,species_counter,C[species_counter][i+1][j]);

                }

        }
        counter++;
    }
}


void CGrid::renormalize_k()
{
	CBTC K = get_kg_btc(0);
	double mu = K.mean();
	double std = K.std();
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			p[i][j].K_gauss[0] = (p[i][j].K_gauss[0] - mu) / std;
	remap_K(0);
}


void CGrid::show_in_window(string s)
{
    #ifdef QT_version
    qDebug()<<QString::fromStdString(s);
    main_window->get_ui()->ShowOutput->append(QString::fromStdString(s));
    QApplication::processEvents();
    #else
    cout<<s<<endl;
    #endif // Qt_version
}

void CGrid::set_progress_value(double s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s*100 << "%                                     ";
}

void CGrid::set_progress_value(string s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r" << s << "                                     ";
}

void CGrid::clear_contents()
{
	this->n_k_dets = 0;
	p.resize(GP.nx);
	for (int i = 0; i < GP.nx; i++)
	{
		p[i].resize(GP.ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].k_det = false;
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K_gauss = CVector(2);
		}
	}
}

CTimeSeries CGrid::GetConcentrationBTCAtX(int species_counter, double x, const string &filename, const string &filename_d)
{
    CTimeSeries output;
    for (int tt=0; tt<p[0][0].C.size(); tt++)
    {
        output.append(tt*dt,GetConcentrationAtX(species_counter, x,tt));
    }
    output.writefile(pathout + filename);
    if (filename_d!="")
        output.derivative().writefile(pathout + filename_d);
    return output;
}

double CGrid::GetConcentrationAtX(int species_counter, double x, int timestep)
{
    int i=x/GP.dx;
    double output = 0;
    for (int j=0; j<GP.ny; j++)
        output += p[i][j].C[timestep][species_counter]/GP.ny;

    return output;
}

double CGrid::mean(double u1, double u2)
{
	if (copula_params.mean_method == "harmonic")
		return 2 * u1*u2 / (u1 + u2);
	if (copula_params.mean_method == "geometric")
		return sqrt(u1*u2);
	if (copula_params.mean_method == "arithmetic")
		return 0.5*(u1 + u2);
	if (copula_params.mean_method == "1")
		return u1;
	if (copula_params.mean_method == "2")
		return u2;
    return 0.5*(u1+u2);
}

CTimeSeries CGrid::GetProfile(int species_id, int timestep, double x_start, double x_end, double interval, const string &filename)
{
    CBTC Profile;
    for (double x=x_start; x<=x_end; x+=interval)
    {
        Profile.append(x,GetConcentrationAtX(species_id, x,timestep));
    }
    Profile.writefile(filename);
    return Profile;
}

CTimeSeries CGrid::GetAllVelocities(const string &dir)
{
    CTimeSeries out;
    if (dir=="x")
    {
        for (int i=0; i<GP.nx; i++)
            for (int j = 0; j < GP.ny - 1; j++)
                out.append(i*(GP.ny-1)+j,vx[i][j]);
    }
    else if (dir=="y")
    {
        for (int i = 0; i<GP.nx-1; i++)
            for (int j = 0; j < GP.ny; j++)
                out.append(i*(GP.ny)+j,vx[i][j]);
    }
    return out;

}

CTimeSeries CGrid::GetResiduals(const string &property)
{
    CTimeSeries residuals;
    for (int i=1; i<GP.nx-1; i++)
        for (int j=1; j<GP.ny-1; j++)
        {
            double avg = 0.25*(getproperty(property,i-1,j)+getproperty(property,i+1,j)+getproperty(property,i,j+1)+getproperty(property,i,j-1));
            double middle = getproperty(property,i,j);
            residuals.append(i+10000*j,avg-middle);
        }
    return residuals;
}

double CGrid::getproperty(const string &property,int i, int j)
{
    if (i>GP.nx-1)
    {
        cout<<"index out of range:"<<i<<","<<j<<endl;
        return 0;
    }
    if (j>GP.ny-1)
    {
        cout<<"index out of range:"<<i<<","<<j<<endl;
        return 0;
    }
    if (property == "k_normal")
        return p[i][j].K_gauss[0];
    if (property == "k")
        return p[i][j].K[0];
    if (property == "v_x")
        return p[i][j].V[0];
    if (property == "v_y")
        return p[i][j].V[1];
    if (property == "u")
        return p[i][j].u;
    if (property == "omega")
        return p[i][j].omega;
    return 0;

}

void CGrid::Assign_Ranks()
{
    CBTC All_V = get_v_btc(0);
    CBTC V_cummulative = All_V.getcummulative_direct(100,false);
    for (int i=0; i<GP.nx; i++)
        for (int j=0; j<GP.ny; j++)
        {
            double rank = V_cummulative.interpol(p[i][j].V[0]);
            p[i][j].u = rank;
            p[i][j].omega = stdnormal_inv(rank);
        }

}

