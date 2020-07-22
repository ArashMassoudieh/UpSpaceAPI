#include "PathwaySet.h"
#include "StringOP.h"
#include "Grid.h"



CPathwaySet::CPathwaySet()
{
}

CPathwaySet::CPathwaySet(const CPathwaySet & P)
{
	paths = P.paths;
    weighted = P.weighted;

}

CPathwaySet &CPathwaySet::operator=(const CPathwaySet & P)
{
	paths = P.paths;
	return *this;
}


CPathwaySet::~CPathwaySet()
{
}

void CPathwaySet::write(string filename)
{
	ofstream file;
	file.open(filename.c_str());
	for (int j = 0; j < paths.size(); j++)
	{
		file << "t_" << numbertostring(j) << ",x_" << numbertostring(j) << ",y_" << numbertostring(j) << ",u_" << numbertostring(j) << ",v_" << numbertostring(j) << ",u_" << numbertostring(j) << ",z_" << numbertostring(j);
	}
	file << endl;
	for (int i = 0; i < max_num_points(); i++)
	{
		for (int j = 0; j < paths.size(); j++)
		{
			if (i < paths[j].positions.size())
				file << paths[j].positions[i].t << "," << paths[j].positions[i].x << "," << paths[j].positions[i].y << "," << paths[j].positions[i].v[0] << "," << paths[j].positions[i].v[1] << "," << paths[j].positions[i].u << "," << paths[j].positions[i].z << ",";
			else
				file << "," << "," << "," << "," << "," << "," << ",";
		}
		file << endl;
	}

file.close();
}

void CPathwaySet::append(const CPathway & P, double weight)
{
   paths.push_back(P);
}

int CPathwaySet::max_num_points()
{
	int max_np = 0;
	for (int i = 0; i < paths.size(); i++)
		max_np = max(max_np, int(paths[i].positions.size()));

	return max_np;
}

void CPathwaySet::create_ou_paths(int n, CDistribution * dist, double x_min, double x_max, double kappa, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_ou(dist, x_min, x_max, kappa, dx);
            P.weight = weight;
            append(P);
	}
}

void CPathwaySet::create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, double r, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_copula(dist, x_min, x_max, epsilon, r, dx);
            P.weight = weight;
            append(P);
	}
}

void CPathwaySet::create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double weight)
{
	for (int i = 0; i < n; i++)
	{
            CPathway P;
            P.create_copula(dist, x_min, x_max, epsilon, copula, dx);
            P.weight = weight;
            append(P);
	}
}


void CPathwaySet::write_vtk(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename)
{


	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

}

vtkSmartPointer<vtkPolyDataMapper> CPathwaySet::pathways_vtk_pdt_vtp(double z_factor, double offset)
{
	vector<vtkSmartPointer<vtkPolyData>> outarray;
	for (int i = 0; i < paths.size(); i++)
		outarray.push_back(paths[i].pathway_vtk_pdt_vtp(z_factor, offset));

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i = 0; i < outarray.size(); i++)
		appendFilter->AddInputData(outarray[i]);
#endif
	appendFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
#endif

	return mapper;
}

CPathway CPathwaySet::snapshotattime(double t)
{
    CPathway Ptwy;
    for (int i = 0; i < paths.size(); i++)
	Ptwy.append(paths[i].get_position_at_t(t));

    return Ptwy;
}

CPathway CPathwaySet::snapshotatlocation(double x)
{
	CPathway Ptwy;
	for (int i = 0; i < paths.size(); i++)
            Ptwy.append(paths[i].get_position_at_x(x));

	return Ptwy;
}

void CPathwaySet::make_uniform_at_x(double dx)
{
	for (int i = 0; i < paths.size(); i++)
		paths[i] = paths[i].make_uniform_x(dx);
}

void CPathwaySet::make_uniform_at_t(double dt)
{
	for (int i = 0; i < paths.size(); i++)
		paths[i] = paths[i].make_uniform_t(dt);

}



CPosition CPathwaySet::get_pair_v_pos(int increment, int num_seq)
{
    CPosition p(num_seq);
    int i = int(unitrandom()*paths.size());
    int j = int((paths[i].positions.size()- (num_seq-1)*increment)*unitrandom());
    for (int ii=0 ; ii<num_seq; ii++)
       p.v[ii] = paths[i].positions[j+ii*increment].v[0];

    p.weight = paths[i].weight;
    return p;
}

CBTCSet CPathwaySet::get_pair_v(int increment, int n, int num_seq)
{
	CBTCSet out(num_seq);
	for (int i = 0; i < n; i++)
	{
		CPosition p = get_pair_v_pos(increment, num_seq);
		if (weighted)
            out.append(double(i), p.v.vec, p.weight);
        else
            out.append(double(i), p.v.vec);
		cout << "\r" << float(i)/float(n)*100 << "%" << std::flush;
	}
	return out;

}

CBTC CPathwaySet::get_BTC(double x, int n_bins, bool velweight, double smoothing_factor)
{
    CBTC BTC;
    if (weighted)
    {
        BTC.weighted = true;
        for (int i = 0; i < paths.size(); i++)
        {
            if (velweight)
                BTC.append(i, paths[i].get_cross_time_vx(x)[0],1.0/paths[i].get_cross_time_vx(x)[1]);
            else
                BTC.append(i, paths[i].get_cross_time(x),paths[i].weight);
        }
    }
    else
    {
        BTC.weighted = false;
        for (int i = 0; i < paths.size(); i++)
        {
            if (velweight)
                BTC.append(i, paths[i].get_cross_time_vx(x)[0],1.0/paths[i].get_cross_time_vx(x)[1]);
            else
                BTC.append(i, paths[i].get_cross_time(x));
        }
    }

    return BTC.distribution(n_bins,(BTC.maxC()-BTC.minC())*smoothing_factor, 0);

}

CBTC CPathwaySet::get_BTC_points(double x, bool vel_inv_weighted)
{
    CBTC BTC;
    if (weighted)
    {   BTC.weighted = true;
        for (int i = 0; i < paths.size(); i++)
        {
            if (vel_inv_weighted)
                BTC.append(i, paths[i].get_cross_time_vx(x)[0],1.0/paths[i].get_cross_time_vx(x)[1]);
            else
                BTC.append(i, paths[i].get_cross_time(x),paths[i].weight);
        }

    }
    else
    {   BTC.weighted = false;
        for (int i = 0; i < paths.size(); i++)
        {
            if (vel_inv_weighted)
                BTC.append(i, paths[i].get_cross_time_vx(x)[0],1.0/paths[i].get_cross_time_vx(x)[1]);
            else
                BTC.append(i, paths[i].get_cross_time(x));
        }
    }


    return BTC;

}

bool CPathwaySet::getfromMODflowfile(const string &filename)
{
    ifstream file;
    file.open (filename, std::fstream::in);
    if (!file.good()) return false;
    int rownum = 0;
    while (!file.eof())
    {
        vector<double> s1 = ATOF(getline(file,' '));
        if (n() == 0)
            if (s1[0]!=0)
                paths.resize(s1[0]);
        if (s1.size()==2)
        {
            int numberofpoints = s1[0];
            double age = s1[1];

            for (int i=0; i<numberofpoints; i++)
            {
                vector<double> s2 = ATOF(getline(file,' '));
                CPosition P;
                P.x = s2[1];
                P.y = s2[2];
                P.z = s2[3];
                P.v = CVector(2);
                P.t = age;
                if (paths[i].size()>0)
                {
                    P.v[0] = (P.x - paths[i].positions[paths[i].size()-1].x)/(P.t - paths[i].positions[paths[i].size()-1].t);
                    P.v[1] = (P.y - paths[i].positions[paths[i].size()-1].y)/(P.t - paths[i].positions[paths[i].size()-1].t);

                }

                P.weight = s2[5];
                paths[(int)s2[0]-1].append(P);
            }

        }
        rownum++;
    }
    file.close();
}



