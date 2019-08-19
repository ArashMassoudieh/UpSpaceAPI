#pragma once
#include "Pathway.h"
#include <vector>
#include <string>
#include "Copula.h"

using namespace std;
class CPathwaySet
{
public:
    bool weighted;
    CPathwaySet();
    CPathwaySet(const CPathwaySet &P);
    CPathwaySet &operator=(const CPathwaySet &P);
    ~CPathwaySet();
    vector<CPathway> paths;
    void write(string filename);
    void append(const CPathway& P, double weight=1);
    int max_num_points();
    void create_ou_paths(int n, CDistribution *dist, double x_min, double x_max, double kappa, double dx, double weight=1);
    void create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, double r, double dx, double weight=1);
    void create_copula_paths(int n, CDistribution * dist, double x_min, double x_max, double epsilon, CCopula *copula, double dx, double weight=1);
    void write_vtk(vtkSmartPointer<vtkPolyDataMapper>, string filename);
    vtkSmartPointer<vtkPolyDataMapper> pathways_vtk_pdt_vtp(double z_factor=1, double offset=0);
    CPathway snapshotattime(double t);
    CPathway snapshotatlocation(double x);
    void make_uniform_at_x(double dx);
    void make_uniform_at_t(double dt);
    CPosition get_pair_v_pos(int increment, int num_seq=2);
    CBTCSet get_pair_v(int increment, int n, int num_seq=2);
    CBTC get_BTC(double x, int n_bins, bool vel_inv_weighted = true, double smoothing_factor=0);
    CBTC get_BTC_points(double x);


};

