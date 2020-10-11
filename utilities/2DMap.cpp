#include "2DMap.h"
#include <fstream>
#include <string>
#include "StringOP.h"
#include "MapAsTimeSeriesSet.h"
#include "BTC.h"
#include "Copula.h"

TDMap::TDMap()
{
    //ctor
}

TDMap::~TDMap()
{
    //dtor
}

TDMap::TDMap(unsigned int number_of_bins, double low_lim, double up_lim)
{
    reset(number_of_bins, number_of_bins, low_lim, up_lim, low_lim, up_lim);
}

TDMap::TDMap(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double low_lim_x, double up_lim_x, double low_lim_y, double up_lim_y)
{
    reset(number_of_bins_x, number_of_bins_y, low_lim_x, up_lim_x, low_lim_y, up_lim_y);
}

TDMap::TDMap(const TDMap& other)
{
    val = other.val;
    x_bin = other.x_bin;
    y_bin = other.y_bin;
    up_lim_x = other.up_lim_x;
    low_lim_x = other.low_lim_x;
    up_lim_y = other.up_lim_y;
    low_lim_y = other.low_lim_y;
}

TDMap& TDMap::operator=(const TDMap& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    val = rhs.val;
    x_bin = rhs.x_bin;
    y_bin = rhs.y_bin;
    up_lim_x = rhs.up_lim_x;
    low_lim_x = rhs.low_lim_x;
    up_lim_y = rhs.up_lim_y;
    low_lim_y = rhs.low_lim_y;
    return *this;
}

void TDMap::reset(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double _low_lim_x, double _up_lim_x, double _low_lim_y, double _up_lim_y)
{
    low_lim_x = _low_lim_x;
    up_lim_x = _up_lim_x;
    low_lim_y = _low_lim_y;
    up_lim_y = _up_lim_y;
    val.resize(number_of_bins_x);
    for (unsigned int i=0; i<number_of_bins_x; i++)
        val[i].resize(number_of_bins_y);

    x_bin.resize(number_of_bins_x+1);
    y_bin.resize(number_of_bins_y+1);

    for (unsigned int i=0; i<number_of_bins_x+1; i++)
        x_bin[i] = low_lim_x + (up_lim_x-low_lim_x)/number_of_bins_x*i;


    for (unsigned int i=0; i<number_of_bins_y+1; i++)
        y_bin[i] = low_lim_y + (up_lim_y-low_lim_y)/number_of_bins_y*i;
}

void TDMap::set_val(unsigned int i, unsigned int j, double value)
{
    val[i][j] = value;
}

void TDMap::add_val(unsigned int i, unsigned int j, double value)
{
    val[i][j] += value;
}

void TDMap::add_val(double x, double y, double value)
{
    int i = int((x-low_lim_x)/(up_lim_x-low_lim_x))*val.size();
    int j = int((y-low_lim_y)/(up_lim_y-low_lim_y))*val.size();
    if (i>=0 && j>=0 && i<val.size() && j<val[0].size())
        val[i][j] += value;
}

double TDMap::sum()
{
    double sum = 0;
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; j<val[i].size(); j++)
            sum+=val[i][j];

    return sum;
}

double TDMap::marginal_x(unsigned int i)
{
    double sum = 0;
    for (unsigned int j=0; j<val[i].size(); j++)
        sum+=val[i][j];

    return sum;
}

double TDMap::marginal_y(unsigned int j)
{
    if (j>val[0].size())
        return 0;
    double sum = 0;
    for (unsigned int i=0; i<val.size(); i++)
        sum+=val[i][j];

    return sum;
}

vector<double> TDMap::marginal_x()
{
    vector<double> out;
    for (unsigned int i=0; i<val.size(); i++)
        out.push_back(marginal_x(i));

    return out;
}

vector<double> TDMap::marginal_y()
{
    vector<double> out;
    for (unsigned int j=0; j<val[0].size(); j++)
        out.push_back(marginal_y(j));

    return out;
}

void TDMap::normalize()
{
    double s = sum();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; j<val[i].size(); j++)
            val[i][j] = val[i][j]/s*val.size()*val[0].size()/(up_lim_x-low_lim_x)/(up_lim_y-low_lim_y);

}

void TDMap::normalise_x()
{
    vector<double> s = marginal_x();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; val[i].size(); j++)
            val[i][j] = val[i][j]/s[i];

}

void TDMap::normalize_y()
{
    vector<double> s = marginal_y();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; val[i].size(); j++)
            val[i][j] = val[i][j]/s[j];
}

MapAsTimeSeriesSet TDMap::getcumulative(string dir)
{
    MapAsTimeSeriesSet out;
    double dy = (up_lim_y-low_lim_y)/val.size();
    if (dir=="x")
    {
        for (int j=0; j<val.size(); j++)
        {
            double dx = (up_lim_x-low_lim_x)/val[j].size();
            CTimeSeries X;
            X.append(0,0);
            for (int i=0; i<val[j].size(); i++)
            {
                X.append(dx*i,val[j][i]*dx+X.C[X.n-1]);
            }
            X = X/X.C[X.n-1];
            out.append(X,dy*j+dy/2);
        }
    }
    if (dir=="sym")
    {
        for (int j=0; j<val.size(); j++)
        {
            double dx = (up_lim_x-low_lim_x)/val[j].size();
            CTimeSeries X;
            X.append(0,0);
            for (int i=0; i<val[j].size(); i++)
            {
                X.append(dx*i+dx,0.5*(val[j][i]+val[i][j])*dx+X.C[X.n-1]);
            }
            X = X/X.C[X.n-1];
            out.append(X,dy*j+dy/2);
        }
    }
    return out;
}

double TDMap::get_val(unsigned int i, unsigned int j)
{
    if (i<val.size() && j<val[0].size())
        return val[i][j];
    return 0;
}

void TDMap::writetofile(string filename)
{
    ofstream file(filename);

    for (unsigned int i=0; i<val.size(); i++)
    {
        file << "," << (x_bin[i] + x_bin[i+1])/2;
    }
    file << endl;
    for (unsigned int j=0; j<val[0].size(); j++)
    {
        file << (y_bin[j] + y_bin[j+1])/2 << ",";
        for (unsigned int i=0; i<val.size(); i++)
            file << val[i][j] << ",";
        file << endl;
    }
    file.close();
}

bool TDMap::readfromfile(const string &filename)
{
    ifstream file(filename);
    if (!file.good())
        {
            cout << "The program was not able to open " << filename << endl;
            return false;
        }
    vector<double> x = ATOF(getline(file));
    reset(x.size()-1,x.size()-1,0,1,0,1);
    for (int i=0; i<x.size()-1; i++)
    {
        vector<double> v = ATOF(getline(file));
        for (int j=0; j<x.size()-1; j++)
        {
            val[i][j] = v[j+1];
        }
    }
    return true;
}

void TDMap::writetheoreticalcopulatofile(string filename, CCopula *copula)
{
    ofstream file(filename);
    double x,y;
    for (unsigned int i=0; i<val.size(); i++)
    {

        file << "," << (x_bin[i] + x_bin[i+1])/2;
    }
    file << endl;
    for (unsigned int j=0; j<val[0].size(); j++)
    {
        file << (y_bin[j] + y_bin[j+1])/2 << ",";
        y = (y_bin[j] + y_bin[j+1])/2.0;
        for (unsigned int i=0; i<val.size(); i++)
        {
            x = (x_bin[i] + x_bin[i+1])/2.0;
            file << copula->evaluate11(x,y) << ",";
        }
        file << endl;
    }
    file.close();
}

void TDMap::writetofile_GNU(string filename, string pngfilename, string xlabel, string ylabel, string title, bool logscale)
{
    if (pngfilename=="")
        pngfilename = split(filename,'.')[0] + ".png";
    ofstream file(filename);
    if (logscale)
    {
        file << "set xrange ["<<low_lim_x+0.49*(up_lim_x-low_lim_x)/val.size()<<":"<<up_lim_x-0.49*(up_lim_x-low_lim_x)/val.size()<<"]"<<endl;
        file << "set yrange ["<<low_lim_y+0.49*(up_lim_y-low_lim_y)/val[0].size()<<":"<<up_lim_y-0.49*(up_lim_y-low_lim_y)/val[0].size()<<"]"<<endl;
    }
    else
    {
        file << "set xrange ["<<low_lim_x+0*0.5*(up_lim_x-low_lim_x)/val.size()<<":"<<up_lim_x-0*0.5*(up_lim_x-low_lim_x)/val.size()<<"]"<<endl;
        file << "set yrange ["<<low_lim_y+0*0.5*(up_lim_y-low_lim_y)/val[0].size()<<":"<<up_lim_y-0*0.5*(up_lim_y-low_lim_y)/val[0].size()<<"]"<<endl;
    }

    file << "show xrange" << endl;
    file << "show yrange" << endl;
    file << "set xlabel \"" << xlabel << "\"" << endl;
    file << "show xlabel" << endl;
    file << "set ylabel \"" << ylabel << "\"" << endl;
    file << "show ylabel" << endl;
    file << "set title \"" << title << "\"" << endl;
    file << "show title" << endl;
    file << "set palette rgb -21,-22,-23" << endl;
    file << "set logscale cb" << endl;
    if (logscale) file << "set logscale xy" << endl;
    file << "set format cb \"10^{%T}\"" << endl;
    file << "set pm3d map interpolate 5,5"<<endl;
    file << "set size square" << endl;
    file << "splot \"-\" matrix using ($1*"<<(up_lim_x-low_lim_x)/val.size()<<"+"<<low_lim_x+0.5*(up_lim_x-low_lim_x)/val.size()<<"):($2*"<<(up_lim_y-low_lim_y)/val.size()<<"+"<<low_lim_x+0.5*(up_lim_y-low_lim_y)/val.size()<<"):3 notitle"<<endl;
    for (unsigned int j=0; j<val[0].size(); j++)
    {
        for (unsigned int i=0; i<val.size(); i++)
            file << val[i][j] << " ";
        file << endl;
    }
    file<<"e"<<endl;
    file<<"e"<<endl;
    file<<"set term png"<<endl;
    file<<"set output \"" << pngfilename << "\""<<endl;
    file<<"replot"<<endl;
    file<<"set term x11"<<endl;
    file.close();
}

double TDMap::interpolate(double x, double y)
{
    double dx = (up_lim_x-low_lim_x)/val[0].size();
    double dy = (up_lim_y-low_lim_y)/val.size();
    int i_1 = max(int( (x - low_lim_x - dx/2)/dx),0);
    int i_2 = min(i_1+1,int(val[0].size()-1));
    int j_1 = max(int( (y - low_lim_y - dy/2)/dy),0);
    int j_2 = min(j_1+1,int(val.size()-1));
    double interpol1 = 0.5*(val[j_1][i_1]+val[i_1][j_1]) + 0.5*(val[j_1][i_2]-val[j_1][i_1]+val[i_2][j_1]-val[i_1][j_1])/dx*(x-i_1*dx-dx/2);
    double interpol2 = 0.5*(val[j_2][i_1]+val[i_1][j_2]) + 0.5*(val[j_2][i_2]-val[j_2][i_1]+val[i_2][j_2]-val[i_1][j_2])/dx*(x-i_1*dx-dx/2);
    double interpol = interpol1 + (interpol2-interpol1)/dy*(y-j_1*dy-dy/2);
    return interpol;
}
