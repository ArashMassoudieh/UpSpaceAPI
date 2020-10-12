#ifndef DMAP_H
#define DMAP_H
#include <vector>
#include <string>


class CCopula;
class MapAsTimeSeriesSet;

using namespace std;

class TDMap
{
    public:
        TDMap();
        TDMap(unsigned int number_of_bins, double low_lim, double up_lim);
        TDMap(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double low_lim_x, double up_lim_x, double low_lim_y, double up_lim_y);
        virtual ~TDMap();
        TDMap(const TDMap& other);
        TDMap& operator=(const TDMap& other);
        void reset(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double low_lim_x, double up_lim_x, double low_lim_y, double up_lim_y);
        void set_val(unsigned int i, unsigned int j, double val);
        void add_val(unsigned int i, unsigned int j, double val);
        void add_val(double x, double y, double val);
        void normalize();
        void normalise_x();
        void normalize_y();
        double get_val(unsigned int i, unsigned int j);
        double sum();
        double marginal_x(unsigned int i);
        double marginal_y(unsigned int j);
        vector<double> marginal_x();
        vector<double> marginal_y();
        void writetofile(string filename);
        void writetofile_GNU(string filename,string pngfilename="", string xlabel="", string ylabel="", string title="",bool logscale=false);
        void writetheoreticalcopulatofile(string filename, CCopula *copula);
        MapAsTimeSeriesSet getcumulative(string dir);
        double interpolate(double x, double y);
        bool readfromfile(const string &filename);
        int nx() {return int(val.size());}
        int ny() {
            if (val.size()==0)
                return 0;
            else
                return val[0].size();
        }
    protected:

    private:
        vector<vector<double>> val;
        vector<double> x_bin;
        vector<double> y_bin;
        double low_lim_x;
        double up_lim_x;
        double low_lim_y;
        double up_lim_y;


};

TDMap operator + (TDMap &m1, TDMap &m2);
TDMap operator * (TDMap &m1, double d);
TDMap operator / (TDMap &m1, double d);
#endif // 2DMAP_H
