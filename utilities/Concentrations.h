#ifndef CONCENTRATIONS_H
#define CONCENTRATIONS_H

#include <vector>

using namespace std;

class Concentrations
{
    public:
        Concentrations();
        Concentrations(vector<vector<double>> vals) {values = vals; }
        virtual ~Concentrations();
        Concentrations(const Concentrations& other);
        Concentrations& operator=(const Concentrations& other);
        void resize(int n) {values.resize(n);}
        void resize(int nt, int nc);
        void append(const vector<double> &value);
        vector<double> &operator [](int i) {return values[i];}
        vector<double> value(int i, int j) {if (i<values.size()) return values[i]; else return vector<double>();  }
        bool setvalue(int i, int j, double value) {if (i<values.size()) values[i][j]=value; else return false; return true;  }
        int size() {return values.size();}
    protected:

    private:
        vector<vector<double>> values;
};

#endif // CONCENTRATIONS_H
