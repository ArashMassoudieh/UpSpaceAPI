#ifndef CONCENTRATIONS_H
#define CONCENTRATIONS_H

#include <vector>

using namespace std;

class Concentrations
{
    public:
        Concentrations();
<<<<<<< HEAD
        Concentrations(vector<double> vals) {values = vals; }
=======
>>>>>>> master
        virtual ~Concentrations();
        Concentrations(const Concentrations& other);
        Concentrations(const double& other);
        Concentrations& operator=(const Concentrations& other);
        Concentrations& operator=(const double& other);
        void resize(int n) {values.resize(n);}
        double &operator [](int i) {return values[i];}
        double value(int i) {if (i<values.size()) return values[i]; else return -999;  }
        bool setvalue(int i, double value) {if (i<values.size()) values[i]=value; else return false; return true;  }
    protected:

    private:
        vector<double> values;
};

#endif // CONCENTRATIONS_H
