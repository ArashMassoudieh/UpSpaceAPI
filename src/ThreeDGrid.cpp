#include "ThreeDGrid.h"
#include "StringOP.h"
#include "fstream"
#include "iostream"
#include "math.h"


ThreeDGrid::ThreeDGrid()
{
    //ctor
}

ThreeDGrid::~ThreeDGrid()
{
    //dtor
}

ThreeDGrid::ThreeDGrid(const ThreeDGrid& other)
{
    nx = other.nx;
    ny = other.ny;
    nz = other.nz;
    values = other.values;
}

ThreeDGrid& ThreeDGrid::operator=(const ThreeDGrid& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    nx = rhs.nx;
    ny = rhs.ny;
    nz = rhs.nz;
    values = rhs.values;
    return *this;
}

bool ThreeDGrid::ReadFromModFlowFile(const string &filename)
{
    ifstream file(filename);
    if (file.good())
    {
        vector<int> dimensions = ATOI(getline(file,' '));
        if (dimensions.size()!=3)
        {
            cout<<"The file format is not correct! "<< endl;
            return false;
        }

        setdimension(dimensions[0],dimensions[1],dimensions[2]);
        vector<double> gridsize = ATOF(getline(file, ' '));
        setincrements(gridsize[0],gridsize[1],gridsize[2]);
        double t = ATOF(getline(file, ' '))[0];
        vector<string> dummy = getline(file, ' ');
        for (int i=0; i<nx; i++)
        {
            set_progress_value(double(i)/double(nx));
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nz; k++)
                {
                    vector<double> vals = ATOF(getline(file, ' '));
                    values[k][j][i] = vals[3];
                    //cout<< i << "," << j << "," << k <<endl;
                }
            }
        }
    }
    else
        {
            cout<<"File '" << filename << "' cannot be opended!";
            return false;
        }
    return true;
}

bool ThreeDGrid::ReadFromModFlowFile(const string &filename, double _time)
{
    ifstream file(filename);
    if (file.good())
    {
        bool read=false;
        vector<int> dimensions = ATOI(getline(file,' '));
        if (dimensions.size()!=3)
        {
            cout<<"The file format is not correct! "<< endl;
            return false;
        }

        setdimension(dimensions[0],dimensions[1],dimensions[2]);
        vector<double> gridsize = ATOF(getline(file, ' '));
        setincrements(gridsize[0],gridsize[1],gridsize[2]);
        while (!read)
        {
            double t = ATOF(getline(file, ' '))[0];

            vector<string> dummy = getline(file, ' ');
            for (int i=0; i<nx; i++)
            {
                set_progress_value(double(i)/double(nx));
                for (int j=0; j<ny; j++)
                {
                    for (int k=0; k<nz; k++)
                    {
                        vector<double> vals = ATOF(getline(file, ' '));
                    }
                }
            }
        }
    }
    else
        {
            cout<<"File '" << filename << "' cannot be opended!";
            return false;
        }
    return true;
}


void ThreeDGrid::setdimension(int _nx, int _ny, int _nz)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
    values.resize(nz);
    for (int k=0; k<nz; k++)
    {
        values[k].resize(ny);
        {
            for (int j=0; j<ny; j++)
                values[k][j].resize(nx);
        }
    }
}

bool ThreeDGrid::create_vts_test(const string &filename)
{
    ofstream file(filename);
    file << "# vtk DataFile Version 3.0" << endl;
    file << "vtk output" << endl;
    file << "ASCII"<< endl;
    file << "DATASET STRUCTURED_GRID"<<endl;
    file << "DIMENSIONS 10 10 10" << endl;
    file << "POINTS 1000 double" << endl;
    for (int i=0; i<10; i++)
    {
        for (int j=0; j<10; j++)
        {
            for (int k=0; k<10; k++)
            {
                file << i*4 << " " << j*2 << " " << k << " ";
            }
        }
    }
    file << "POINT_DATA " << 1000 << endl;
    file << "FIELD concentration 1" << endl;
    file << "Concentration " << 1 << " " << 1000 << " double" << endl;
     for (int i=0; i<10; i++)
    {
        for (int j=0; j<10; j++)
        {
            for (int k=0; k<10; k++)
            {
                file << sin(double(i+2*j+4*k)/double(100)) << " ";
            }
        }
    }


    file.close();
}


bool ThreeDGrid::create_vts(const string &filename, int x_limit, int y_limit, int z_limit, int x_interval, int y_interval, int z_interval)
{
    if (x_limit == 0 || x_limit>nx) x_limit = nx;
    if (y_limit == 0 || y_limit>ny) y_limit = ny;
    if (z_limit == 0 || z_limit>nz) z_limit = nz;
    if (x_interval == 0) x_interval = 1;
    if (y_interval == 0) y_interval = 1;
    if (z_interval == 0) z_interval = 1;

    cout << x_limit << ", " << y_limit << ", " << z_limit << endl;
    ofstream file(filename);
    file << "# vtk DataFile Version 3.0" << endl;
    file << "vtk output" << endl;
    file << "ASCII"<< endl;
    file << "DATASET STRUCTURED_GRID"<<endl;
    file << "DIMENSIONS " << x_limit/x_interval << " " << y_limit/y_interval << " " << z_limit/z_interval << endl;
    file << "POINTS " << x_limit*y_limit*z_limit/(x_interval*y_interval*z_interval) << " double" << endl;
    for (int i=0; i<x_limit; i+=x_interval)
    {
        for (int j=0; j<y_limit; j+=y_interval)
        {
            for (int k=0; k<z_limit; k+=z_interval)
            {
                file << i*dx + dx/2 << " " << j*dy + dy/2 << " " << k*dz+dz/2 << " ";
            }
        file << endl;
        }
    }
    file << "POINT_DATA " << x_limit*y_limit*z_limit/(x_interval*y_interval*z_interval) << endl;
    file << "FIELD concentration 1" << endl;
    file << "Concentration " << 1 << " " << x_limit*y_limit*z_limit/(x_interval*y_interval*z_interval) << " double" << endl;
    for (int i=0; i<x_limit; i+=x_interval)
    {
        set_progress_value(double(i)/double(nx));
        for (int j=0; j<y_limit; j+=y_interval)
        {
            for (int k=0; k<z_limit; k+=z_interval)
            {
                file << values[k][j][i] << " ";
            }
            file << endl;
        }
    }


    file.close();
    cout<<endl;
    cout<<"Writing the VTK file, finished!" << endl;
    return true;
}
