#include "Concentrations.h"

Concentrations::Concentrations()
{
    values.resize(1);
}

Concentrations::~Concentrations()
{
    //dtor
}

Concentrations::Concentrations(const Concentrations& other)
{
    values = other.values;
}


Concentrations& Concentrations::operator=(const Concentrations& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    values = rhs.values;
    return *this;
}

void Concentrations::resize(int nt, int nc)
{
    values.resize(nt);
    for (int i=0; i<nt; i++)
        values[i].resize(nc);
}

