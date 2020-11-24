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


Concentrations::Concentrations(const double& other)
{
    values[0] = other;
}

Concentrations& Concentrations::operator=(const Concentrations& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    values = rhs.values;
    return *this;
}

Concentrations& Concentrations::operator=(const double& rhs)
{
    values[0] = rhs;
    return *this;
}

