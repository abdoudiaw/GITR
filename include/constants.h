
#ifndef CONSTANTS_H
#define CONSTANTS_H


#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

namespace gitr_constants
{
const gitr_precision pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
const gitr_precision boltz = 1.3806504e-23;
const gitr_precision hplanck = 6.62606896e-34;
const gitr_precision e_mass = 9.1093837015e-31;
const gitr_precision p_mass = 1.67262192369e-27;
const gitr_precision e_charge = 1.6021765e-19;
const gitr_precision VACUUM_PERMITTIVITY = 8.854187e-12;
}

#endif

