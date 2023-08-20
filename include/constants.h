#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

namespace gitr_constants
{
    gitr_precision boltz = 1.3806504e-23;
    gitr_precision hplanck = 6.62606896e-34;
    gitr_precision e_mass = 9.1093837015e-31;
    gitr_precision p_mass = 1.67262192369e-27; 
    gitr_precision e_charge = 1.6021765e-19;
    gitr_precision VACUUM_PERMITTIVITY = 8.854187e-12;

}
