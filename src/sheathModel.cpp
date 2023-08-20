#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "sheathModel.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef gitr_precision gitr_precision;
#endif

gitr_precision _boltz = 1.3806504e-23;
gitr_precision e_charge = 1.6021765e-19;

gitr_precision CouletteManfredi(gitr_precision distance, gitr_precision alpha ) {
    gitr_precision C1_alpha = -0.00281407f * alpha - 2.31655435f;
    gitr_precision C2_alpha = 0.00640402f * alpha + 0.01023915f;
    gitr_precision result = C1_alpha * C2_alpha * std::exp(-C2_alpha * distance);


    return result;
}

    gitr_precision BrooksModel(gitr_precision distance, gitr_precision alpha, gitr_precision larmorRadius, gitr_precision Te) {
    gitr_precision Vmps = log( abs(alpha )* 180.0 / M_PI) * _boltz * 11600.0* Te / e_charge;
    gitr_precision result = Vmps * std::exp( - distance / larmorRadius) / larmorRadius;

    return result;


}
