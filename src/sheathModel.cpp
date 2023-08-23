#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "sheathModel.h"
#include "constants.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef gitr_precision gitr_precision;
#endif


gitr_precision CouletteManfredi(gitr_precision distance, gitr_precision alpha ) {
    gitr_precision C1_alpha = -0.00281407f * alpha - 2.31655435f;
    gitr_precision C2_alpha = 0.00640402f * alpha + 0.01023915f;
    gitr_precision result = C1_alpha * C2_alpha * std::exp(-C2_alpha * distance);


    return result;
}

    gitr_precision BrooksModel(gitr_precision distance, gitr_precision alpha, gitr_precision larmorRadius, gitr_precision Te) {
     gitr_precision Vmps = 0.0;
    gitr_precision eps = 1.e-3; 
    if (alpha < eps) {
        Vmps = 0.0;
    }
    else {
        Vmps = log( abs(alpha )* 180.0 / M_PI) * gitr_constants::boltz * 11600.0* Te / gitr_constants::e_charge;
    }
    gitr_precision result = Vmps * std::exp( - distance / larmorRadius) / larmorRadius;

    return result;


}
