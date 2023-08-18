#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "sheathModel.h"

float CouletteManfredi(float distance, float alpha) {
    float C1_alpha = -0.00281407f * alpha - 2.31655435f;
    float C2_alpha = 0.00640402f * alpha + 0.01023915f;
    float result = C1_alpha * std::exp(-C2_alpha * distance);

    if (result > 0) {
        result = 0;
    }
    return result;
}

