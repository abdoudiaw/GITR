#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

float ThompsonAngularDistribution(float theta) {
    float b = 1.0;
    float f_theta = (1 + b) * std::pow(std::cos(theta * M_PI / 180), b) * std::sin(theta * M_PI / 180);
    return f_theta;
}

float ThompsonEnergyDistribution(float E, float Emax, float Eb) {
    float Ec = Emax / Eb;
    float C = 2 * Eb / (1 - (1 + 2 * Ec) / std::pow(1 + Ec, 2));
    float f_E = C * E / std::pow(Eb + E, 3);
    return f_E;
}

std::pair<float, float> sampleThompsonDistribution(float Emax, float Eb) {
    float randomE = static_cast<float>(rand()) / RAND_MAX;
    float randomTheta = static_cast<float>(rand()) / RAND_MAX;

    float E = randomE * Emax;
    float theta = randomTheta * M_PI;

    return std::make_pair(E, theta);
}