#include <iostream>
#include <tuple>
#include <vector>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_errno.h>



std::tuple<gsl_interp2d*, gsl_interp2d*> build_interpolated_function(const std::vector<double>& energy_grid,
                                                                    const std::vector<double>& angle_grid,
                                                                    const std::vector<double>& sputter_table,
                                                                    const std::vector<double>& reflection_table)
{
    // Create GSL interpolation objects
    gsl_interp2d* sputter_interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, energy_grid.size(), angle_grid.size());
    gsl_interp2d* reflection_interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, energy_grid.size(), angle_grid.size());

    // Set GSL interpolation data
    gsl_interp2d_init(sputter_interp, energy_grid.data(), angle_grid.data(), sputter_table.data(), energy_grid.size(), angle_grid.size());
    gsl_interp2d_init(reflection_interp, energy_grid.data(), angle_grid.data(), reflection_table.data(), energy_grid.size(), angle_grid.size());

    return std::make_tuple(sputter_interp, reflection_interp);
}

int main()
{
    std::vector<double> energy_grid = {10., 13.05537869, 17.04429127, 22.25196771, 29.05078651,
                                       37.92690191, 49.51500669, 64.64371632, 84.39481966, 110.18063301,
                                       143.84498883, 187.79508019, 245.1735888, 320.08340466, 417.88100603,
                                       545.55947812, 712.24855849, 929.86746526, 1213.97718907, 1584.89319246};
    std::vector<double> angle_grid = {0., 11.2375, 22.475, 33.7125, 44.95, 56.1875, 67.425, 78.6625, 89.9};
    std::vector<double> sputter_table = {
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.200e-03, 6.900e-03, 2.100e-02,
        4.340e-02, 7.580e-02, 1.096e-01, 1.407e-01, 1.931e-01, 2.309e-01, 2.807e-01, 3.401e-01, 4.060e-01, 4.301e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.600e-03, 8.500e-03, 1.940e-02,
        4.350e-02, 7.500e-02, 1.102e-01, 1.509e-01, 1.972e-01, 2.344e-01, 3.013e-01, 3.487e-01, 4.029e-01, 4.464e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.000e-03, 6.300e-03, 2.220e-02,
        4.320e-02, 7.990e-02, 1.115e-01, 1.549e-01, 2.072e-01, 2.479e-01, 3.086e-01, 3.778e-01, 4.313e-01, 5.164e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.800e-03, 7.200e-03, 2.200e-02,
        4.870e-02, 7.460e-02, 1.231e-01, 1.774e-01, 2.250e-01, 2.929e-01, 3.660e-01, 4.313e-01, 5.178e-01, 5.982e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 6.000e-04, 5.900e-03, 1.780e-02,
        4.300e-02, 8.190e-02, 1.325e-01, 1.904e-01, 2.728e-01, 3.504e-01, 4.289e-01, 5.205e-01, 6.073e-01, 6.991e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 6.000e-04, 5.200e-03, 1.400e-02,
        4.100e-02, 8.200e-02, 1.375e-01, 2.141e-01, 3.120e-01, 4.079e-01, 5.105e-01, 6.254e-01, 7.324e-01, 8.675e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 1.000e-04, 1.400e-03, 9.800e-03,
        2.770e-02, 6.640e-02, 1.267e-01, 2.039e-01, 3.022e-01, 4.174e-01, 5.389e-01, 6.685e-01, 8.187e-01, 9.829e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 2.000e-04, 2.400e-03,
        1.100e-02, 2.480e-02, 5.440e-02, 9.030e-02, 1.399e-01, 2.031e-01, 2.978e-01, 4.033e-01, 5.559e-01, 6.878e-01,
        0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 2.000e-04,
        8.000e-04, 5.000e-04, 3.000e-04, 5.000e-04, 4.000e-04, 4.000e-04, 4.000e-04, 0.000e+00, 0.000e+00, 1.000e-04};

    std::vector<double> reflection_table = {
        0.539, 0.58, 0.5926, 0.5706, 0.5608, 0.5504, 0.5238, 0.5076, 0.4996, 0.4772,
        0.4581, 0.4517, 0.4373, 0.4237, 0.4162, 0.3967, 0.383, 0.3819, 0.3817, 0.3588,
        0.5473, 0.5873, 0.6065, 0.5847, 0.5755, 0.5548, 0.5387, 0.5177, 0.5032, 0.4867,
        0.4682, 0.4564, 0.4448, 0.4257, 0.4164, 0.4124, 0.3982, 0.3847, 0.3816, 0.3753,
        0.5838, 0.6177, 0.6298, 0.6133, 0.6081, 0.5879, 0.5676, 0.5525, 0.5271, 0.5206,
        0.491, 0.482, 0.4651, 0.4563, 0.4418, 0.4272, 0.4183, 0.4092, 0.3918, 0.3874,
        0.6234, 0.6534, 0.6674, 0.6618, 0.6512, 0.6322, 0.6219, 0.5902, 0.5839, 0.5604,
        0.5458, 0.5153, 0.504, 0.5012, 0.4814, 0.472, 0.4537, 0.4517, 0.4388, 0.422,
        0.6788, 0.72, 0.7346, 0.736, 0.7187, 0.7089, 0.6857, 0.6812, 0.6556, 0.6336,
        0.6159, 0.5888, 0.5757, 0.5586, 0.5344, 0.529, 0.5154, 0.4949, 0.49, 0.4793,
        0.7551, 0.79, 0.8032, 0.8094, 0.8068, 0.7944, 0.7821, 0.765, 0.7512, 0.7294,
        0.7139, 0.6866, 0.6652, 0.6538, 0.6254, 0.604, 0.5906, 0.5758, 0.5582, 0.543,
        0.8304, 0.861, 0.8889, 0.8959, 0.8994, 0.8934, 0.8859, 0.8831, 0.8695, 0.8581,
        0.8373, 0.8162, 0.794, 0.7708, 0.7469, 0.7204, 0.7072, 0.6764, 0.6681, 0.6433,
        0.8856, 0.9231, 0.9484, 0.9547, 0.964, 0.9676, 0.9679, 0.9705, 0.969, 0.9688,
        0.9659, 0.9576, 0.9497, 0.9405, 0.9303, 0.9118, 0.8928, 0.8779, 0.8469, 0.8226,
        0.904, 0.9423, 0.9575, 0.9777, 0.9852, 0.9895, 0.9917, 0.9963, 0.9977, 0.9987,
        0.999, 0.9997, 0.9999, 0.9998, 1., 0.9999, 0.9999, 1., 1., 1.};

    gsl_interp2d* sputter_interp;
    gsl_interp2d* reflection_interp;

    std::tie(sputter_interp, reflection_interp) = build_interpolated_function(energy_grid, angle_grid, sputter_table, reflection_table);

    double sputter_value = 0.0;
    double reflection_value = 0.0;
    double energy = 10.2;
    double angle = 1.5;

    if (energy >= energy_grid.front() && energy <= energy_grid.back() &&
        angle >= angle_grid.front() && angle <= angle_grid.back()) {
    sputter_value = gsl_interp2d_eval(sputter_interp, energy_grid.data(), angle_grid.data(), sputter_table.data(), energy, angle, nullptr, nullptr);
    reflection_value = gsl_interp2d_eval(reflection_interp, energy_grid.data(), angle_grid.data(), reflection_table.data(), energy, angle, nullptr, nullptr);
    }

    std::cout << "Interpolated Sputter Value: " << sputter_value << std::endl;
    std::cout << "Interpolated Reflection Value: " << reflection_value << std::endl;

    // Free GSL interpolation objects
    gsl_interp2d_free(sputter_interp);
    gsl_interp2d_free(reflection_interp);



std::cout << "Interpolated Sputter Value: " << sputter_value << std::endl;
std::cout << "Interpolated Reflection Value: " << reflection_value << std::endl;



    return 0;
}
