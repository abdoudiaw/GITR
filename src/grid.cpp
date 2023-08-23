#include <iostream>
#include <vector>
#include <cmath>
#include "Particles.h"


class Grid1D {
private:
    double x_min, x_max;
    int n_cells;
    std::vector<int> cells;

public:
    Grid1D(double xmin, double xmax, int n) : x_min(xmin), x_max(xmax), n_cells(n), cells(n, 0) {}

    int getCellIndex(double x) {
        int i = floor((x - x_min) / (x_max - x_min) * n_cells);
        return std::min(std::max(i, 0), n_cells - 1);
    }

    void addParticle(const Particle& p) {
        int idx = getCellIndex(p.x);
        cells[idx]++;
    }

    void printGrid() {
        for (int i = 0; i < n_cells; i++) {
            std::cout << "Cell " << i << ": " << cells[i] << " particles" << std::endl;
        }
    }
};
