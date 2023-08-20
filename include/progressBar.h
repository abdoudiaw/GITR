#include <iostream>
#include <iomanip>

class ProgressBar {
private:
    int total;
    int step;
    int barWidth;

public:
    ProgressBar(int total, int barWidth = 50) : total(total), step(0), barWidth(barWidth) {}

    void increment() {
        step++;
        int progress = static_cast<int>((static_cast<float>(step) / total) * barWidth);

        std::cout << "[";
        for (int i = 0; i < barWidth; i++) {
            if (i < progress) std::cout << "=";
            else if (i == progress) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::fixed << std::setprecision(2) << (static_cast<float>(step) / total) * 100.0 << "%\r";
        std::cout.flush();
    }

    void finish() {
        std::cout << std::endl;
    }
};
