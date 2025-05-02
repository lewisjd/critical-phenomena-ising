#include <iostream>
#include <vector>
#include <ctime>
#include "IsingSystem.h"

using namespace std;

static const vector<double> betas = {
    0.02, 0.06, 0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38,
    0.40, 0.405, 0.410, 0.415, 0.420, 0.425, 0.430, 0.435, 0.440, 
    0.445, 0.450, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.50,
    0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.92, 0.96, 1.00
};

int main() {
    int seed = static_cast<int>(time(nullptr));
    cout << "Setting seed: " << seed << endl;
    IsingSystem sys;
    sys.setSeed(seed);

    for (int j = 0; j < 10; j++) {
        cout << "Starting loop " << j + 1 << "\n";
        for (size_t i = 0; i < betas.size(); ++i) {
            double beta = betas[i];
            cout << "Starting Beta = " << beta << endl;
            sys.Reset();
            sys.setTemperature(1.0 / beta);
            sys.resetSweepCount();

            int newSeed = static_cast<int>(time(nullptr)) + static_cast<int>(i);
            cout << "Seed for this run: " << newSeed << endl;
            sys.setSeed(newSeed);
            sys.setRunning();

            // Run the simulation for this beta
            while (sys.isRunning()) {
                sys.Update();
            }
            cout << "Finished Beta = " << beta << "\n";
        }
        cout << "Finished loop " << j + 1 << "\n";
        std::cout << "\a";
    }

    cout << "All betas complete. Exiting." << endl;
    return 0;
}
