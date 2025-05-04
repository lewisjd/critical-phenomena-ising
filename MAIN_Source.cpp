#include <iostream>
#include <vector>
#include <ctime>
#include "IsingSystem.h"

using namespace std;

static const vector<double> betas = {
    /*0.02 , 0.04 , 0.06 , 0.08 , 0.1  , 0.12 , 0.14 , 0.16 , 0.18 ,
    0.185, 0.19 , 0.195, 0.2  , 0.205, 0.21 , 0.215, 0.22 , 0.225,
    0.23 , 0.235, 0.24 , 0.245, 0.25 , 0.255, 0.26 , 0.28 , 0.3  ,          //SC
    0.32 , 0.34 , 0.36 , 0.38 , 0.4  , 0.42 , 0.44 , 0.46 , 0.48 ,
    0.5  , 0.52 , 0.54 , 0.56 , 0.58 , 0.6*/

    0.02, 0.04, 0.06, 0.08, 0.10, 0.12,
    0.125, 0.13, 0.135, 0.14, 0.145, 0.15,
    0.155, 0.16, 0.165, 0.17, 0.175, 0.18,
    0.185, 0.19, 0.195, 0.20,
    0.22, 0.24, 0.26, 0.28, 0.30, 0.32,         //BCC
    0.34, 0.36, 0.38, 0.40, 0.42, 0.44,
    0.46, 0.48, 0.50, 0.52, 0.54, 0.56,
    0.58, 0.60,
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
