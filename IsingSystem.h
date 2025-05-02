#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <sstream>

#include "rnd.h"

using namespace std;

// IsingSystem class: this is where the Ising simulation happens
class IsingSystem {
private:
    // size of grid: fixed constant
    static const int gridSize = 128;
    int** grid;  // 2D array storing each site's spin (+1 or -1)

    rnd rgen;    // random number generator

    int isActive;       // 1 if system is running, 0 if paused
    int slowNotFast;    // 1 if in "slow" mode, 0 in fast mode
    double inverseTemperatureBeta; // inverse temperature (beta)

    ofstream logfile;   // output file for logging

    int sweepCount = 0;
    static constexpr int maxSweeps = 150;

    int acceptedFlips = 0;  // count of accepted single-spin flips per sweep
    std::vector<double> calculateCorrelation();

    int attemptedClusters = 0; // clusters found in SW sweep
    int acceptedClusters = 0;  // clusters flipped in SW sweep

public:
    // constructor
    IsingSystem();
    // destructor
    ~IsingSystem();

    // set the random seed
    void setSeed(int s) { rgen.setSeed(s); }

    // functions to update the system
    double computeLocalField(int pos[]);
    void attemptSpinFlip();
    void MCsweep();

    // adjust temperature
    void Hotter() { inverseTemperatureBeta -= 0.05; }
    void Colder() { inverseTemperatureBeta += 0.05; }

    void setTemperature(double TT) { inverseTemperatureBeta = 1.0 / TT; }

    // access grid entries via functions
    void setGrid(int pos[], int val);
    int readGrid(int pos[]);
    void flipSpin(int pos[]);

    // reset temperature and re-initialize grid
    void Reset();

    // interface for isActive
    int isRunning() { return isActive; }
    void setRunning() { isActive = 1; }
    void pauseRunning() { isActive = 0; }

    // interface for slowNotFast
    int isSlow() { return slowNotFast; }
    void setSlow() { slowNotFast = 1; }
    void setFast() { slowNotFast = 0; }

    // set setpos to the position of a neighbour of pos
    // val = 0,1,2,3 determines which neighbour (with periodic BC)
    void setPosNeighbour(int setpos[], int pos[], int val);

    // perform a Monte Carlo or Swendsen-Wang sweep
    void Update();

    // draw the system (no-op after removing graphics)
    void DrawSquares();

    double calculateMagnetisation();
    double calculateEnergy();

    void SWsweep(); // one Swendsen-Wang sweep

    void resetSweepCount() { sweepCount = 0; }
};
