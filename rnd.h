#pragma once

#include <random>

class rnd {
private:
    std::default_random_engine generator;
    int max;
    std::uniform_int_distribution<int>* intmax;
    std::uniform_real_distribution<double>* real01;

public:
    // constructor
    rnd() {
        max = 0x7fffffff;
        intmax = new std::uniform_int_distribution<int>(0, max);
        real01 = new std::uniform_real_distribution<double>(0.0, 1.0);
    }
    // destructor
    ~rnd() { delete intmax; delete real01; }

    // set the random seed
    void setSeed(int seed) { generator.seed(seed); }
    // random double in [0,1) and random integer in [0,max-1]
    double random01() { return (*real01)(generator); }
    int randomInt(int max) { return (*intmax)(generator) % max; }
};
