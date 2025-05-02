#include "IsingSystem.h"
#include <vector>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <unordered_set>

// constructor
IsingSystem::IsingSystem() {
    cout << "creating system, gridSize " << gridSize << endl;
    inverseTemperatureBeta = 1.0 / 4.0;
    slowNotFast = 1;
    isActive = 0;

    // Allocate memory for the grid
    grid = new int* [gridSize];
    for (int i = 0; i < gridSize; i++) {
        grid[i] = new int[gridSize];
    }

    // initialize temperature and grid
    Reset();

    logfile.open("SWsweepdataL128.txt");
    std::vector<int> rs = { 1,2,3,4,5,6,7,8,9,10,
                           50,75,100,125,150,175,200,225,250 };

    // write header for output file
    logfile << "beta,sweep,M,absM,M2,E,E2,acceptedFlips,attemptedClusters,acceptedClusters,fixedPairCorr";
    for (int d : rs) {
        logfile << ",corr_r" << d;
    }
    logfile << "\n";
}

void IsingSystem::Reset() {
    double initialTemp = 4.0;
    setTemperature(initialTemp);
    // set all spins to -1
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int pos[2] = { i, j };
            setGrid(pos, -1);
        }
    }
}

// destructor
IsingSystem::~IsingSystem() {
    if (logfile.is_open())
        logfile.close();

    // delete the grid
    for (int i = 0; i < gridSize; i++)
        delete[] grid[i];
    delete[] grid;
}

// this draws the system (no-op after removing graphics)
void IsingSystem::DrawSquares() {
    // Drawing code removed: no graphical output.
}

// attempt N spin flips (N = gridSize*gridSize)
void IsingSystem::MCsweep() {
    acceptedFlips = 0;
    attemptedClusters = 0;  // no clusters in MC
    acceptedClusters = 0;
    for (int i = 0; i < gridSize * gridSize; i++) {
        attemptSpinFlip();
    }
}

// attempt to flip one spin with Metropolis criterion
void IsingSystem::attemptSpinFlip() {
    int pos[2];
    // pick random site
    pos[0] = rgen.randomInt(gridSize);
    pos[1] = rgen.randomInt(gridSize);

    double hloc = computeLocalField(pos);
    double dE = 2.0 * hloc * readGrid(pos);
    if (dE < 0) {
        flipSpin(pos);
        ++acceptedFlips;
    }
    else if (rgen.random01() < exp(-dE)) {
        flipSpin(pos);
        ++acceptedFlips;
    }
}

// compute local field * (Beta)
double IsingSystem::computeLocalField(int pos[]) {
    double result = 0.0;
    for (int i = 0; i < 4; i++) {
        int nborPos[2];
        setPosNeighbour(nborPos, pos, i);
        result += readGrid(nborPos);
    }
    result *= inverseTemperatureBeta;
    return result;
}

// set a grid cell
void IsingSystem::setGrid(int pos[], int val) {
    grid[pos[0]][pos[1]] = val;
}

// read a grid cell
int IsingSystem::readGrid(int pos[]) {
    return grid[pos[0]][pos[1]];
}

// flip a spin
void IsingSystem::flipSpin(int pos[]) {
    grid[pos[0]][pos[1]] = -grid[pos[0]][pos[1]];
}

// set setpos to the neighbour of pos (val = 0,1,2,3 for +x, -x, +y, -y)
void IsingSystem::setPosNeighbour(int setpos[], int pos[], int val) {
    switch (val) {
    case 0:
        setpos[0] = (pos[0] + 1) % gridSize;
        setpos[1] = pos[1];
        break;
    case 1:
        setpos[0] = (pos[0] - 1 + gridSize) % gridSize;
        setpos[1] = pos[1];
        break;
    case 2:
        setpos[0] = pos[0];
        setpos[1] = (pos[1] + 1) % gridSize;
        break;
    case 3:
        setpos[0] = pos[0];
        setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
        break;
    }
}

// update the system: one Monte Carlo or SW sweep, log data, check sweeps
void IsingSystem::Update() {
    SWsweep();
    //MCsweep();

    ++sweepCount;

    double M = calculateMagnetisation();
    double E = calculateEnergy();
    auto corr = calculateCorrelation();
    int p0[2] = { 0, 0 };
    int p1[2] = { 1, 0 };

    // write one line to log file
    logfile << inverseTemperatureBeta << ','
        << sweepCount << ','
        << M << ','
        << std::abs(M) << ','
        << M * M << ','
        << E << ','
        << E * E << ','
        << acceptedFlips << ','
        << attemptedClusters << ','
        << acceptedClusters << ','
        << (readGrid(p0) * readGrid(p1));

    for (size_t r = 0; r < corr.size(); ++r) {
        logfile << ',' << corr[r];
    }
    logfile << '\n';

    cout << "Sweep " << sweepCount << "/" << maxSweeps << "\n";
    if (sweepCount >= maxSweeps) {
        pauseRunning();
    }
}

// magnetisation per spin
double IsingSystem::calculateMagnetisation() {
    int totalSpin = 0;
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int pos[2] = { i, j };
            totalSpin += readGrid(pos);
        }
    }
    return static_cast<double>(totalSpin) / (gridSize * gridSize);
}

// energy per spin (dimensionless E/(N*J))
double IsingSystem::calculateEnergy() {
    double energy = 0.0;
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int pos[2] = { i, j };
            int neighborPos[2];
            // right neighbor
            setPosNeighbour(neighborPos, pos, 2);
            energy += -readGrid(pos) * readGrid(neighborPos);
            // down neighbor
            setPosNeighbour(neighborPos, pos, 0);
            energy += -readGrid(pos) * readGrid(neighborPos);
        }
    }
    return energy / (gridSize * gridSize);
}

// one Swendsen–Wang sweep
void IsingSystem::SWsweep() {
    attemptedClusters = 0;
    acceptedClusters = 0;
    const int N = gridSize * gridSize;
    std::vector<int> parent(N);
    for (int i = 0; i < N; ++i) parent[i] = i;

    std::function<int(int)> find = [&](int x) -> int {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
        };
    auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
        };

    // bond probability for J=1: p = 1 - exp(-2?J)
    double p = 1.0 - std::exp(-2.0 * inverseTemperatureBeta);

    // build clusters
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int idx = i * gridSize + j;
            int s = grid[i][j];
            // right neighbor
            int ni = i;
            int nj = (j + 1) % gridSize;
            if (grid[ni][nj] == s && rgen.random01() < p)
                unite(idx, ni * gridSize + nj);
            // down neighbor
            ni = (i + 1) % gridSize;
            nj = j;
            if (grid[ni][nj] == s && rgen.random01() < p)
                unite(idx, ni * gridSize + nj);
        }
    }

    // decide for each cluster whether to flip it
    std::unordered_map<int, bool> flipCluster;
    flipCluster.reserve(N);
    for (int i = 0; i < N; ++i) {
        int r = find(i);
        if (!flipCluster.count(r))
            flipCluster[r] = (rgen.random01() < 0.5);
    }

    // apply cluster flips
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int idx = i * gridSize + j;
            if (flipCluster[find(idx)])
                grid[i][j] = -grid[i][j];
        }
    }

    // count distinct clusters
    {
        std::unordered_set<int> roots;
        roots.reserve(N);
        for (int i = 0; i < N; ++i)
            roots.insert(find(i));
        attemptedClusters = static_cast<int>(roots.size());
    }

    // count how many clusters flipped
    for (auto& kv : flipCluster) {
        if (kv.second)
            ++acceptedClusters;
    }
}

// calculate correlation at set distances
std::vector<double> IsingSystem::calculateCorrelation() {
    std::vector<int> rs = {
        1,2,3,4,5,6,7,8,9,10,
        50,75,100,125,150,175,200,225,250
    };

    int nPoints = rs.size();
    std::vector<double> G(nPoints, 0.0);
    std::vector<int> counts(nPoints, 0);

    int pos[2], nbr[2];
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            pos[0] = i; pos[1] = j;
            for (int k = 0; k < nPoints; ++k) {
                int d = rs[k];
                // x-direction
                nbr[0] = (i + d) % gridSize;  nbr[1] = j;
                G[k] += readGrid(pos) * readGrid(nbr);
                counts[k]++;
                // y-direction
                nbr[0] = i; nbr[1] = (j + d) % gridSize;
                G[k] += readGrid(pos) * readGrid(nbr);
                counts[k]++;
            }
        }
    }
    // normalize
    for (int k = 0; k < nPoints; ++k) {
        if (counts[k] > 0) G[k] /= counts[k];
    }
    return G;
}
