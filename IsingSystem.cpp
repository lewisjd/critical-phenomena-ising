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

    // Allocate grid vector for the chosen lattice type
    int N;
    if (currentLattice == LATTICE_2D) {
        N = gridSize * gridSize;
    }
    else {
        N = gridSize * gridSize * gridSize;
    }
    grid.resize(N);

    // initialize temperature and grid
    Reset();

    logfile.open("MCsweepdata.txt");
    std::vector<int> rs = {
        1,2,3,4,5,6,7,8,9,10,
        50,75,100,125,150,175,200,225,250
    };

    // write header for output file
    logfile << "beta,sweep,M,absM,M2,E,E2,acceptedFlips,attemptedClusters,acceptedClusters,fixedPairCorr";
    for (int d : rs) {
        logfile << ",corr_r" << d;
    }
    logfile << "\n";
}

IsingSystem::~IsingSystem() {
    if (logfile.is_open())
        logfile.close();
    // grid (vector) is automatically deallocated
}

// reset temperature and re-initialize grid
void IsingSystem::Reset() {
    double initialTemp = 4.0;
    setTemperature(initialTemp);
    // set all spins to -1
    int N = grid.size();
    for (int i = 0; i < N; ++i)
        grid[i] = -1;
}

// draw the system (no-op after removing graphics)
void IsingSystem::DrawSquares() {
    // Drawing code removed: no graphical output.
}

// calculate magnetisation per spin
double IsingSystem::calculateMagnetisation() {
    long long totalSpin = 0;
    if (currentLattice == LATTICE_2D) {
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                int idx = i * gridSize + j;
                totalSpin += grid[idx];
            }
        }
        return static_cast<double>(totalSpin) / (gridSize * gridSize);
    }
    else {
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    totalSpin += grid[idx];
                }
            }
        }
        return static_cast<double>(totalSpin) / (gridSize * gridSize * gridSize);
    }
}

// energy per spin (dimensionless E/(N*J))
double IsingSystem::calculateEnergy() {
    double energy = 0.0;
    if (currentLattice == LATTICE_2D) {
        // sum over 2D nearest neighbors (right and down)
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                int idx = i * gridSize + j;
                int s = grid[idx];
                // right neighbor
                int idxR = i * gridSize + ((j + 1) % gridSize);
                energy += -s * grid[idxR];
                // down neighbor
                int idxD = ((i + 1) % gridSize) * gridSize + j;
                energy += -s * grid[idxD];
            }
        }
        return energy / (gridSize * gridSize);
    }
    else if (currentLattice == LATTICE_3D_SC) {
        // simple cubic: sum over positive directions (x, y, z)
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    int s = grid[idx];
                    // x-direction
                    int in = (i + 1) % gridSize;
                    int idxX = (in * gridSize + j) * gridSize + k;
                    energy += -s * grid[idxX];
                    // y-direction
                    int jn = (j + 1) % gridSize;
                    int idxY = (i * gridSize + jn) * gridSize + k;
                    energy += -s * grid[idxY];
                    // z-direction
                    int kn = (k + 1) % gridSize;
                    int idxZ = (i * gridSize + j) * gridSize + kn;
                    energy += -s * grid[idxZ];
                }
            }
        }
        return energy / (gridSize * gridSize * gridSize);
    }
    else { // LATTICE_3D_BCC
        // body-centered cubic: sum over diagonal neighbors (8 directions)
        // we use offsets with dx=1 and all combinations of dy,dz = ?1
        const int bccOffsets[4][3] = {
            { 0,  0,  1},
            { 1,  0,  1},
            { 0, -1,  1},
            { 1, -1,  1}
        };

        double energy = 0.0;

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    int s = grid[idx];

                    for (int m = 0; m < 4; ++m) {
                        int ni = (i + bccOffsets[m][0]) % gridSize;
                        int nj = (j + bccOffsets[m][1] + gridSize) % gridSize;
                        int nk = (k + bccOffsets[m][2]) % gridSize;

                        int nidx = (ni * gridSize + nj) * gridSize + nk;
                        energy += -s * grid[nidx];
                    }
                }
            }
        }

        return energy / (gridSize * gridSize * gridSize);
    }
}

// one Swendsen?Wang sweep
void IsingSystem::SWsweep() {
    acceptedFlips = 0;
    attemptedClusters = 0;
    acceptedClusters = 0;
    int N;
    if (currentLattice == LATTICE_2D) {
        N = gridSize * gridSize;
    }
    else {
        N = gridSize * gridSize * gridSize;
    }
    std::vector<int> parent(N);
    for (int i = 0; i < N; ++i) parent[i] = i;

    std::function<int(int)> find = [&](int x) -> int {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
        };
    auto unite = [&](int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
        };

    // bond probability for J=1: p = 1 - exp(-2*beta)
    double p = 1.0 - std::exp(-2.0 * inverseTemperatureBeta);

    // build clusters by adding bonds between equal spins
    if (currentLattice == LATTICE_2D) {
        // 2D: bonds to the right and down neighbors
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                int idx = i * gridSize + j;
                int s = grid[idx];
                // right neighbor
                int ni = i;
                int nj = (j + 1) % gridSize;
                int idxR = ni * gridSize + nj;
                if (grid[idxR] == s && rgen.random01() < p)
                    unite(idx, idxR);
                // down neighbor
                ni = (i + 1) % gridSize;
                nj = j;
                int idxD = ni * gridSize + nj;
                if (grid[idxD] == s && rgen.random01() < p)
                    unite(idx, idxD);
            }
        }
    }
    else if (currentLattice == LATTICE_3D_SC) {
        // SC: bonds along +x, +y, +z
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    int s = grid[idx];
                    // +x neighbor
                    {
                        int ni = (i + 1) % gridSize;
                        int nj = j;
                        int nk = k;
                        int idxX = (ni * gridSize + nj) * gridSize + nk;
                        if (grid[idxX] == s && rgen.random01() < p)
                            unite(idx, idxX);
                    }
                    // +y neighbor
                    {
                        int ni = i;
                        int nj = (j + 1) % gridSize;
                        int nk = k;
                        int idxY = (ni * gridSize + nj) * gridSize + nk;
                        if (grid[idxY] == s && rgen.random01() < p)
                            unite(idx, idxY);
                    }
                    // +z neighbor
                    {
                        int ni = i;
                        int nj = j;
                        int nk = (k + 1) % gridSize;
                        int idxZ = (ni * gridSize + nj) * gridSize + nk;
                        if (grid[idxZ] == s && rgen.random01() < p)
                            unite(idx, idxZ);
                    }
                }
            }
        }
    }
    else { // LATTICE_3D_BCC
        // BCC: bonds along diagonal offsets (dx=1, dy,dz=?1)
        const int bccOffsets[4][3] = {
            { 0,  0,  1},
            { 1,  0,  1},
            { 0, -1,  1},
            { 1, -1,  1}
        };

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    int s = grid[idx];

                    for (int m = 0; m < 4; ++m) {
                        int ni = (i + bccOffsets[m][0]) % gridSize;
                        int nj = (j + bccOffsets[m][1] + gridSize) % gridSize;
                        int nk = (k + bccOffsets[m][2]) % gridSize;

                        int nidx = (ni * gridSize + nj) * gridSize + nk;
                        if (grid[nidx] == s && rgen.random01() < p)
                            unite(idx, nidx);
                    }
                }
            }
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
    if (currentLattice == LATTICE_2D) {
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                int idx = i * gridSize + j;
                if (flipCluster[find(idx)])
                    grid[idx] = -grid[idx];
            }
        }
    }
    else {
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    if (flipCluster[find(idx)])
                        grid[idx] = -grid[idx];
                }
            }
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

    if (currentLattice == LATTICE_2D) {
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                int idx = i * gridSize + j;
                for (int k = 0; k < nPoints; ++k) {
                    int d = rs[k];
                    // x-direction
                    int idxX = ((i + d) % gridSize) * gridSize + j;
                    G[k] += grid[idx] * grid[idxX];
                    counts[k]++;
                    // y-direction
                    int idxY = i * gridSize + ((j + d) % gridSize);
                    G[k] += grid[idx] * grid[idxY];
                    counts[k]++;
                }
            }
        }
    }
    else {
        // 3D: include z-direction
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                for (int k = 0; k < gridSize; ++k) {
                    int idx = (i * gridSize + j) * gridSize + k;
                    for (int m = 0; m < nPoints; ++m) {
                        int d = rs[m];
                        // x-direction
                        int idxX = (((i + d) % gridSize) * gridSize + j) * gridSize + k;
                        G[m] += grid[idx] * grid[idxX];
                        counts[m]++;
                        // y-direction
                        int idxY = (i * gridSize + (j + d) % gridSize) * gridSize + k;
                        G[m] += grid[idx] * grid[idxY];
                        counts[m]++;
                        // z-direction
                        int idxZ = (i * gridSize + j) * gridSize + ((k + d) % gridSize);
                        G[m] += grid[idx] * grid[idxZ];
                        counts[m]++;
                    }
                }
            }
        }
    }

    // normalize
    for (int k = 0; k < nPoints; ++k) {
        if (counts[k] > 0) G[k] /= counts[k];
    }
    return G;
}

// set a grid cell
void IsingSystem::setGrid(int pos[], int val) {
    if (currentLattice == LATTICE_2D) {
        int idx = pos[0] * gridSize + pos[1];
        grid[idx] = val;
    }
    else {
        int idx = (pos[0] * gridSize + pos[1]) * gridSize + pos[2];
        grid[idx] = val;
    }
}

// read a grid cell
int IsingSystem::readGrid(int pos[]) {
    if (currentLattice == LATTICE_2D) {
        return grid[pos[0] * gridSize + pos[1]];
    }
    else {
        int idx = (pos[0] * gridSize + pos[1]) * gridSize + pos[2];
        return grid[idx];
    }
}

// flip a spin
void IsingSystem::flipSpin(int pos[]) {
    if (currentLattice == LATTICE_2D) {
        int idx = pos[0] * gridSize + pos[1];
        grid[idx] = -grid[idx];
    }
    else {
        int idx = (pos[0] * gridSize + pos[1]) * gridSize + pos[2];
        grid[idx] = -grid[idx];
    }
}

// set setpos to the position of a neighbour of pos
// For 2D: val = 0..3 for +x, -x, +y, -y
// For 3D SC: val = 0..5 for +x,-x,+y,-y,+z,-z
// For 3D BCC: val = 0..7 for the 8 diagonal neighbors
void IsingSystem::setPosNeighbour(int setpos[], int pos[], int val) {
    if (currentLattice == LATTICE_2D) {
        switch (val) {
        case 0: // +x
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = 0;
            break;
        case 1: // -x
            setpos[0] = (pos[0] - 1 + gridSize) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = 0;
            break;
        case 2: // +y
            setpos[0] = pos[0];
            setpos[1] = (pos[1] + 1) % gridSize;
            setpos[2] = 0;
            break;
        case 3: // -y
            setpos[0] = pos[0];
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = 0;
            break;
        }
    }
    else if (currentLattice == LATTICE_3D_SC) {
        switch (val) {
        case 0: // +x
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = pos[2];
            break;
        case 1: // -x
            setpos[0] = (pos[0] - 1 + gridSize) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = pos[2];
            break;
        case 2: // +y
            setpos[0] = pos[0];
            setpos[1] = (pos[1] + 1) % gridSize;
            setpos[2] = pos[2];
            break;
        case 3: // -y
            setpos[0] = pos[0];
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = pos[2];
            break;
        case 4: // +z
            setpos[0] = pos[0];
            setpos[1] = pos[1];
            setpos[2] = (pos[2] + 1) % gridSize;
            break;
        case 5: // -z
            setpos[0] = pos[0];
            setpos[1] = pos[1];
            setpos[2] = (pos[2] - 1 + gridSize) % gridSize;
            break;
        }
    }
    else { // LATTICE_3D_BCC
        // diagonal offsets for BCC neighbors
        switch (val) {
        case 0: // 0 0 +1
            setpos[0] = pos[0];
            setpos[1] = pos[1];
            setpos[2] = (pos[2] + 1) % gridSize;
            break;
        case 1: // +1 0 +1
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = (pos[2] + 1) % gridSize;
            break;
        case 2: // 0 -1 +1 
            setpos[0] = pos[0];
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = (pos[2] + 1) % gridSize;
            break;
        case 3: // +1 -1 +1
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = (pos[2] + 1) % gridSize;
            break;
        case 4: // 0 0 -1
            setpos[0] = pos[0];
            setpos[1] = pos[1];
            setpos[2] = (pos[2] - 1 + gridSize) % gridSize;
            break;
        case 5: // +1 0 -1
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = pos[1];
            setpos[2] = (pos[2] - 1 + gridSize) % gridSize;
            break;
        case 6: // 0 -1 -1
            setpos[0] = pos[0];
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = (pos[2] - 1 + gridSize) % gridSize;
            break;
        case 7: // +1 -1 -1
            setpos[0] = (pos[0] + 1) % gridSize;
            setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
            setpos[2] = (pos[2] - 1 + gridSize) % gridSize;
            break;
        }
    }
}

// update the system: one Swendsen-Wang sweep, log data, check sweeps
void IsingSystem::Update() {
    MCsweep();
    //SWsweep();
    ++sweepCount;

    double M = calculateMagnetisation();
    double E = calculateEnergy();
    auto corr = calculateCorrelation();
    int p0[3] = { 0, 0, 0 };
    int p1[3] = { 1, 0, 0 };

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

// one Metropolis sweep (2D only)
void IsingSystem::MCsweep() {
    acceptedFlips = 0;
    int N = gridSize * gridSize;     // 2D only
    for (int idx = 0; idx < N; ++idx) {
        // pick spin at random
        int i = idx / gridSize, j = idx % gridSize;
        int s = grid[idx];
        
        // sum over 4 neighbours (periodic BC)
        int up = grid[((i - 1 + gridSize) % gridSize) * gridSize + j];
        int down = grid[((i + 1) % gridSize) * gridSize + j];
        int left = grid[i * gridSize + (j - 1 + gridSize) % gridSize];
        int right = grid[i * gridSize + (j + 1) % gridSize];
        
        int deltaE = 2 * s * (up + down + left + right);
        if (deltaE <= 0 || rgen.random01() < std::exp(-deltaE * inverseTemperatureBeta)) {
            grid[idx] = -s;
            ++acceptedFlips;
            
        }
    }
}