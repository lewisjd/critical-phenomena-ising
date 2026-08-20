# Critical Phenomena in the Ising Model

A C++ Monte Carlo simulator for the Ising model on **2D square**, **3D simple-cubic (SC)**, and **3D body-centred-cubic (BCC)** lattices, implementing both the local **Metropolis** algorithm and the cluster-based **Swendsen–Wang** algorithm. It logs sweep-by-sweep thermodynamic observables to disk for downstream analysis of critical behaviour (magnetisation, energy, heat capacity, susceptibility, spin–spin correlations, and Binder cumulants).

This code underpins the accompanying report, [*Investigation of Critical Behaviour of the Ising Model*](IsingReport.pdf) — a Year 3 physics project (Bath, PH30056) that scored **96/100**.

## Physics background

The Ising Hamiltonian is

```
E = -J * Σ⟨ij⟩ σᵢσⱼ,   σᵢ ∈ {+1, -1}
```

summed over nearest-neighbour pairs `⟨ij⟩`, with `J = 1` (ferromagnetic) and `kB = 1`. Three lattice geometries are supported, differing in coordination number `z`:

| Lattice | Dimension | Coordination number `z` |
|---|---|---|
| Square | 2D | 4 |
| Simple cubic (SC) | 3D | 6 |
| Body-centred cubic (BCC) | 3D | 8 |

The 2D case has an exact solution (Onsager, 1944), used as a benchmark; the 3D cases do not, and their critical points are estimated numerically here via finite-size scaling of the Binder cumulant.

## Repository structure

```
critical-phenomena-ising/
├── MAIN_Source.cpp   # entry point: list of β values to scan, main simulation loop
├── IsingSystem.h      # IsingSystem class interface (grid, RNG, observables)
├── IsingSystem.cpp    # lattice setup, Metropolis + Swendsen–Wang updates, observable/CSV logging
├── rnd.h               # thin wrapper around <random> for seeding and sampling
└── IsingReport.pdf     # write-up: methodology, results, figures, references
```

## Algorithms

**Metropolis** (local, single-spin updates): for each site, propose a flip, accept it if `ΔE ≤ 0`, otherwise accept with probability `exp(-βΔE)`. Implemented in `MCsweep()`.

**Swendsen–Wang** (non-local, cluster updates): connect same-spin nearest-neighbour pairs with bond probability `p = 1 - exp(-2βJ)`, identify the resulting clusters with a union–find structure, then flip each cluster as a whole with probability 1/2. This suppresses critical slowing-down near the transition. Implemented in `SWsweep()`, and supports all three lattice geometries.

`IsingSystem::Update()` currently calls `MCsweep()` by default, with a call to `SWsweep()` present but commented out — swap which one is called (and recompile) to switch update schemes.

## Configuration

There's no command-line interface — everything is a compile-time constant, so changing settings means editing the source and rebuilding:

| Constant | File | Default | Purpose |
|---|---|---|---|
| `currentLattice` | `IsingSystem.h` | `LATTICE_2D` | selects `LATTICE_2D`, `LATTICE_3D_SC`, or `LATTICE_3D_BCC` |
| `gridSize` | `IsingSystem.h` | `500` | linear lattice size `L` (grid has `L^d` sites) |
| `maxSweeps` | `IsingSystem.h` | `1000` | sweeps recorded per `β` |
| `betas` | `MAIN_Source.cpp` | list of ~43 values | inverse temperatures to scan |

**Note:** `gridSize = 500` is sized for the 2D production runs described in the report. For 3D lattices this becomes `500³ ≈ 1.25×10⁸` sites, which is impractical — reduce `gridSize` (the report used `L = 16, 32, 64` for 3D, and `L = 16, 32, 64, 128` for the 2D finite-size scaling analysis) before switching `currentLattice`.

## Output format

`MCsweepdata.txt` is a CSV file with header:

```
beta,sweep,M,absM,M2,E,E2,acceptedFlips,attemptedClusters,acceptedClusters,fixedPairCorr,corr_r1,corr_r2,...,corr_r250
```

- `M`, `absM`, `M2` — magnetisation, its magnitude, and its square (per-spin values feed `χ`)
- `E`, `E2` — energy per spin and its square (feed the heat capacity `C`)
- `acceptedFlips` / `attemptedClusters` / `acceptedClusters` — update-scheme diagnostics
- `fixedPairCorr` — spin–spin correlation at separation `r = 1` between two fixed sites
- `corr_r*` — spin–spin correlation `⟨σᵢσⱼ⟩` at each separation in `{1,2,3,4,5,6,7,8,9,10,50,75,100,125,150,175,200,225,250}`

Susceptibility, heat capacity, and the Binder cumulant are derived from these columns in post-processing (not computed in this codebase):

```
χ  = (N/T)·(⟨m²⟩ - ⟨m⟩²)
C  = (N/T²)·(⟨e²⟩ - ⟨e⟩²)
U₄ = 1 - ⟨M⁴⟩ / (3⟨M²⟩²)
```

## Results summary

Using Swendsen–Wang sampling and Binder-cumulant finite-size scaling:

| Lattice | βc (this work) | Reference value |
|---|---|---|
| 2D square | 0.44051 ± 0.00035 | 0.44068 (Onsager, exact) |
| 3D simple cubic | 0.22429 ± 0.00049 | ≈0.22165 |
| 3D BCC | 0.15885 ± 0.00079 | ≈0.15737 |

Mean-field theory reproduces qualitative trends but overestimates `βc` and misses the non-analytic features at the transition — most severely in 2D, less so in 3D as coordination number increases. Full derivations, figures, and discussion are in `IsingReport.pdf`.

## Key references

- E. Ising, *Zeitschrift für Physik* **31**, 253 (1925)
- L. Onsager, *Physical Review* **65**, 117 (1944)
- N. Metropolis et al., *J. Chem. Phys.* **21**, 1087 (1953)
- R. H. Swendsen & J.-S. Wang, *Phys. Rev. Lett.* **58**, 86 (1987)
- K. Binder, *Z. Phys. B* **43**, 119 (1981)

(Full reference list in the report.)

## License

No license file is currently included. Add one (e.g. MIT) if you'd like to permit reuse.
