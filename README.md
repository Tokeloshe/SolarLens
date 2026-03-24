# SolarLens

**Mission control system for a Starlink-scale CubeSat swarm operating at the Solar Gravitational Lens focal region (~650 AU), enabling direct imaging and spectroscopic analysis of exoplanets.**

The Sun's gravity bends light from distant stars, creating a natural telescope with amplification of **~10^11** at optical wavelengths. SolarLens implements the complete software stack for a swarm of thousands of CubeSats that would fly to this focal region and work together as a distributed interferometric array to resolve exoplanet surfaces at kilometer-scale resolution from hundreds of light-years away.

---

## Architecture

```
solarlens/
├── include/solarlens/
│   ├── physics/
│   │   ├── constants.hpp           Physical constants (CODATA 2018, SI)
│   │   ├── gravitational_lens.hpp  SGL optics: PSF, gain, corona model
│   │   ├── orbital_mechanics.hpp   Kepler solver, transfers, solar sail
│   │   └── interferometry.hpp      UV-plane synthesis, CLEAN deconvolution
│   ├── spacecraft/
│   │   ├── formation_control.hpp   7 formation geometries, collision avoidance
│   │   ├── navigation.hpp          X-ray pulsar navigation (XNAV)
│   │   └── power_thermal.hpp       RTG + solar + battery, thermal equilibrium
│   ├── communication/
│   │   └── deep_space_relay.hpp    Reed-Solomon GF(2^8), link budgets, relay chain
│   └── imaging/
│       └── exoplanet_detector.hpp  Detection pipeline, biosignature scoring
├── src/                            Implementations (~3500 lines)
├── tests/                          Comprehensive test suite
└── CMakeLists.txt                  C++17, space-grade compile flags
```

## Key Features

### Solar Gravitational Lens Physics
- Turyshev & Toth (2017) wave-optical treatment
- Bessel J_0 PSF (not Gaussian) via Abramowitz & Stegun rational approximations
- Peak amplification: `mu = 4 pi^2 r_g / lambda` (~2.12 x 10^11 at 550 nm)
- Baumbach-Allen corona electron density and brightness model
- Plasma-corrected focal distance calculation
- Full SNR model with photon statistics, corona noise, dark current, read noise

### Orbital Mechanics & Trajectory
- Kepler equation solver (Newton-Raphson, <1e-12 convergence)
- Hohmann and bi-elliptic transfer orbit computation
- Gravity assist delta-V (hyperbolic flyby geometry)
- Solar sail trajectory integration (RK4)
- Complete Earth -> Jupiter assist -> SGL mission planning

### Starlink-Scale Swarm Formation
- 7 formation geometries: hexagonal grid, Golay sparse, ring array, logarithmic spiral, random uniform, VLA Y-array, Reuleaux triangle
- Collision detection and resolution
- Station keeping with PD control
- Formation quality metrics (max baseline, error RMS, health fraction)

### Aperture Synthesis Interferometry
- Baseline computation from satellite positions
- UV-plane coverage metrics (filling factor, max baseline)
- Direct Fourier Transform visibility measurement
- Hogbom CLEAN deconvolution algorithm
- Dirty beam/image generation

### Deep Space Communications
- Full GF(2^8) Reed-Solomon codec: RS(255, 223) with 16-symbol error correction
  - Berlekamp-Massey decoding, Chien search, Forney algorithm
- Friis equation link budget (Ka-band to DSN, laser inter-satellite, SGL-amplified)
- CCSDS-like frame encoding with sync patterns
- Relay chain planning with geometric spacing optimization

### X-ray Pulsar Navigation (XNAV)
- 8 real millisecond pulsars (PSR B1937+21, J0437-4715, etc.)
- Weighted least-squares position solver
- Sub-km positioning accuracy at 650 AU
- GDOP computation, 4x4 matrix inversion via Cramer's rule

### Power & Thermal Management
- RTG power decay (Pu-238 exponential, 87.7 yr half-life)
- Solar panel inverse-square scaling with radiation degradation
- Radiative thermal equilibrium (Stefan-Boltzmann)
- Battery state-of-charge tracking with charge/discharge efficiency
- Operational and survival temperature limit monitoring

### Exoplanet Detection Pipeline
- Richardson-Lucy deconvolution
- Matched-filter source detection with sub-pixel centroid
- Atmospheric composition analysis (O2, CH4, H2O, CO2, O3, NH3, SO2)
- Biosignature scoring (O2+CH4 thermodynamic disequilibrium)
- Chen & Kipping (2017) mass-radius relation
- Kopparapu et al. (2013) habitable zone boundaries
- Synthetic Earth-like spectrum and planet image generation

## Building

```bash
# Requirements: C++17 compiler, CMake 3.16+
cd solarlens
mkdir build && cd build
cmake ..
make -j$(nproc)

# Run the mission simulation demo
./solarlens

# Run the test suite
./solarlens_tests
```

## Sample Output

```
SOLARLENS — Solar Gravitational Lens Mission
Starlink-Scale CubeSat Swarm Control System

Physical Constants:
  Schwarzschild radius: 2954.01 m
  SGL focal minimum: 1097.3 AU
  SGL nominal distance: 650.0 AU

SGL peak amplification: 2.120e+11
Formation deployed: 16 satellites, max baseline 47 km

Exoplanet Detection:
  Temperature: 255 K
  Orbital radius: 1.00 AU (habitable zone: YES)
  Atmospheric O2: 39.4%, H2O: 6.7%, O3: 4.4%
  Biosignature Score: 55%

Reed-Solomon RS(255,223): 3 symbol errors corrected
Relay Chain: 7 nodes, 90 hour latency to Earth
```

## Physics References

1. **Turyshev & Toth (2017)** — "Diffraction of electromagnetic waves in the gravitational field of the Sun" — Wave-optical SGL PSF and gain formulas
2. **Turyshev et al. (2020)** — "Direct multipixel imaging and spectroscopy of an exoplanet with a solar gravitational lens mission" — Image reconstruction framework
3. **Allen (1973)** — Baumbach-Allen coronal electron density model
4. **Kopparapu et al. (2013)** — Habitable zone boundaries as function of stellar properties
5. **Chen & Kipping (2017)** — Probabilistic mass-radius forecaster for exoplanets
6. **Hogbom (1974)** — CLEAN deconvolution algorithm for radio interferometry

## License

MIT
