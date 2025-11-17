#pragma once

#include <cstdint>

namespace solarlens {
namespace constants {

// ============================================================================
// FUNDAMENTAL PHYSICAL CONSTANTS (All SI units)
// ============================================================================

// Physical constants
constexpr double G = 6.67430e-11;              // Gravitational constant (m³/kg·s²)
constexpr double C = 299792458.0;              // Speed of light (m/s)
constexpr double H = 6.62607015e-34;           // Planck constant (J·s)
constexpr double K_B = 1.380649e-23;           // Boltzmann constant (J/K)

// Solar parameters
constexpr double M_SUN = 1.98847e30;           // Solar mass (kg)
constexpr double R_SUN = 6.95700e8;            // Solar radius (m)
constexpr double L_SUN = 3.828e26;             // Solar luminosity (W)
constexpr double T_SUN = 5778.0;               // Solar temperature (K)

// Distances
constexpr double AU = 1.495978707e11;          // Astronomical unit (m)
constexpr double LY = 9.4607304725808e15;      // Light year (m)
constexpr double PC = 3.0857e16;               // Parsec (m)

// Mission parameters
constexpr double FOCAL_MIN_AU = 547.8;         // Minimum focal distance
constexpr double FOCAL_OPTIMAL_AU = 650.0;     // Optimal for visible light
constexpr double FOCAL_MAX_AU = 900.0;         // Maximum useful distance

// Spacecraft constraints
constexpr uint32_t MAX_SWARM_SIZE = 256;       // Maximum satellites
constexpr double MIN_SEPARATION_M = 1000.0;    // Minimum satellite separation
constexpr double MAX_BASELINE_KM = 100000.0;   // Maximum array size

} // namespace constants
} // namespace solarlens
