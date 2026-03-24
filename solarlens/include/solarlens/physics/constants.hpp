#pragma once

#include <cstdint>
#include <cmath>

namespace solarlens {
namespace constants {

// ============================================================================
// FUNDAMENTAL PHYSICAL CONSTANTS (CODATA 2018, SI units)
// ============================================================================

constexpr double G          = 6.67430e-11;          // Gravitational constant (m³ kg⁻¹ s⁻²)
constexpr double C          = 299792458.0;           // Speed of light (m/s)
constexpr double H          = 6.62607015e-34;        // Planck constant (J·s)
constexpr double H_BAR      = 1.054571817e-34;       // Reduced Planck constant (J·s)
constexpr double K_B        = 1.380649e-23;          // Boltzmann constant (J/K)
constexpr double SIGMA_SB   = 5.670374419e-8;        // Stefan-Boltzmann constant (W m⁻² K⁻⁴)
constexpr double E_CHARGE   = 1.602176634e-19;       // Elementary charge (C)
constexpr double M_ELECTRON = 9.1093837015e-31;      // Electron mass (kg)
constexpr double EPSILON_0  = 8.8541878128e-12;      // Vacuum permittivity (F/m)
constexpr double WIEN_B     = 2.897771955e-3;        // Wien displacement constant (m·K)

// ============================================================================
// SOLAR PARAMETERS
// ============================================================================

constexpr double M_SUN   = 1.98892e30;              // Solar mass (kg)
constexpr double R_SUN   = 6.9634e8;                // Solar radius (m) — IAU 2015
constexpr double L_SUN   = 3.828e26;                // Solar luminosity (W)
constexpr double T_SUN   = 5772.0;                  // Solar effective temperature (K) — IAU 2015

// Derived solar constants
constexpr double R_SCHWARZSCHILD = 2.0 * G * M_SUN / (C * C);  // ~2953.25 m

// ============================================================================
// ASTRONOMICAL DISTANCES
// ============================================================================

constexpr double AU = 1.495978707e11;               // Astronomical unit (m) — exact IAU
constexpr double LY = 9.4607304725808e15;           // Light year (m)
constexpr double PC = 3.0856775814913673e16;         // Parsec (m)

// ============================================================================
// SOLAR GRAVITATIONAL LENS PARAMETERS (Turyshev & Toth 2017, 2020)
// ============================================================================

// Minimum focal distance: z_0 = R_sun² / r_g  (r_g = Schwarzschild radius)
// z_0 = (6.9634e8)² / 2953.25 ≈ 1.641e14 m ≈ 547.6 AU
constexpr double SGL_FOCAL_MIN_M  = (R_SUN * R_SUN) / R_SCHWARZSCHILD;
constexpr double SGL_FOCAL_MIN_AU = SGL_FOCAL_MIN_M / AU;  // ~547.6 AU

// Operational range
constexpr double SGL_FOCAL_NOMINAL_AU = 650.0;      // Design point for visible light
constexpr double SGL_FOCAL_MAX_AU     = 900.0;      // Maximum useful distance

// SGL amplification for monochromatic point source:
// μ_max = 4π² r_g / λ  (dimensionless)
// At λ=550nm: μ ≈ 2.13 × 10^11
// We store r_g for runtime calculation
constexpr double R_G = R_SCHWARZSCHILD;  // Gravitational radius for gain formula

// Einstein ring angular radius at distance z from Sun:
// θ_E = sqrt(r_g / z)  [radians]
// Physical Einstein ring radius at observer plane:
// b_E = sqrt(r_g * z) = R_sun * sqrt(z / z_0)

// ============================================================================
// STARLINK-SCALE SWARM PARAMETERS
// ============================================================================

constexpr uint32_t MAX_SWARM_SIZE       = 10000;     // Starlink-scale constellation
constexpr uint32_t NOMINAL_SWARM_SIZE   = 4096;      // Nominal operational swarm
constexpr double   MIN_SEPARATION_M     = 100.0;     // Minimum inter-satellite distance (m)
constexpr double   MAX_BASELINE_M       = 1.0e8;     // 100,000 km max baseline
constexpr double   SAT_APERTURE_M       = 0.3;       // 30 cm telescope per satellite
constexpr double   SAT_MASS_KG          = 260.0;     // ~Starlink v2 mass
constexpr double   SAT_SOLAR_PANEL_M2   = 12.0;      // Solar panel area (inner solar system)
constexpr double   LASER_LINK_POWER_W   = 2.0;       // Inter-satellite laser link power
constexpr double   LASER_LINK_RATE_BPS  = 1.0e9;     // 1 Gbps inter-satellite link

// ============================================================================
// PROPULSION PARAMETERS
// ============================================================================

constexpr double ION_ISP_S           = 3000.0;       // Ion engine specific impulse (s)
constexpr double ION_THRUST_N        = 0.2;          // Ion engine thrust per sat (N)
constexpr double SOLAR_SAIL_SIGMA    = 5.0e-3;       // Sail area loading (kg/m²)
constexpr double SOLAR_SAIL_AREA_M2  = 1600.0;       // 40m × 40m sail
constexpr double SOLAR_PRESSURE_1AU  = 4.56e-6;      // Solar radiation pressure at 1 AU (N/m²)

// ============================================================================
// SENSOR PARAMETERS
// ============================================================================

constexpr double CCD_QUANTUM_EFF   = 0.90;           // Peak quantum efficiency
constexpr double CCD_READ_NOISE_E  = 3.0;            // Read noise (electrons RMS)
constexpr double CCD_DARK_CURRENT  = 0.001;          // Dark current at -100°C (e⁻/pix/s)
constexpr double CCD_PIXEL_SIZE_UM = 10.0;           // Pixel size (micrometers)
constexpr uint32_t CCD_PIXELS_X    = 4096;           // Detector width
constexpr uint32_t CCD_PIXELS_Y    = 4096;           // Detector height
constexpr double CCD_FULL_WELL     = 100000.0;       // Full well capacity (electrons)

// ============================================================================
// MATH HELPERS
// ============================================================================

constexpr double PI      = 3.14159265358979323846;
constexpr double TWO_PI  = 2.0 * PI;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double MAS2RAD = PI / (180.0 * 3600.0 * 1000.0);  // Milliarcsecond to radians

} // namespace constants
} // namespace solarlens
