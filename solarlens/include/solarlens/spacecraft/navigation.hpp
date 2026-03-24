#pragma once

#include "solarlens/physics/constants.hpp"
#include "solarlens/physics/orbital_mechanics.hpp"
#include <array>
#include <cstdint>

namespace solarlens {

// ============================================================================
// DEEP SPACE NAVIGATION
//
// X-ray pulsar-based navigation (XNAV) using millisecond pulsars as
// natural GPS satellites. Uses weighted least-squares to solve for
// 3D position from pulse time-of-arrival measurements.
//
// Based on SEXTANT/NICER demonstration (2018) which achieved
// ~7 km accuracy in low Earth orbit. At 650 AU with longer integration
// times, sub-km accuracy is achievable.
// ============================================================================

class PulsarNavigator {
public:
    static constexpr uint32_t NUM_PULSARS = 8;

    struct Pulsar {
        const char* name;
        double ra_rad;            // Right ascension (J2000)
        double dec_rad;           // Declination (J2000)
        double period_s;          // Spin period (seconds)
        double period_dot;        // Period derivative (s/s)
        double distance_kpc;      // Distance (kiloparsecs)
        double flux_photons;      // X-ray flux (photons/cm²/s)
        // Direction unit vector (pre-computed)
        double dir[3];
    };

    PulsarNavigator();

    struct NavigationFix {
        Vec3 position_au;         // Position in AU (heliocentric ecliptic)
        Vec3 velocity_km_s;       // Velocity in km/s
        double position_error_km; // 1-sigma position uncertainty
        double time_error_ns;     // Clock offset
        uint8_t pulsars_used;     // Number of pulsars in solution
        double gdop;              // Geometric dilution of precision
        double chi_squared;       // Goodness of fit
        bool valid;               // Solution converged
    };

    // Compute navigation fix from pulse time-of-arrival residuals
    // toa_residuals_ns: measured - predicted pulse arrival times (nanoseconds)
    //                   for each pulsar (array of NUM_PULSARS)
    // weights: measurement weights (1/sigma² for each pulsar)
    // initial_guess: approximate position for linearization (AU)
    NavigationFix solve_position(
        const double toa_residuals_ns[NUM_PULSARS],
        const double weights[NUM_PULSARS],
        const Vec3& initial_guess
    ) const;

    // Get pulsar catalog
    const std::array<Pulsar, NUM_PULSARS>& catalog() const { return pulsars_; }

    // Predict pulse time-of-arrival at a given position
    // Returns predicted TOA relative to solar system barycenter (ns)
    double predict_toa_ns(uint32_t pulsar_idx, const Vec3& position_au, double epoch_s) const;

    // Compute GDOP (geometric dilution of precision) for current pulsar geometry
    // from a given position
    double compute_gdop(const Vec3& position_au) const;

private:
    std::array<Pulsar, NUM_PULSARS> pulsars_;

    // 4x4 matrix inversion for least-squares (position + clock offset)
    static bool invert_4x4(const double m[16], double inv[16]);
};

} // namespace solarlens
