#pragma once

#include "solarlens/physics/constants.hpp"
#include <array>
#include <cstdint>

namespace solarlens {

// ============================================================================
// DEEP SPACE NAVIGATION
// ============================================================================
class DeepSpaceNavigator {
private:
    // Pulsar timing for navigation (like GPS in space)
    struct Pulsar {
        double ra_rad;           // Right ascension
        double dec_rad;          // Declination
        double period_ms;        // Millisecond period
        double period_derivative; // Spin-down rate
        uint32_t last_pulse_time;
    };

    // Known millisecond pulsars for navigation
    std::array<Pulsar, 6> navigation_pulsars;

public:
    DeepSpaceNavigator();

    struct NavigationSolution {
        double position_au[3];   // Position in AU from Sun
        double velocity_km_s[3]; // Velocity in km/s
        double position_error_km;
        double time_error_ns;
        uint8_t pulsars_used;
        float gdop;              // Geometric dilution of precision
    };

    // Calculate position from pulsar timing (X-ray navigation)
    [[nodiscard]] NavigationSolution calculate_position(
        const uint32_t pulse_times[6],
        uint64_t current_time_ns
    );
};

} // namespace solarlens
