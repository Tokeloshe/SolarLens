#pragma once

#include "solarlens/physics/constants.hpp"
#include <array>
#include <cstdint>
#include <cmath>

namespace solarlens {

// ============================================================================
// SWARM FORMATION CONTROL
// ============================================================================
class SwarmController {
private:
    struct Spacecraft {
        uint8_t id;
        double position[3];      // m from sun
        double velocity[3];      // m/s
        double quaternion[4];     // Attitude
        float fuel_kg;
        float battery_wh;
        float temperature_k;
        uint8_t status;          // Bit flags
    };

    std::array<Spacecraft, constants::MAX_SWARM_SIZE> swarm;
    uint8_t active_count;

public:
    SwarmController() : active_count(0) {}

    // Formation types for different observations
    enum class Formation : uint8_t {
        HEXAGONAL_GRID,     // For imaging
        LINEAR_ARRAY,       // For interferometry
        CIRCULAR_RING,      // For coronagraphy
        DISPERSED_CLOUD,    // For wide field
        EINSTEIN_RING       // For maximum magnification
    };

    // Optimize formation with collision avoidance
    [[nodiscard]] bool optimize_formation(
        Formation type,
        const double target_vector[3],
        double baseline_km
    );

    void set_active_count(uint8_t count) {
        active_count = (count < constants::MAX_SWARM_SIZE) ? count : constants::MAX_SWARM_SIZE;
    }

private:
    bool form_hexagonal_grid(const double target[3], double spacing_km);
    bool form_linear_array(const double target[3], double length_km);
    bool form_einstein_ring(const double target[3]);
};

} // namespace solarlens
