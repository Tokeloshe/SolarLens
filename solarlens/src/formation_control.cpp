#include "solarlens/spacecraft/formation_control.hpp"
#include <cmath>

namespace solarlens {

bool SwarmController::optimize_formation(
    Formation type,
    const double target_vector[3],
    double baseline_km
) {
    switch (type) {
        case Formation::HEXAGONAL_GRID:
            return form_hexagonal_grid(target_vector, baseline_km);

        case Formation::LINEAR_ARRAY:
            return form_linear_array(target_vector, baseline_km);

        case Formation::EINSTEIN_RING:
            return form_einstein_ring(target_vector);

        default:
            return false;
    }
}

bool SwarmController::form_hexagonal_grid(const double target[3], double spacing_km) {
    // Hexagonal close packing for maximum coverage
    (void)target; // Reserved for target-oriented formation
    const double spacing_m = spacing_km * 1000.0;

    uint8_t sat_idx = 0;
    for (int ring = 0; ring <= 10 && sat_idx < active_count; ++ring) {
        const int sats_in_ring = (ring == 0) ? 1 : 6 * ring;

        for (int i = 0; i < sats_in_ring && sat_idx < active_count; ++i) {
            const double angle = 2.0 * M_PI * i / sats_in_ring;
            const double radius = ring * spacing_m;

            swarm[sat_idx].position[0] = constants::FOCAL_OPTIMAL_AU * constants::AU;
            swarm[sat_idx].position[1] = radius * std::cos(angle);
            swarm[sat_idx].position[2] = radius * std::sin(angle);

            // Check collision avoidance
            for (uint8_t j = 0; j < sat_idx; ++j) {
                const double dx = swarm[sat_idx].position[0] - swarm[j].position[0];
                const double dy = swarm[sat_idx].position[1] - swarm[j].position[1];
                const double dz = swarm[sat_idx].position[2] - swarm[j].position[2];
                const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (distance < constants::MIN_SEPARATION_M) {
                    return false;  // Collision risk
                }
            }

            sat_idx++;
        }
    }

    return true;
}

bool SwarmController::form_linear_array(const double target[3], double length_km) {
    // Linear interferometer array
    (void)target; // Reserved for target-oriented formation
    const double spacing_m = (length_km * 1000.0) / active_count;

    for (uint8_t i = 0; i < active_count; ++i) {
        swarm[i].position[0] = constants::FOCAL_OPTIMAL_AU * constants::AU;
        swarm[i].position[1] = (i - active_count/2.0) * spacing_m;
        swarm[i].position[2] = 0;
    }

    return true;
}

bool SwarmController::form_einstein_ring(const double target[3]) {
    // Form ring at Einstein radius for maximum magnification
    (void)target; // Reserved for target-oriented formation
    const double einstein_radius_km = 5000.0;  // Typical at 650 AU

    for (uint8_t i = 0; i < active_count; ++i) {
        const double angle = 2.0 * M_PI * i / active_count;
        swarm[i].position[0] = constants::FOCAL_OPTIMAL_AU * constants::AU;
        swarm[i].position[1] = einstein_radius_km * 1000.0 * std::cos(angle);
        swarm[i].position[2] = einstein_radius_km * 1000.0 * std::sin(angle);
    }

    return true;
}

} // namespace solarlens
