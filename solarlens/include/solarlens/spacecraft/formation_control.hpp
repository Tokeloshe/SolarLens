#pragma once

#include "solarlens/physics/constants.hpp"
#include "solarlens/physics/orbital_mechanics.hpp"
#include <vector>
#include <cstdint>
#include <cmath>
#include <array>

namespace solarlens {

// ============================================================================
// STARLINK-SCALE SWARM FORMATION CONTROL
//
// Manages a constellation of thousands of mass-produced satellites,
// each with a 30cm telescope. The swarm operates as a distributed
// interferometric array at the solar gravitational lens focal region.
//
// Formations are oriented toward the target star system and optimized
// for UV-plane coverage in the aperture synthesis.
// ============================================================================

struct SatelliteState {
    uint32_t id;
    Vec3 position;            // Meters from swarm center
    Vec3 velocity;            // m/s relative to swarm center
    double quaternion[4];     // Attitude quaternion [w, x, y, z]
    float fuel_kg;
    float battery_wh;
    float temperature_k;
    uint8_t health;           // 0=dead, 1=degraded, 2=nominal, 3=optimal
    bool camera_active;
    bool laser_link_active;

    SatelliteState() : id(0), fuel_kg(0), battery_wh(0), temperature_k(0),
                       health(2), camera_active(true), laser_link_active(true) {
        quaternion[0] = 1.0; quaternion[1] = 0.0;
        quaternion[2] = 0.0; quaternion[3] = 0.0;
    }
};

class SwarmFormation {
public:
    // Formation geometries optimized for different science objectives
    enum class Geometry : uint8_t {
        HEXAGONAL_GRID,       // Dense UV coverage for imaging
        GOLAY_SPARSE,         // Golay-optimized sparse array (max UV coverage per sat)
        RING_ARRAY,           // Ring on Einstein ring for max gain sampling
        SPIRAL_ARM,           // Logarithmic spiral for multi-scale baselines
        RANDOM_UNIFORM,       // Random positions in disk (surprisingly good UV coverage)
        Y_SHAPED,             // VLA-style Y-array (good snapshot UV coverage)
        REULEAUX_TRIANGLE     // Optimal for uniform UV coverage
    };

    SwarmFormation();

    // Set the number of active satellites
    void set_swarm_size(uint32_t count);
    uint32_t get_swarm_size() const { return active_count_; }

    // Access satellite states
    const std::vector<SatelliteState>& satellites() const { return swarm_; }
    std::vector<SatelliteState>& satellites() { return swarm_; }

    // ---------- Formation Deployment ----------

    // Deploy formation oriented toward target direction
    // target_dir: unit vector toward target star
    // baseline_m: maximum extent of the array (m)
    // geometry: formation type
    bool deploy_formation(
        Geometry geometry,
        const Vec3& target_dir,
        double baseline_m
    );

    // ---------- Collision Avoidance ----------

    // Check all pairs for minimum separation violations
    // Returns number of violations found
    uint32_t check_collisions() const;

    // Resolve collision risks by adjusting violating satellites
    void resolve_collisions();

    // ---------- Station Keeping ----------

    // Compute delta-V commands to maintain formation against perturbations
    // drift_tolerance_m: maximum allowed position error
    // Returns total delta-V consumed across swarm (m/s)
    double station_keeping_step(double dt_seconds, double drift_tolerance_m);

    // ---------- Formation Quality Metrics ----------

    // Get 2D projected positions perpendicular to target direction
    std::vector<std::array<double, 2>> get_projected_positions(const Vec3& target_dir) const;

    // Maximum baseline length in current formation
    double max_baseline() const;

    // RMS of satellite distances from planned positions
    double formation_error_rms() const;

    // Fraction of satellites with health >= 2 (nominal)
    double health_fraction() const;

private:
    std::vector<SatelliteState> swarm_;
    std::vector<Vec3> planned_positions_;  // Target positions for each satellite
    uint32_t active_count_;

    // Formation generators
    void generate_hexagonal(double baseline_m);
    void generate_golay_sparse(double baseline_m);
    void generate_ring(double baseline_m);
    void generate_spiral(double baseline_m);
    void generate_random_uniform(double baseline_m);
    void generate_y_shaped(double baseline_m);
    void generate_reuleaux(double baseline_m);

    // Rotate all positions so z-axis aligns with target direction
    void orient_to_target(const Vec3& target_dir);

    // Simple LCG pseudo-random for deterministic formations
    uint32_t rng_state_;
    double random_uniform();
};

} // namespace solarlens
