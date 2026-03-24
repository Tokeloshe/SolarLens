#include "solarlens/spacecraft/formation_control.hpp"
#include <cmath>
#include <algorithm>

namespace solarlens {

// ============================================================================
// SwarmFormation
// ============================================================================

SwarmFormation::SwarmFormation() : active_count_(0), rng_state_(42) {}

void SwarmFormation::set_swarm_size(uint32_t count) {
    active_count_ = std::min(count, constants::MAX_SWARM_SIZE);
    swarm_.resize(active_count_);
    planned_positions_.resize(active_count_);
    for (uint32_t i = 0; i < active_count_; ++i) {
        swarm_[i].id = i;
    }
}

double SwarmFormation::random_uniform() {
    // LCG: simple deterministic PRNG
    rng_state_ = rng_state_ * 1664525u + 1013904223u;
    return static_cast<double>(rng_state_) / 4294967296.0;
}

// ============================================================================
// Formation Deployment
// ============================================================================

bool SwarmFormation::deploy_formation(
    Geometry geometry,
    const Vec3& target_dir,
    double baseline_m
) {
    if (active_count_ == 0) return false;

    // Generate positions in the XY plane (will be rotated to face target)
    switch (geometry) {
        case Geometry::HEXAGONAL_GRID:    generate_hexagonal(baseline_m); break;
        case Geometry::GOLAY_SPARSE:      generate_golay_sparse(baseline_m); break;
        case Geometry::RING_ARRAY:        generate_ring(baseline_m); break;
        case Geometry::SPIRAL_ARM:        generate_spiral(baseline_m); break;
        case Geometry::RANDOM_UNIFORM:    generate_random_uniform(baseline_m); break;
        case Geometry::Y_SHAPED:          generate_y_shaped(baseline_m); break;
        case Geometry::REULEAUX_TRIANGLE: generate_reuleaux(baseline_m); break;
    }

    // Rotate the formation to face the target
    orient_to_target(target_dir);

    // Copy planned positions to actual positions
    for (uint32_t i = 0; i < active_count_; ++i) {
        swarm_[i].position = planned_positions_[i];
        swarm_[i].velocity = Vec3(0, 0, 0);
    }

    return check_collisions() == 0;
}

// ============================================================================
// Formation Generators
// ============================================================================

void SwarmFormation::generate_hexagonal(double baseline_m) {
    // Hexagonal close-packed grid filling a disk of diameter = baseline_m
    const double radius = baseline_m / 2.0;
    const double spacing = baseline_m / std::sqrt(static_cast<double>(active_count_) * 1.15);

    uint32_t idx = 0;
    for (int ring = 0; idx < active_count_; ++ring) {
        if (ring == 0) {
            planned_positions_[idx++] = Vec3(0, 0, 0);
            continue;
        }

        int n_in_ring = 6 * ring;
        for (int k = 0; k < n_in_ring && idx < active_count_; ++k) {
            double angle = constants::TWO_PI * k / n_in_ring;
            double r = ring * spacing;
            if (r > radius) r = radius;
            planned_positions_[idx++] = Vec3(r * std::cos(angle), r * std::sin(angle), 0);
        }
    }
}

void SwarmFormation::generate_golay_sparse(double baseline_m) {
    // Golay-inspired non-redundant array: maximize unique baselines
    // Uses a power-law radial distribution with golden-angle azimuthal spacing
    const double golden_angle = constants::PI * (3.0 - std::sqrt(5.0));  // ~137.5°
    const double max_r = baseline_m / 2.0;

    for (uint32_t i = 0; i < active_count_; ++i) {
        // Sunflower pattern: sqrt distribution for uniform area coverage
        double r = max_r * std::sqrt(static_cast<double>(i + 1) / active_count_);
        double theta = i * golden_angle;
        planned_positions_[i] = Vec3(r * std::cos(theta), r * std::sin(theta), 0);
    }
}

void SwarmFormation::generate_ring(double baseline_m) {
    // Satellites arranged on the Einstein ring radius
    double radius = baseline_m / 2.0;
    for (uint32_t i = 0; i < active_count_; ++i) {
        double angle = constants::TWO_PI * i / active_count_;
        planned_positions_[i] = Vec3(radius * std::cos(angle), radius * std::sin(angle), 0);
    }
}

void SwarmFormation::generate_spiral(double baseline_m) {
    // Logarithmic spiral arms — multi-scale baseline coverage
    const double max_r = baseline_m / 2.0;
    const uint32_t n_arms = 5;
    const double turns = 3.0;

    for (uint32_t i = 0; i < active_count_; ++i) {
        uint32_t arm = i % n_arms;
        uint32_t idx_in_arm = i / n_arms;
        uint32_t sats_per_arm = (active_count_ + n_arms - 1) / n_arms;

        double t = static_cast<double>(idx_in_arm) / sats_per_arm;
        double r = max_r * t;
        double theta = turns * constants::TWO_PI * t + arm * constants::TWO_PI / n_arms;

        planned_positions_[i] = Vec3(r * std::cos(theta), r * std::sin(theta), 0);
    }
}

void SwarmFormation::generate_random_uniform(double baseline_m) {
    // Random positions within a disk — surprisingly good UV coverage for large N
    double radius = baseline_m / 2.0;
    rng_state_ = 42;  // Deterministic seed

    for (uint32_t i = 0; i < active_count_; ++i) {
        // Uniform disk sampling: r = R * sqrt(U), θ = 2π V
        double r = radius * std::sqrt(random_uniform());
        double theta = constants::TWO_PI * random_uniform();
        planned_positions_[i] = Vec3(r * std::cos(theta), r * std::sin(theta), 0);
    }
}

void SwarmFormation::generate_y_shaped(double baseline_m) {
    // VLA-style Y-array: three arms at 120° angles
    double arm_length = baseline_m / 2.0;
    uint32_t per_arm = active_count_ / 3;
    uint32_t remainder = active_count_ % 3;

    uint32_t idx = 0;
    for (uint32_t arm = 0; arm < 3; ++arm) {
        double angle = arm * constants::TWO_PI / 3.0;
        uint32_t n = per_arm + (arm < remainder ? 1 : 0);

        for (uint32_t k = 0; k < n && idx < active_count_; ++k) {
            // Logarithmic spacing: more sats near center for short baselines
            double t = static_cast<double>(k + 1) / n;
            double r = arm_length * t * t;  // Quadratic spacing
            planned_positions_[idx++] = Vec3(r * std::cos(angle), r * std::sin(angle), 0);
        }
    }
}

void SwarmFormation::generate_reuleaux(double baseline_m) {
    // Reuleaux triangle: optimal for uniform UV coverage
    // Three circular arcs forming a curved triangle
    double R = baseline_m / 2.0;
    // The three vertices of the underlying equilateral triangle
    double tri_vertices[3][2];
    for (int k = 0; k < 3; ++k) {
        double angle = k * constants::TWO_PI / 3.0 - constants::PI / 2.0;
        tri_vertices[k][0] = R * std::cos(angle);
        tri_vertices[k][1] = R * std::sin(angle);
    }

    // Distribute satellites along the three arcs
    uint32_t per_arc = active_count_ / 3;
    uint32_t remainder = active_count_ % 3;
    uint32_t idx = 0;

    for (int arc = 0; arc < 3; ++arc) {
        // Each arc is centered on vertex[arc], from vertex[(arc+1)%3] to vertex[(arc+2)%3]
        double cx = tri_vertices[arc][0];
        double cy = tri_vertices[arc][1];

        double start_angle = std::atan2(tri_vertices[(arc + 1) % 3][1] - cy,
                                        tri_vertices[(arc + 1) % 3][0] - cx);
        double end_angle = std::atan2(tri_vertices[(arc + 2) % 3][1] - cy,
                                      tri_vertices[(arc + 2) % 3][0] - cx);

        // Ensure we go the short way around
        double arc_length_rad = end_angle - start_angle;
        if (arc_length_rad < 0) arc_length_rad += constants::TWO_PI;
        if (arc_length_rad > constants::PI) arc_length_rad -= constants::TWO_PI;

        double arc_radius = R * std::sqrt(3.0);  // Distance between vertices = R√3

        uint32_t n = per_arc + (static_cast<uint32_t>(arc) < remainder ? 1 : 0);
        for (uint32_t k = 0; k < n && idx < active_count_; ++k) {
            double t = static_cast<double>(k) / std::max(n - 1, 1u);
            double angle = start_angle + t * arc_length_rad;
            planned_positions_[idx++] = Vec3(
                cx + arc_radius * std::cos(angle),
                cy + arc_radius * std::sin(angle),
                0
            );
        }
    }
}

// ============================================================================
// Orientation
// ============================================================================

void SwarmFormation::orient_to_target(const Vec3& target_dir) {
    // Rotate all positions so that the formation's z-axis (normal to formation plane)
    // aligns with target_dir
    Vec3 z_hat = target_dir.normalized();
    if (z_hat.norm() < 0.5) z_hat = Vec3(1, 0, 0);  // Default to +x if zero

    // Construct orthonormal basis using Gram-Schmidt
    Vec3 up(0, 0, 1);
    if (std::abs(z_hat.z) > 0.99) up = Vec3(1, 0, 0);  // Avoid parallel

    Vec3 x_hat = up.cross(z_hat).normalized();
    Vec3 y_hat = z_hat.cross(x_hat).normalized();

    // Rotate each position: originally in XY plane (z=0)
    // New position = x_hat * old_x + y_hat * old_y + z_hat * old_z
    for (uint32_t i = 0; i < active_count_; ++i) {
        Vec3 old = planned_positions_[i];
        planned_positions_[i] = x_hat * old.x + y_hat * old.y + z_hat * old.z;
    }
}

// ============================================================================
// Collision Avoidance
// ============================================================================

uint32_t SwarmFormation::check_collisions() const {
    uint32_t violations = 0;
    for (uint32_t i = 0; i < active_count_; ++i) {
        for (uint32_t j = i + 1; j < active_count_; ++j) {
            Vec3 diff = swarm_[i].position - swarm_[j].position;
            if (diff.norm() < constants::MIN_SEPARATION_M) {
                ++violations;
            }
        }
    }
    return violations;
}

void SwarmFormation::resolve_collisions() {
    // Push apart any satellites that are too close
    for (uint32_t iter = 0; iter < 100; ++iter) {
        bool any_violation = false;
        for (uint32_t i = 0; i < active_count_; ++i) {
            for (uint32_t j = i + 1; j < active_count_; ++j) {
                Vec3 diff = swarm_[j].position - swarm_[i].position;
                double dist = diff.norm();
                if (dist < constants::MIN_SEPARATION_M && dist > 1e-6) {
                    Vec3 push = diff.normalized() * ((constants::MIN_SEPARATION_M - dist) * 0.6);
                    swarm_[i].position = swarm_[i].position - push;
                    swarm_[j].position = swarm_[j].position + push;
                    any_violation = true;
                }
            }
        }
        if (!any_violation) break;
    }
}

// ============================================================================
// Station Keeping
// ============================================================================

double SwarmFormation::station_keeping_step(double dt_seconds, double drift_tolerance_m) {
    double total_dv = 0;

    for (uint32_t i = 0; i < active_count_; ++i) {
        Vec3 error = planned_positions_[i] - swarm_[i].position;
        double err_mag = error.norm();

        if (err_mag > drift_tolerance_m) {
            // PD controller: apply correction thrust
            double kp = 0.001;  // Position gain
            double kd = 0.1;    // Damping gain
            Vec3 correction = error * kp - swarm_[i].velocity * kd;

            double dv = correction.norm() * dt_seconds;
            total_dv += dv;

            swarm_[i].velocity = swarm_[i].velocity + correction * dt_seconds;
            swarm_[i].position = swarm_[i].position + swarm_[i].velocity * dt_seconds;

            // Fuel consumption: Tsiolkovsky with ion engine ISP
            double exhaust_vel = constants::ION_ISP_S * 9.80665;
            double mass_ratio = 1.0 - std::exp(-dv / exhaust_vel);
            swarm_[i].fuel_kg -= static_cast<float>(swarm_[i].fuel_kg * mass_ratio);
        } else {
            swarm_[i].position = swarm_[i].position + swarm_[i].velocity * dt_seconds;
        }
    }

    return total_dv;
}

// ============================================================================
// Metrics
// ============================================================================

std::vector<std::array<double, 2>> SwarmFormation::get_projected_positions(const Vec3& target_dir) const {
    std::vector<std::array<double, 2>> projected(active_count_);

    Vec3 z_hat = target_dir.normalized();
    Vec3 up(0, 0, 1);
    if (std::abs(z_hat.z) > 0.99) up = Vec3(1, 0, 0);
    Vec3 x_hat = up.cross(z_hat).normalized();
    Vec3 y_hat = z_hat.cross(x_hat).normalized();

    for (uint32_t i = 0; i < active_count_; ++i) {
        projected[i][0] = swarm_[i].position.dot(x_hat);
        projected[i][1] = swarm_[i].position.dot(y_hat);
    }

    return projected;
}

double SwarmFormation::max_baseline() const {
    double max_b = 0;
    for (uint32_t i = 0; i < active_count_; ++i) {
        for (uint32_t j = i + 1; j < active_count_; ++j) {
            double d = (swarm_[i].position - swarm_[j].position).norm();
            if (d > max_b) max_b = d;
        }
    }
    return max_b;
}

double SwarmFormation::formation_error_rms() const {
    if (active_count_ == 0) return 0;
    double sum_sq = 0;
    for (uint32_t i = 0; i < active_count_; ++i) {
        Vec3 err = swarm_[i].position - planned_positions_[i];
        sum_sq += err.dot(err);
    }
    return std::sqrt(sum_sq / active_count_);
}

double SwarmFormation::health_fraction() const {
    if (active_count_ == 0) return 0;
    uint32_t healthy = 0;
    for (uint32_t i = 0; i < active_count_; ++i) {
        if (swarm_[i].health >= 2) ++healthy;
    }
    return static_cast<double>(healthy) / active_count_;
}

} // namespace solarlens
