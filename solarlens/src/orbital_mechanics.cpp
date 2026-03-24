#include "solarlens/physics/orbital_mechanics.hpp"
#include <algorithm>

namespace solarlens {

// ============================================================================
// Kepler's Equation Solver
// ============================================================================

double OrbitalMechanics::solve_kepler(double M, double e, double tol) {
    // Newton-Raphson iteration for E - e*sin(E) = M
    // Initial guess: E = M for small e, better starter for large e
    double E = (e < 0.8) ? M : constants::PI;

    for (int iter = 0; iter < 50; ++iter) {
        double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < tol) break;
    }
    return E;
}

// ============================================================================
// Orbital Element Conversions
// ============================================================================

void OrbitalMechanics::elements_to_state(const OrbitalElements& oe, Vec3& pos, Vec3& vel) {
    // Solve Kepler's equation for eccentric anomaly
    double E = solve_kepler(oe.M, oe.e);

    // True anomaly
    double cos_nu = (std::cos(E) - oe.e) / (1.0 - oe.e * std::cos(E));
    double sin_nu = (std::sqrt(1.0 - oe.e * oe.e) * std::sin(E)) / (1.0 - oe.e * std::cos(E));
    double nu = std::atan2(sin_nu, cos_nu);

    // Distance
    double r = oe.a * (1.0 - oe.e * std::cos(E));

    // Position and velocity in orbital plane
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);

    double p = oe.a * (1.0 - oe.e * oe.e);
    double h = std::sqrt(oe.mu * p);
    double vx_orb = -oe.mu / h * std::sin(nu);
    double vy_orb = oe.mu / h * (oe.e + std::cos(nu));

    // Rotation matrix elements (3-1-3: Omega, i, omega)
    double cO = std::cos(oe.Omega), sO = std::sin(oe.Omega);
    double ci = std::cos(oe.i),     si = std::sin(oe.i);
    double cw = std::cos(oe.omega), sw = std::sin(oe.omega);

    // Rotation from orbital plane to inertial frame
    double R11 = cO * cw - sO * sw * ci;
    double R12 = -cO * sw - sO * cw * ci;
    double R21 = sO * cw + cO * sw * ci;
    double R22 = -sO * sw + cO * cw * ci;
    double R31 = sw * si;
    double R32 = cw * si;

    pos.x = R11 * x_orb + R12 * y_orb;
    pos.y = R21 * x_orb + R22 * y_orb;
    pos.z = R31 * x_orb + R32 * y_orb;

    vel.x = R11 * vx_orb + R12 * vy_orb;
    vel.y = R21 * vx_orb + R22 * vy_orb;
    vel.z = R31 * vx_orb + R32 * vy_orb;
}

OrbitalElements OrbitalMechanics::state_to_elements(const Vec3& pos, const Vec3& vel, double mu) {
    OrbitalElements oe{};
    oe.mu = mu;

    double r = pos.norm();
    double v = vel.norm();

    // Specific angular momentum
    Vec3 h = pos.cross(vel);
    double h_mag = h.norm();

    // Node vector
    Vec3 n(-h.y, h.x, 0.0);
    double n_mag = n.norm();

    // Eccentricity vector
    Vec3 e_vec = (pos * (v * v - mu / r) - vel * pos.dot(vel)) / mu;
    oe.e = e_vec.norm();

    // Semi-major axis
    double energy = 0.5 * v * v - mu / r;
    if (std::abs(oe.e - 1.0) > 1e-10) {
        oe.a = -mu / (2.0 * energy);
    } else {
        oe.a = h_mag * h_mag / mu;  // Parabolic
    }

    // Inclination
    oe.i = std::acos(std::clamp(h.z / h_mag, -1.0, 1.0));

    // RAAN
    if (n_mag > 1e-10) {
        oe.Omega = std::acos(std::clamp(n.x / n_mag, -1.0, 1.0));
        if (n.y < 0) oe.Omega = constants::TWO_PI - oe.Omega;
    }

    // Argument of periapsis
    if (n_mag > 1e-10 && oe.e > 1e-10) {
        oe.omega = std::acos(std::clamp(n.dot(e_vec) / (n_mag * oe.e), -1.0, 1.0));
        if (e_vec.z < 0) oe.omega = constants::TWO_PI - oe.omega;
    }

    // True anomaly
    if (oe.e > 1e-10) {
        double cos_nu = std::clamp(e_vec.dot(pos) / (oe.e * r), -1.0, 1.0);
        double nu = std::acos(cos_nu);
        if (pos.dot(vel) < 0) nu = constants::TWO_PI - nu;

        // Convert true anomaly to mean anomaly
        double E = std::atan2(std::sqrt(1.0 - oe.e * oe.e) * std::sin(nu),
                              oe.e + std::cos(nu));
        oe.M = E - oe.e * std::sin(E);
        if (oe.M < 0) oe.M += constants::TWO_PI;
    }

    return oe;
}

OrbitalElements OrbitalMechanics::propagate_kepler(const OrbitalElements& oe, double dt_seconds) {
    OrbitalElements result = oe;

    // Mean motion
    double n = std::sqrt(oe.mu / (std::abs(oe.a) * oe.a * oe.a));
    result.M = std::fmod(oe.M + n * dt_seconds, constants::TWO_PI);
    if (result.M < 0) result.M += constants::TWO_PI;

    return result;
}

// ============================================================================
// Transfer Orbits
// ============================================================================

OrbitalMechanics::HohmannTransfer OrbitalMechanics::hohmann_transfer(double r1, double r2, double mu) {
    HohmannTransfer h{};

    // Transfer orbit semi-major axis
    h.a_transfer = (r1 + r2) / 2.0;

    // Velocities at r1 and r2 on circular orbits
    double v1_circ = std::sqrt(mu / r1);
    double v2_circ = std::sqrt(mu / r2);

    // Velocities on transfer orbit at r1 and r2
    double v1_transfer = std::sqrt(mu * (2.0 / r1 - 1.0 / h.a_transfer));
    double v2_transfer = std::sqrt(mu * (2.0 / r2 - 1.0 / h.a_transfer));

    // Delta-V
    h.dv1 = std::abs(v1_transfer - v1_circ);
    h.dv2 = std::abs(v2_circ - v2_transfer);
    h.dv_total = h.dv1 + h.dv2;

    // Transfer time (half orbital period of transfer orbit)
    h.transfer_time = constants::PI * std::sqrt(h.a_transfer * h.a_transfer * h.a_transfer / mu);

    return h;
}

OrbitalMechanics::BiEllipticTransfer OrbitalMechanics::bi_elliptic_transfer(
    double r1, double r2, double r_intermediate, double mu
) {
    BiEllipticTransfer b{};

    // First transfer ellipse: r1 to r_intermediate
    double a1 = (r1 + r_intermediate) / 2.0;
    double v1_circ = std::sqrt(mu / r1);
    double v1_trans = std::sqrt(mu * (2.0 / r1 - 1.0 / a1));
    b.dv1 = std::abs(v1_trans - v1_circ);

    // At r_intermediate: velocity on first ellipse
    double v_mid_1 = std::sqrt(mu * (2.0 / r_intermediate - 1.0 / a1));

    // Second transfer ellipse: r_intermediate to r2
    double a2 = (r_intermediate + r2) / 2.0;
    double v_mid_2 = std::sqrt(mu * (2.0 / r_intermediate - 1.0 / a2));
    b.dv2 = std::abs(v_mid_2 - v_mid_1);

    // At r2: velocity on second ellipse vs circular
    double v2_trans = std::sqrt(mu * (2.0 / r2 - 1.0 / a2));
    double v2_circ = std::sqrt(mu / r2);
    b.dv3 = std::abs(v2_circ - v2_trans);

    b.dv_total = b.dv1 + b.dv2 + b.dv3;
    b.transfer_time = constants::PI * (std::sqrt(a1 * a1 * a1 / mu) +
                                       std::sqrt(a2 * a2 * a2 / mu));

    return b;
}

// ============================================================================
// Gravity Assist
// ============================================================================

OrbitalMechanics::GravityAssist OrbitalMechanics::gravity_assist(
    const Vec3& v_inf_in,
    const Vec3& planet_vel,
    double periapsis,
    double mu_planet,
    double r_planet
) {
    GravityAssist result{};
    result.periapsis_r = periapsis;
    result.valid = periapsis > r_planet * 1.05;  // 5% margin above surface

    if (!result.valid) {
        return result;
    }

    // Incoming hyperbolic excess velocity
    double v_inf = v_inf_in.norm();

    // Turn angle from hyperbolic geometry
    // δ = 2 × arcsin(1 / e_hyp)  where e_hyp = 1 + r_p v_inf² / μ
    double e_hyp = 1.0 + periapsis * v_inf * v_inf / mu_planet;
    result.turn_angle = 2.0 * std::asin(1.0 / e_hyp);

    // Rotate incoming velocity by turn angle in the flyby plane
    // The flyby plane contains v_inf_in and the planet center
    // For simplicity, we rotate in the plane defined by v_inf_in and planet_vel
    Vec3 v_inf_hat = v_inf_in.normalized();
    Vec3 normal = v_inf_in.cross(planet_vel).normalized();
    // Handle case where v_inf is parallel to planet_vel
    if (normal.norm() < 1e-10) {
        normal = Vec3(0, 0, 1);
    }
    Vec3 perp = normal.cross(v_inf_hat).normalized();

    // Rotated outgoing v_infinity (relative to planet)
    double cos_delta = std::cos(result.turn_angle);
    double sin_delta = std::sin(result.turn_angle);
    Vec3 v_inf_out = v_inf_hat * cos_delta + perp * sin_delta;
    v_inf_out = v_inf_out * v_inf;

    // Heliocentric outgoing velocity
    result.v_out = planet_vel + v_inf_out;
    result.dv_gained = result.v_out.norm() - (v_inf_in + planet_vel).norm();

    return result;
}

// ============================================================================
// RK4 Integration
// ============================================================================

void OrbitalMechanics::rk4_step(Vec3& pos, Vec3& vel, double dt, double mu_sun) {
    // Acceleration function: a = -μ r / |r|³
    auto accel = [mu_sun](const Vec3& r) -> Vec3 {
        double r3 = r.norm() * r.norm() * r.norm();
        if (r3 < 1e-30) return Vec3{};
        return r * (-mu_sun / r3);
    };

    // RK4 stages
    Vec3 k1v = accel(pos);
    Vec3 k1r = vel;

    Vec3 k2v = accel(pos + k1r * (dt * 0.5));
    Vec3 k2r = vel + k1v * (dt * 0.5);

    Vec3 k3v = accel(pos + k2r * (dt * 0.5));
    Vec3 k3r = vel + k2v * (dt * 0.5);

    Vec3 k4v = accel(pos + k3r * dt);
    Vec3 k4r = vel + k3v * dt;

    pos = pos + (k1r + k2r * 2.0 + k3r * 2.0 + k4r) * (dt / 6.0);
    vel = vel + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (dt / 6.0);
}

// ============================================================================
// Solar Sail Trajectory
// ============================================================================

std::vector<OrbitalMechanics::SailState> OrbitalMechanics::integrate_solar_sail(
    const Vec3& pos0,
    const Vec3& vel0,
    double sail_area,
    double sail_mass,
    double cone_angle,
    double dt,
    double duration
) {
    std::vector<SailState> trajectory;

    const double mu_sun = constants::G * constants::M_SUN;
    Vec3 pos = pos0;
    Vec3 vel = vel0;
    double t = 0;

    // Characteristic acceleration at 1 AU:
    // a_c = 2 P_sun A cos²(α) / (c m)
    // where P_sun = L_sun / (4π r²) is radiation pressure
    // At 1 AU: P = 4.56 μN/m²
    const double cos_alpha = std::cos(cone_angle);

    // Reserve space for trajectory
    const uint32_t num_steps = static_cast<uint32_t>(duration / dt) + 1;
    trajectory.reserve(std::min(num_steps, static_cast<uint32_t>(100000)));

    uint32_t save_interval = std::max(num_steps / 10000, static_cast<uint32_t>(1));
    uint32_t step = 0;

    while (t < duration) {
        double r = pos.norm();
        if (r < 1e6) break;  // Too close to Sun

        // Solar radiation pressure force
        // F = (L_sun / (4π r² c)) × A × 2cos²α × r_hat
        // (factor of 2 for perfect reflection)
        double pressure = constants::L_SUN / (4.0 * constants::PI * r * r * constants::C);
        double force = pressure * sail_area * 2.0 * cos_alpha * cos_alpha;
        double sail_accel = force / sail_mass;

        // Gravitational acceleration
        double grav_accel = mu_sun / (r * r);

        // Net radial acceleration (sail pushes outward, gravity inward)
        Vec3 r_hat = pos.normalized();
        Vec3 accel = r_hat * (sail_accel - grav_accel);

        // Symplectic Euler integration (more stable for orbits)
        vel = vel + accel * dt;
        pos = pos + vel * dt;
        t += dt;
        step++;

        // Save state at intervals
        if (step % save_interval == 0) {
            SailState state;
            state.pos = pos;
            state.vel = vel;
            state.time = t;
            state.distance_au = r / constants::AU;
            state.speed_km_s = vel.norm() / 1000.0;
            trajectory.push_back(state);
        }
    }

    return trajectory;
}

// ============================================================================
// SGL Mission Planning
// ============================================================================

OrbitalMechanics::SGLTrajectory OrbitalMechanics::plan_sgl_mission(
    double target_distance_au,
    bool use_solar_sail,
    bool use_jupiter_assist
) {
    SGLTrajectory result{};
    const double mu_sun = constants::G * constants::M_SUN;

    // Jupiter parameters
    const double r_jupiter = 5.2 * constants::AU;
    const double mu_jupiter = 1.266865e17;  // GM_Jupiter (m³/s²)
    const double r_jupiter_body = 7.1492e7; // Jupiter radius (m)
    const double v_jupiter = circular_velocity(r_jupiter, mu_sun);

    if (use_jupiter_assist) {
        // Phase 1: Earth to Jupiter Hohmann transfer
        auto h = hohmann_transfer(constants::AU, r_jupiter, mu_sun);
        result.dv_earth_departure = h.dv1;

        // Phase 2: Jupiter gravity assist
        // Incoming v_infinity relative to Jupiter
        double v_arrive = std::sqrt(mu_sun * (2.0 / r_jupiter - 1.0 / h.a_transfer));
        double v_inf = std::abs(v_arrive - v_jupiter);

        Vec3 v_inf_vec(v_inf, 0, 0);
        Vec3 v_jup(0, v_jupiter, 0);
        auto ga = gravity_assist(v_inf_vec, v_jup, r_jupiter_body * 3.0, mu_jupiter, r_jupiter_body);

        result.dv_jupiter_assist = ga.dv_gained;
        double v_after_assist = ga.v_out.norm();

        // Phase 3: Coast to target distance
        // Time to reach target_distance_au from Jupiter at v_after_assist
        double coast_distance = (target_distance_au - 5.2) * constants::AU;
        double coast_time = coast_distance / v_after_assist;
        result.cruise_time_years = (h.transfer_time + coast_time) / (365.25 * 86400.0);
        result.arrival_speed_km_s = v_after_assist / 1000.0;
    }

    if (use_solar_sail) {
        // Solar sail from 1 AU outward
        // Start on circular orbit at 1 AU
        Vec3 pos0(constants::AU, 0, 0);
        Vec3 vel0(0, circular_velocity(constants::AU, mu_sun), 0);

        double dt = 86400.0;  // 1 day time step
        double max_time = 40.0 * 365.25 * 86400.0;  // 40 years max

        auto traj = integrate_solar_sail(
            pos0, vel0,
            constants::SOLAR_SAIL_AREA_M2,
            constants::SAT_MASS_KG,
            35.0 * constants::DEG2RAD,  // Optimal cone angle ~35°
            dt, max_time
        );

        result.trajectory = traj;

        // Find when we reach target distance
        for (const auto& state : traj) {
            if (state.distance_au >= target_distance_au) {
                result.cruise_time_years = state.time / (365.25 * 86400.0);
                result.arrival_speed_km_s = state.speed_km_s;
                break;
            }
        }

        if (result.cruise_time_years == 0 && !traj.empty()) {
            result.cruise_time_years = traj.back().time / (365.25 * 86400.0);
            result.arrival_speed_km_s = traj.back().speed_km_s;
        }
    }

    result.dv_course_correct = 200.0;  // 200 m/s budget for course corrections
    result.dv_total = result.dv_earth_departure + result.dv_course_correct;

    return result;
}

} // namespace solarlens
