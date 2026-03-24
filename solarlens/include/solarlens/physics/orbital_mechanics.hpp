#pragma once

#include "constants.hpp"
#include <array>
#include <vector>
#include <cmath>

namespace solarlens {

// ============================================================================
// ORBITAL MECHANICS & TRAJECTORY CALCULATIONS
//
// Provides Keplerian orbit propagation, delta-V calculations,
// gravity assist modeling, and solar sail trajectory integration.
// ============================================================================

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3 operator/(double s) const { return {x/s, y/s, z/s}; }
    double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    Vec3 normalized() const { double n = norm(); return (n > 0) ? *this / n : Vec3{}; }
};

// Classical orbital elements
struct OrbitalElements {
    double a;      // Semi-major axis (m)
    double e;      // Eccentricity
    double i;      // Inclination (rad)
    double omega;  // Argument of periapsis (rad)
    double Omega;  // Longitude of ascending node (rad)
    double M;      // Mean anomaly (rad)
    double mu;     // Gravitational parameter GM (m³/s²)
};

class OrbitalMechanics {
public:
    // ---------- Keplerian Orbits ----------

    // Convert orbital elements to state vector (position, velocity)
    static void elements_to_state(const OrbitalElements& oe, Vec3& pos, Vec3& vel);

    // Convert state vector to orbital elements
    static OrbitalElements state_to_elements(const Vec3& pos, const Vec3& vel, double mu);

    // Propagate orbit by time dt using Kepler's equation
    static OrbitalElements propagate_kepler(const OrbitalElements& oe, double dt_seconds);

    // Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E
    // Uses Newton-Raphson iteration
    static double solve_kepler(double M, double e, double tol = 1e-12);

    // ---------- Transfer Orbits ----------

    // Hohmann transfer between circular orbits
    struct HohmannTransfer {
        double dv1;           // Departure delta-V (m/s)
        double dv2;           // Arrival delta-V (m/s)
        double dv_total;      // Total delta-V (m/s)
        double transfer_time; // Transfer time (seconds)
        double a_transfer;    // Transfer orbit semi-major axis (m)
    };

    static HohmannTransfer hohmann_transfer(double r1, double r2, double mu);

    // Bi-elliptic transfer (more efficient for r2/r1 > 11.94)
    struct BiEllipticTransfer {
        double dv1, dv2, dv3;
        double dv_total;
        double transfer_time;
    };

    static BiEllipticTransfer bi_elliptic_transfer(double r1, double r2, double r_intermediate, double mu);

    // ---------- Gravity Assist ----------

    struct GravityAssist {
        Vec3 v_out;           // Post-flyby velocity (heliocentric)
        double dv_gained;     // Speed gained (m/s)
        double periapsis_r;   // Closest approach distance (m)
        double turn_angle;    // Deflection angle (rad)
        bool valid;           // Whether periapsis is above surface
    };

    // Calculate gravity assist at a planet
    // v_inf_in: incoming hyperbolic excess velocity relative to planet
    // planet_vel: planet's heliocentric velocity
    // periapsis: closest approach distance from planet center (m)
    // mu_planet: planet's gravitational parameter (m³/s²)
    // r_planet: planet's radius (for validity check)
    static GravityAssist gravity_assist(
        const Vec3& v_inf_in,
        const Vec3& planet_vel,
        double periapsis,
        double mu_planet,
        double r_planet
    );

    // ---------- Solar Sail Trajectory ----------

    struct SailState {
        Vec3 pos;             // Position from Sun (m)
        Vec3 vel;             // Velocity (m/s)
        double time;          // Mission time (s)
        double distance_au;   // Distance from Sun (AU)
        double speed_km_s;    // Speed (km/s)
    };

    // Integrate solar sail trajectory from initial state
    // sail_area: sail area (m²)
    // sail_mass: total spacecraft mass (kg)
    // cone_angle: sail cone angle from radial (rad), 0 = face-on
    // dt: time step (s)
    // duration: total integration time (s)
    static std::vector<SailState> integrate_solar_sail(
        const Vec3& pos0,
        const Vec3& vel0,
        double sail_area,
        double sail_mass,
        double cone_angle,
        double dt,
        double duration
    );

    // ---------- Multi-body Propagation ----------

    // RK4 integration step for N-body problem (Sun + perturbations)
    static void rk4_step(Vec3& pos, Vec3& vel, double dt, double mu_sun);

    // ---------- Mission-specific ----------

    // Compute delta-V budget for Earth -> Jupiter gravity assist -> SGL focal line
    struct SGLTrajectory {
        double dv_earth_departure;  // Earth departure (m/s)
        double dv_jupiter_assist;   // Delta-V from Jupiter assist (m/s)
        double dv_course_correct;   // Mid-course corrections (m/s)
        double dv_total;            // Total delta-V (m/s)
        double cruise_time_years;   // Time to reach 650 AU
        double arrival_speed_km_s;  // Speed at SGL focal point
        std::vector<SailState> trajectory;  // Full trajectory
    };

    static SGLTrajectory plan_sgl_mission(
        double target_distance_au,
        bool use_solar_sail,
        bool use_jupiter_assist
    );

    // Orbital velocity at distance r from Sun
    static double circular_velocity(double r, double mu = constants::G * constants::M_SUN) {
        return std::sqrt(mu / r);
    }

    // Escape velocity at distance r from Sun
    static double escape_velocity(double r, double mu = constants::G * constants::M_SUN) {
        return std::sqrt(2.0 * mu / r);
    }
};

} // namespace solarlens
