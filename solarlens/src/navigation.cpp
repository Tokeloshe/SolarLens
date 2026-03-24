#include "solarlens/spacecraft/navigation.hpp"
#include <cmath>

namespace solarlens {

// ============================================================================
// Pulsar Catalog — Real millisecond pulsars used for XNAV
// ============================================================================

PulsarNavigator::PulsarNavigator() {
    // Initialize with real millisecond pulsar parameters
    // RA/Dec in radians (J2000), period in seconds
    pulsars_ = {{
        {"B1937+21",     5.0276,  0.3665,  0.001558,  1.05e-19, 3.6,  1.2e-4, {}},
        {"J1744-1134",   4.6141, -0.2046,  0.004075,  8.93e-21, 0.4,  3.5e-4, {}},
        {"J1909+3744",   5.0155,  0.6582,  0.002947,  1.40e-20, 1.1,  2.1e-4, {}},
        {"J0437-4715",   1.1662, -0.8288,  0.005757,  5.73e-20, 0.16, 5.5e-3, {}},
        {"J2124-3358",   5.5373, -0.5910,  0.004931,  2.06e-20, 0.3,  4.8e-4, {}},
        {"J0613-0200",   1.6330, -0.0349,  0.003062,  9.60e-21, 0.5,  1.8e-4, {}},
        {"J1012+5307",   2.6444,  0.9236,  0.005256,  1.71e-20, 0.5,  3.2e-4, {}},
        {"J1713+0747",   4.5094,  0.1346,  0.004570,  8.55e-21, 1.1,  2.5e-4, {}}
    }};

    // Pre-compute direction unit vectors
    for (auto& p : pulsars_) {
        double cos_dec = std::cos(p.dec_rad);
        p.dir[0] = std::cos(p.ra_rad) * cos_dec;
        p.dir[1] = std::sin(p.ra_rad) * cos_dec;
        p.dir[2] = std::sin(p.dec_rad);
    }
}

// ============================================================================
// Position Solution via Weighted Least Squares
// ============================================================================

PulsarNavigator::NavigationFix PulsarNavigator::solve_position(
    const double toa_residuals_ns[NUM_PULSARS],
    const double weights[NUM_PULSARS],
    const Vec3& initial_guess
) const {
    NavigationFix fix{};
    fix.valid = false;

    // Weighted Least Squares solution
    // Model: Δτ_i = (1/c) * (n_i · Δr) + Δt_clock
    // where n_i is direction to pulsar i, Δr is position correction, Δt_clock is clock offset
    //
    // This is a 4-parameter fit: Δx, Δy, Δz, Δt_clock
    // Design matrix H (NUM_PULSARS × 4):
    //   H[i] = [ n_ix/c, n_iy/c, n_iz/c, 1 ]
    // Observation: Δτ_i (TOA residuals in seconds)
    //
    // Solution: x = (H^T W H)^(-1) H^T W y

    // Build H^T W H (4×4) and H^T W y (4×1)
    double htwh[16] = {};  // 4×4 normal matrix
    double htwy[4] = {};   // 4×1 right-hand side

    uint8_t pulsars_used = 0;
    const double c_inv = 1.0 / constants::C;
    const double ns_to_s = 1e-9;
    const double au_to_m = constants::AU;

    for (uint32_t i = 0; i < NUM_PULSARS; ++i) {
        if (weights[i] <= 0) continue;

        double w = weights[i];
        // Direction cosines (already in cartesian)
        double h[4] = {
            pulsars_[i].dir[0] * c_inv * au_to_m,  // Convert AU offset to time delay
            pulsars_[i].dir[1] * c_inv * au_to_m,
            pulsars_[i].dir[2] * c_inv * au_to_m,
            1.0  // Clock offset (seconds)
        };

        double y = toa_residuals_ns[i] * ns_to_s;  // Convert ns to seconds

        // Accumulate normal equations
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                htwh[j * 4 + k] += w * h[j] * h[k];
            }
            htwy[j] += w * h[j] * y;
        }

        pulsars_used++;
    }

    if (pulsars_used < 4) {
        return fix;  // Need at least 4 pulsars for 3D + clock
    }

    // Solve via 4×4 matrix inversion
    double inv[16];
    if (!invert_4x4(htwh, inv)) {
        return fix;  // Singular matrix
    }

    // x = inv * htwy
    double dx[4] = {};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            dx[i] += inv[i * 4 + j] * htwy[j];
        }
    }

    // Apply correction to initial guess
    fix.position_au.x = initial_guess.x + dx[0];
    fix.position_au.y = initial_guess.y + dx[1];
    fix.position_au.z = initial_guess.z + dx[2];
    fix.time_error_ns = dx[3] * 1e9;  // Convert back to ns

    // Position error from covariance diagonal
    double pos_var = 0;
    for (int i = 0; i < 3; ++i) {
        pos_var += inv[i * 4 + i];
    }
    fix.position_error_km = std::sqrt(pos_var) * constants::AU / 1000.0;

    // GDOP from covariance trace
    double trace = 0;
    for (int i = 0; i < 4; ++i) {
        trace += inv[i * 4 + i];
    }
    fix.gdop = static_cast<float>(std::sqrt(trace));

    // Chi-squared
    fix.chi_squared = 0;
    for (uint32_t i = 0; i < NUM_PULSARS; ++i) {
        if (weights[i] <= 0) continue;
        double predicted = 0;
        for (int j = 0; j < 3; ++j) {
            predicted += pulsars_[i].dir[j] * c_inv * au_to_m * dx[j];
        }
        predicted += dx[3];
        double residual = toa_residuals_ns[i] * ns_to_s - predicted;
        fix.chi_squared += weights[i] * residual * residual;
    }

    fix.pulsars_used = pulsars_used;
    fix.valid = true;

    return fix;
}

// ============================================================================
// TOA Prediction
// ============================================================================

double PulsarNavigator::predict_toa_ns(uint32_t pulsar_idx, const Vec3& position_au, double epoch_s) const {
    if (pulsar_idx >= NUM_PULSARS) return 0;

    const auto& p = pulsars_[pulsar_idx];

    // Geometric delay: time for pulse to travel from barycenter to spacecraft
    // τ_geo = (1/c) * (n · r)  where n is direction to pulsar, r is position
    double delay_s = 0;
    double pos_m[3] = {position_au.x * constants::AU,
                       position_au.y * constants::AU,
                       position_au.z * constants::AU};

    for (int i = 0; i < 3; ++i) {
        delay_s += p.dir[i] * pos_m[i] / constants::C;
    }

    // Pulsar rotational phase at epoch
    // Phase = f * (t - t0) + 0.5 * f_dot * (t - t0)²
    // TOA of next pulse = when phase is integer
    double f = 1.0 / p.period_s;
    double f_dot = -p.period_dot / (p.period_s * p.period_s);

    double phase = f * epoch_s + 0.5 * f_dot * epoch_s * epoch_s;
    double next_pulse_phase = std::ceil(phase);
    double dt_to_pulse = (next_pulse_phase - phase) / f;

    return (delay_s + dt_to_pulse) * 1e9;  // Convert to nanoseconds
}

// ============================================================================
// GDOP Computation
// ============================================================================

double PulsarNavigator::compute_gdop(const Vec3& position_au) const {
    (void)position_au;  // GDOP depends on geometry, not position for distant pulsars

    // Build H^T H (4×4)
    double hth[16] = {};
    const double c_inv = 1.0 / constants::C;

    for (uint32_t i = 0; i < NUM_PULSARS; ++i) {
        double h[4] = {
            pulsars_[i].dir[0] * c_inv,
            pulsars_[i].dir[1] * c_inv,
            pulsars_[i].dir[2] * c_inv,
            1.0
        };
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                hth[j * 4 + k] += h[j] * h[k];
            }
        }
    }

    double inv[16];
    if (!invert_4x4(hth, inv)) return 999.0;

    double trace = 0;
    for (int i = 0; i < 4; ++i) trace += inv[i * 4 + i];
    return std::sqrt(std::abs(trace));
}

// ============================================================================
// 4×4 Matrix Inversion (Cramer's rule)
// ============================================================================

bool PulsarNavigator::invert_4x4(const double m[16], double inv[16]) {
    // Cofactor expansion
    double s0 = m[0]*m[5] - m[1]*m[4];
    double s1 = m[0]*m[6] - m[2]*m[4];
    double s2 = m[0]*m[7] - m[3]*m[4];
    double s3 = m[1]*m[6] - m[2]*m[5];
    double s4 = m[1]*m[7] - m[3]*m[5];
    double s5 = m[2]*m[7] - m[3]*m[6];

    double c5 = m[10]*m[15] - m[11]*m[14];
    double c4 = m[9]*m[15]  - m[11]*m[13];
    double c3 = m[9]*m[14]  - m[10]*m[13];
    double c2 = m[8]*m[15]  - m[11]*m[12];
    double c1 = m[8]*m[14]  - m[10]*m[12];
    double c0 = m[8]*m[13]  - m[9]*m[12];

    double det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0;
    if (std::abs(det) < 1e-30) return false;

    double invdet = 1.0 / det;

    inv[0]  = ( m[5]*c5 - m[6]*c4 + m[7]*c3) * invdet;
    inv[1]  = (-m[1]*c5 + m[2]*c4 - m[3]*c3) * invdet;
    inv[2]  = ( m[13]*s5 - m[14]*s4 + m[15]*s3) * invdet;
    inv[3]  = (-m[9]*s5 + m[10]*s4 - m[11]*s3) * invdet;

    inv[4]  = (-m[4]*c5 + m[6]*c2 - m[7]*c1) * invdet;
    inv[5]  = ( m[0]*c5 - m[2]*c2 + m[3]*c1) * invdet;
    inv[6]  = (-m[12]*s5 + m[14]*s2 - m[15]*s1) * invdet;
    inv[7]  = ( m[8]*s5 - m[10]*s2 + m[11]*s1) * invdet;

    inv[8]  = ( m[4]*c4 - m[5]*c2 + m[7]*c0) * invdet;
    inv[9]  = (-m[0]*c4 + m[1]*c2 - m[3]*c0) * invdet;
    inv[10] = ( m[12]*s4 - m[13]*s2 + m[15]*s0) * invdet;
    inv[11] = (-m[8]*s4 + m[9]*s2 - m[11]*s0) * invdet;

    inv[12] = (-m[4]*c3 + m[5]*c1 - m[6]*c0) * invdet;
    inv[13] = ( m[0]*c3 - m[1]*c1 + m[2]*c0) * invdet;
    inv[14] = (-m[12]*s3 + m[13]*s1 - m[14]*s0) * invdet;
    inv[15] = ( m[8]*s3 - m[9]*s1 + m[10]*s0) * invdet;

    return true;
}

} // namespace solarlens
