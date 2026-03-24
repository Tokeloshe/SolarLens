#include "solarlens/physics/gravitational_lens.hpp"
#include <cmath>
#include <algorithm>

namespace solarlens {

// ============================================================================
// Bessel functions J0 and J1 — rational approximations (Abramowitz & Stegun)
// Accurate to ~1e-8 over all x
// ============================================================================

double GravitationalLens::bessel_j0(double x) {
    x = std::abs(x);
    if (x < 8.0) {
        double y = x * x;
        double num = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                    + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        double den = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
        return num / den;
    } else {
        double z = 8.0 / x;
        double y = z * z;
        double xx = x - 0.785398164;
        double p0 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                   + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        double q0 = -0.1562499995e-1 + y * (0.1430488765e-3
                   + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                   + y * (-0.934935152e-7))));
        return std::sqrt(0.636619772 / x) * (p0 * std::cos(xx) - z * q0 * std::sin(xx));
    }
}

double GravitationalLens::bessel_j1(double x) {
    double sign = (x < 0) ? -1.0 : 1.0;
    x = std::abs(x);
    if (x < 8.0) {
        double y = x * x;
        double num = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                    + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
        double den = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                    + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
        return sign * num / den;
    } else {
        double z = 8.0 / x;
        double y = z * z;
        double xx = x - 2.356194491;
        double p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                   + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
        double q1 = 0.04687499995 + y * (-0.2002690873e-3
                   + y * (0.8449199096e-5 + y * (-0.88228987e-6
                   + y * 0.105787412e-6)));
        return sign * std::sqrt(0.636619772 / x) * (p1 * std::cos(xx) - z * q1 * std::sin(xx));
    }
}

// ============================================================================
// Constructor
// ============================================================================

GravitationalLens::GravitationalLens()
    : r_g_(constants::R_SCHWARZSCHILD)
    , z0_(constants::SGL_FOCAL_MIN_M)
    , z0_au_(constants::SGL_FOCAL_MIN_AU)
{}

// ============================================================================
// Focal Geometry
// ============================================================================

double GravitationalLens::focal_distance_au(double wavelength_nm) const {
    // Base geometric focal distance (achromatic): z_0 = R_sun² / r_g
    // Plasma correction: the solar corona refracts light, adding a
    // wavelength-dependent shift. The plasma frequency at the solar limb
    // (r = R_sun) with electron density ~10^8 cm^-3 is:
    // f_p = 8.98 kHz * sqrt(n_e) ≈ 8.98e3 * 1e4 ≈ 90 MHz
    //
    // The corrected focal distance:
    // z_min(λ) = z_0 / (1 - (f_p/f)²)
    // For optical wavelengths, (f_p/f)² ~ 10^-13, so correction is tiny.

    const double wavelength_m = wavelength_nm * 1e-9;
    const double light_freq = constants::C / wavelength_m;

    // Electron density at the solar limb (cm^-3) — Baumbach-Allen at r=1.0
    const double n_e_limb = 1.0e8;  // cm^-3
    // Plasma frequency
    const double f_plasma = 8.98e3 * std::sqrt(n_e_limb);  // Hz

    const double ratio_sq = (f_plasma * f_plasma) / (light_freq * light_freq);
    const double z_corrected = z0_ / (1.0 - ratio_sq);

    return z_corrected / constants::AU;
}

double GravitationalLens::einstein_angle_rad(double observer_distance_au) const {
    // θ_E = sqrt(r_g / z)  where z is observer distance from Sun
    const double z = observer_distance_au * constants::AU;
    return std::sqrt(r_g_ / z);
}

double GravitationalLens::einstein_ring_radius_m(double observer_distance_au) const {
    // b_E = sqrt(r_g * z) = R_sun * sqrt(z / z_0)
    const double z = observer_distance_au * constants::AU;
    return std::sqrt(r_g_ * z);
}

// ============================================================================
// SGL Amplification
// ============================================================================

double GravitationalLens::peak_amplification(double wavelength_nm) const {
    // μ_max = 4π² r_g / λ  (Turyshev & Toth 2017, eq. 57)
    // At 550nm: μ ≈ 4 × 9.8696 × 2953.25 / 5.5e-7 ≈ 2.13 × 10^11
    const double lambda = wavelength_nm * 1e-9;
    return 4.0 * constants::PI * constants::PI * r_g_ / lambda;
}

double GravitationalLens::amplification_at_offset(double wavelength_nm, double offset_angle_rad) const {
    // μ(β) = μ_peak × J_0²(k r_g β)
    // where k = 2π/λ and β is the angular offset from optical axis
    //
    // More precisely: μ(β) = (4π² r_g / λ) × J_0²(2π r_g β / λ)
    const double lambda = wavelength_nm * 1e-9;
    const double arg = constants::TWO_PI * r_g_ * offset_angle_rad / lambda;
    const double j0 = bessel_j0(arg);
    return peak_amplification(wavelength_nm) * j0 * j0;
}

// ============================================================================
// Point Spread Function
// ============================================================================

GravitationalLens::PSFResult GravitationalLens::compute_psf(
    double wavelength_nm,
    double observer_distance_au,
    uint32_t grid_size,
    double field_of_view_mas
) const {
    PSFResult result;
    result.size = grid_size;

    const double lambda = wavelength_nm * 1e-9;
    const double z = observer_distance_au * constants::AU;
    const double fov_rad = field_of_view_mas * constants::MAS2RAD;
    result.pixel_scale_rad = fov_rad / grid_size;

    const double mu_peak = peak_amplification(wavelength_nm);
    result.peak_gain = mu_peak;

    result.data.resize(static_cast<size_t>(grid_size) * grid_size);

    const double center = grid_size / 2.0;
    double max_val = 0;
    double half_max_radius = 0;
    double total_energy = 0;

    // Generate the SGL PSF: proportional to J_0²(2π R_sun ρ / (λz))
    // where ρ is the physical offset at the observer plane
    // ρ = θ × z  (small angle)
    for (uint32_t iy = 0; iy < grid_size; ++iy) {
        for (uint32_t ix = 0; ix < grid_size; ++ix) {
            const double dx = (ix - center) * result.pixel_scale_rad;
            const double dy = (iy - center) * result.pixel_scale_rad;
            const double theta = std::sqrt(dx * dx + dy * dy);

            // Physical offset at observer plane
            const double rho = theta * z;

            // Argument to J0: 2π R_sun ρ / (λ z)
            // Equivalent to: 2π r_g θ / λ  (since r_g/z = θ_E², and R_sun²=r_g*z_0)
            const double arg = constants::TWO_PI * constants::R_SUN * rho / (lambda * z);
            const double j0_val = bessel_j0(arg);

            float psf_val = static_cast<float>(mu_peak * j0_val * j0_val);
            result.data[iy * grid_size + ix] = psf_val;

            if (psf_val > max_val) max_val = psf_val;
            total_energy += psf_val;
        }
    }

    // Normalize to peak = 1
    if (max_val > 0) {
        for (auto& v : result.data) v /= static_cast<float>(max_val);
    }

    // Find FWHM (first radius where PSF drops to 0.5)
    for (uint32_t ir = 0; ir < grid_size / 2; ++ir) {
        const double theta = ir * result.pixel_scale_rad;
        const double rho = theta * z;
        const double arg = constants::TWO_PI * constants::R_SUN * rho / (lambda * z);
        const double j0_val = bessel_j0(arg);
        if (j0_val * j0_val < 0.5) {
            half_max_radius = theta;
            break;
        }
    }
    result.fwhm_rad = 2.0 * half_max_radius;

    // Encircled energy at 50% — find radius containing half the total
    double cumulative = 0;
    result.encircled_energy_50 = fov_rad / 2.0;  // Default
    for (uint32_t ir = 1; ir < grid_size / 2; ++ir) {
        const double theta = ir * result.pixel_scale_rad;
        // Sum annular ring
        double ring_sum = 0;
        for (uint32_t iy = 0; iy < grid_size; ++iy) {
            for (uint32_t ix = 0; ix < grid_size; ++ix) {
                const double dx = (ix - center) * result.pixel_scale_rad;
                const double dy = (iy - center) * result.pixel_scale_rad;
                const double r = std::sqrt(dx * dx + dy * dy);
                if (r >= (ir - 1) * result.pixel_scale_rad && r < theta) {
                    ring_sum += result.data[iy * grid_size + ix];
                }
            }
        }
        cumulative += ring_sum;
        if (cumulative >= total_energy * 0.5) {
            result.encircled_energy_50 = theta;
            break;
        }
    }

    return result;
}

// ============================================================================
// Solar Corona Model
// ============================================================================

double GravitationalLens::corona_electron_density(double r_solar_radii) const {
    // Baumbach-Allen model for coronal electron density (cm^-3):
    // n_e(r) = 10^8 × (2.99 r^-16 + 1.55 r^-6 + 0.036 r^-1.5)
    // Valid for r > 1.0 solar radii
    if (r_solar_radii < 1.0) return 1e14;  // Inside the Sun — very high

    const double r = r_solar_radii;
    return 1e8 * (2.99 * std::pow(r, -16.0) +
                  1.55 * std::pow(r, -6.0) +
                  0.036 * std::pow(r, -1.5));
}

double GravitationalLens::corona_brightness(double r_solar_radii, double wavelength_nm) const {
    if (r_solar_radii < 1.0) return 1e20;  // Solar disk

    const double r = r_solar_radii;

    // K-corona (Thomson scattering of photospheric light by coronal electrons)
    // Surface brightness relative to solar disk center:
    // B_K = 10^(-6) × (3.67 r^-2.5 + 0.06 r^-4)  (Allen 1973)
    const double b_k = 1e-6 * (3.67 * std::pow(r, -2.5) + 0.06 * std::pow(r, -4.0));

    // F-corona (zodiacal dust scattering)
    // B_F = 10^(-6) × (1.2 r^-2.2)
    const double b_f = 1e-6 * 1.2 * std::pow(r, -2.2);

    // Total corona brightness relative to solar disk center
    const double b_total = b_k + b_f;

    // Convert to absolute photon flux (photons/m²/s/sr/nm)
    // Solar spectral irradiance at 550nm ≈ 1.9e18 photons/m²/s/nm
    const double solar_flux = 1.9e18 * std::pow(550.0 / wavelength_nm, 1.0);
    // Solar disk solid angle ≈ 6.8e-5 sr (at 1 AU)
    const double solar_solid_angle = constants::PI * std::pow(constants::R_SUN / constants::AU, 2.0);

    // Corona brightness per steradian per nm
    return b_total * solar_flux / solar_solid_angle;
}

// ============================================================================
// Signal-to-Noise Computation
// ============================================================================

GravitationalLens::SNRResult GravitationalLens::compute_snr(
    double planet_flux_photons_m2_s,
    double wavelength_nm,
    double observer_distance_au,
    double aperture_m,
    uint32_t num_satellites,
    double integration_time_s
) const {
    SNRResult result;
    result.integration_time_s = integration_time_s;

    // Collecting area of the swarm
    const double r = aperture_m / 2.0;
    const double area_per_sat = constants::PI * r * r;
    const double total_area = area_per_sat * num_satellites;

    // SGL gain at on-axis
    const double gain = peak_amplification(wavelength_nm);

    // Planet signal photons
    result.signal_photons = planet_flux_photons_m2_s * total_area * gain *
                           constants::CCD_QUANTUM_EFF * integration_time_s;

    // Corona background
    // At the SGL focal distance, the Einstein ring has angular radius θ_E
    // The corona at the Einstein ring (in solar radii from Sun center, as seen from observer)
    // θ_E corresponds to looking at the Sun's limb, so r ≈ 1 solar radius + small offset
    // More precisely: the Einstein ring grazes the solar limb at z = z_0
    // At z > z_0, the ring is slightly outside the limb
    const double r_ring_solar_radii = 1.0 + (observer_distance_au - z0_au_) / observer_distance_au;
    const double corona_flux = corona_brightness(r_ring_solar_radii, wavelength_nm);

    // Corona photons per satellite aperture per integration
    // Solid angle of one pixel ≈ (λ/D)²
    const double pixel_solid_angle = std::pow(wavelength_nm * 1e-9 / aperture_m, 2.0);
    result.corona_photons = corona_flux * pixel_solid_angle * total_area *
                           constants::CCD_QUANTUM_EFF * integration_time_s;

    // Detector noise
    const double n_pixels = static_cast<double>(num_satellites) * 100.0;  // ~100 pixels per source per sat
    result.dark_current_electrons = constants::CCD_DARK_CURRENT * n_pixels * integration_time_s;
    result.read_noise_electrons = constants::CCD_READ_NOISE_E * constants::CCD_READ_NOISE_E * n_pixels;

    // SNR (photon-counting regime)
    const double noise_variance = result.signal_photons + result.corona_photons +
                                  result.dark_current_electrons + result.read_noise_electrons;
    result.snr = result.signal_photons / std::sqrt(noise_variance + 1e-30);

    return result;
}

} // namespace solarlens
