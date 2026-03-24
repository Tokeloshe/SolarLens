#pragma once

#include "constants.hpp"
#include <vector>
#include <cmath>
#include <array>

namespace solarlens {

// ============================================================================
// SOLAR GRAVITATIONAL LENS (SGL) PHYSICS
//
// Based on Turyshev & Toth (2017): "Wave-optical treatment of the shadow
// region with the solar gravitational lens" and subsequent papers.
//
// The Sun acts as a convex lens with focal distance z_0 ≈ 547.6 AU.
// Beyond z_0, light from a distant source forms an Einstein ring around
// the Sun. The ring's gain is enormous: μ ~ 4π² r_g/λ ~ 10^11.
//
// Key physics: the PSF is NOT Gaussian. It's an Airy-like pattern
// modulated by the Bessel function J_0, arising from diffraction
// around the Sun's limb.
// ============================================================================
class GravitationalLens {
public:
    GravitationalLens();

    // ---------- Focal geometry ----------

    // Minimum focal distance for wavelength λ (geometric optics, achromatic)
    // z_0 = R_sun² / r_g ≈ 547.6 AU
    // Plasma adds a chromatic shift: z_min(λ) = z_0 / (1 - (ω_p/ω)²)
    // where ω_p is the plasma frequency at the solar limb
    double focal_distance_au(double wavelength_nm) const;

    // Einstein ring angular radius at observer distance z: θ_E = sqrt(r_g / z)
    double einstein_angle_rad(double observer_distance_au) const;

    // Physical Einstein ring radius at observer plane: b_E = sqrt(r_g * z)
    double einstein_ring_radius_m(double observer_distance_au) const;

    // ---------- SGL amplification ----------

    // Peak amplification for monochromatic point source on-axis:
    // μ = 4π² r_g / λ  (Turyshev & Toth 2017 eq. 57)
    double peak_amplification(double wavelength_nm) const;

    // Amplification as function of angular offset from optical axis:
    // μ(β) = μ_peak * |J_0(2π r_g β / λ)|²  (cylindrical symmetry)
    // β = angular separation of source from optical axis
    double amplification_at_offset(double wavelength_nm, double offset_angle_rad) const;

    // ---------- Point Spread Function ----------

    struct PSFResult {
        std::vector<float> data;     // 2D PSF kernel (row-major, size×size)
        uint32_t size;               // Pixel grid dimension
        double pixel_scale_rad;      // Angular size per pixel
        double fwhm_rad;             // Full-width half-maximum
        double peak_gain;            // Central amplification
        double encircled_energy_50;  // Radius containing 50% energy (rad)
    };

    // Compute the SGL PSF using the correct Bessel-function form
    // PSF(ρ) ∝ μ_peak * J_0²(2π R_sun ρ / (λ z))
    // where ρ is distance from optical axis at observer plane
    PSFResult compute_psf(
        double wavelength_nm,
        double observer_distance_au,
        uint32_t grid_size = 256,
        double field_of_view_mas = 100.0
    ) const;

    // ---------- Solar corona noise ----------

    // Electron density in solar corona (Baumbach-Allen model):
    // n_e(r) = 10^8 * (2.99 r^-16 + 1.55 r^-6 + 0.036 r^-1.5) cm^-3
    // r in solar radii
    double corona_electron_density(double r_solar_radii) const;

    // Corona surface brightness in photons/m²/s/sr/nm at wavelength λ
    // from Thomson scattering + F-corona dust scattering
    double corona_brightness(double r_solar_radii, double wavelength_nm) const;

    // Signal-to-noise ratio for an exoplanet observation
    // accounting for: SGL gain, corona noise, detector noise, photon statistics
    struct SNRResult {
        double signal_photons;        // Detected planet photons
        double corona_photons;        // Corona background photons
        double dark_current_electrons;
        double read_noise_electrons;
        double snr;                   // Final SNR
        double integration_time_s;    // Time used
    };

    SNRResult compute_snr(
        double planet_flux_photons_m2_s,  // Planet photon flux at Earth (ph/m²/s)
        double wavelength_nm,
        double observer_distance_au,
        double aperture_m,                 // Single satellite aperture
        uint32_t num_satellites,            // Number of collecting satellites
        double integration_time_s
    ) const;

private:
    // Bessel function J_0 (first kind, order 0)
    static double bessel_j0(double x);

    // Bessel function J_1 (first kind, order 1) — for Airy pattern
    static double bessel_j1(double x);

    // Pre-computed values
    double r_g_;            // Schwarzschild radius / 2 = gravitational radius
    double z0_;             // Minimum focal distance (m)
    double z0_au_;          // Minimum focal distance (AU)
};

} // namespace solarlens
