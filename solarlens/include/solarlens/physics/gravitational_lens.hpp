#pragma once

#include "constants.hpp"
#include <cmath>
#include <array>

namespace solarlens {

// ============================================================================
// GRAVITATIONAL LENS CALCULATIONS
// ============================================================================
class GravitationalLensPhysics {
private:
    // Pre-computed values for efficiency
    const double schwarzschild_radius;
    const double einstein_radius_1au;

public:
    GravitationalLensPhysics();

    // Calculate exact focal distance for given wavelength
    // Based on: f = R_sun² / (4 * R_schwarzschild) * (1 + λ/λ_0)
    // Where λ_0 = 1 μm reference wavelength
    [[nodiscard]] double calculate_focal_distance_au(double wavelength_nm) const;

    // Calculate magnification factor
    // μ = (θ_E / θ_S)² where θ_E is Einstein ring angle, θ_S is source angle
    [[nodiscard]] double calculate_magnification(
        double source_distance_ly,
        double observer_distance_au,
        double impact_parameter_km
    ) const;

    // Calculate point spread function for gravitational lens
    struct PSF {
        static constexpr uint32_t SIZE = 256;
        std::array<std::array<float, SIZE>, SIZE> kernel;
        double fwhm_mas;  // Full width half maximum in milliarcseconds
    };

    [[nodiscard]] PSF calculate_psf(
        double wavelength_nm,
        double observer_distance_au
    ) const;

    // Solar corona noise model
    [[nodiscard]] double calculate_corona_brightness(
        double angular_distance_solar_radii,
        double wavelength_nm
    ) const;
};

} // namespace solarlens
