#include "solarlens/physics/gravitational_lens.hpp"
#include <cmath>

namespace solarlens {

GravitationalLensPhysics::GravitationalLensPhysics()
    : schwarzschild_radius(2.0 * constants::G * constants::M_SUN / (constants::C * constants::C))
    , einstein_radius_1au(std::sqrt(4.0 * constants::G * constants::M_SUN * constants::AU / (constants::C * constants::C)))
{}

double GravitationalLensPhysics::calculate_focal_distance_au(double wavelength_nm) const {
    const double wavelength_m = wavelength_nm * 1e-9;
    const double r_s = 2.0 * constants::G * constants::M_SUN / (constants::C * constants::C);

    // Base focal distance (achromatic)
    const double f_base = (constants::R_SUN * constants::R_SUN) / (4.0 * r_s);

    // Chromatic aberration correction (plasma dispersion)
    // Based on solar corona electron density model
    const double n_e = 1e8;  // electrons/cm³ at 1 AU
    const double plasma_freq = 8.98e3 * std::sqrt(n_e);  // Hz
    const double light_freq = constants::C / wavelength_m;
    const double dispersion_factor = 1.0 - (plasma_freq * plasma_freq) / (light_freq * light_freq);

    // Wavelength-dependent focal length
    const double f_chromatic = f_base * std::sqrt(dispersion_factor);

    // Convert to AU
    return f_chromatic / constants::AU;
}

double GravitationalLensPhysics::calculate_magnification(
    double source_distance_ly,
    double observer_distance_au,
    double impact_parameter_km
) const {
    const double d_s = source_distance_ly * constants::LY;  // Source distance
    const double d_l = observer_distance_au * constants::AU; // Lens distance

    // Einstein ring radius at observer position
    const double theta_E = std::sqrt(
        (4.0 * constants::G * constants::M_SUN / (constants::C * constants::C)) *
        (d_s - d_l) / (d_l * d_s)
    );

    // Physical Einstein ring radius
    const double r_E = theta_E * d_l;

    // Magnification from lens equation
    const double u = (impact_parameter_km * 1000.0) / r_E;  // Normalized impact

    if (u < 1e-6) {  // Perfect alignment
        return 1e12;  // Theoretical maximum
    }

    // Standard magnification formula
    const double mu = (u*u + 2.0) / (u * std::sqrt(u*u + 4.0));

    // Account for solar corona scattering (reduces magnification)
    const double corona_factor = std::exp(-0.1 / 500.0);

    return mu * corona_factor;
}

GravitationalLensPhysics::PSF GravitationalLensPhysics::calculate_psf(
    double wavelength_nm,
    double observer_distance_au
) const {
    PSF psf{};

    // Angular resolution (Rayleigh criterion modified for gravitational lens)
    const double lambda = wavelength_nm * 1e-9;
    const double baseline = observer_distance_au * constants::AU;
    const double theta_resolution = 1.22 * lambda / baseline;

    // Convert to milliarcseconds
    psf.fwhm_mas = theta_resolution * 206265000.0;

    // Generate 2D Gaussian PSF kernel
    const double sigma = PSF::SIZE / 6.0;  // Standard deviation in pixels
    const double center = PSF::SIZE / 2.0;

    for (uint32_t i = 0; i < PSF::SIZE; ++i) {
        for (uint32_t j = 0; j < PSF::SIZE; ++j) {
            const double r2 = (i - center) * (i - center) + (j - center) * (j - center);
            psf.kernel[i][j] = std::exp(-r2 / (2.0 * sigma * sigma));
        }
    }

    return psf;
}

double GravitationalLensPhysics::calculate_corona_brightness(
    double angular_distance_solar_radii,
    double wavelength_nm
) const {
    // F-corona (dust) + K-corona (electrons) model
    // Based on Allen's Astrophysical Quantities

    const double r = angular_distance_solar_radii;

    if (r < 1.0) {
        return 1e10;  // Solar disk - saturated
    }

    // K-corona (Thomson scattering)
    const double k_corona = 1e6 * std::pow(r, -2.5);

    // F-corona (zodiacal light)
    const double f_corona = 1e5 * std::pow(r, -2.2);

    // Wavelength dependence
    const double lambda_factor = std::pow(wavelength_nm / 550.0, -1.2);

    return (k_corona + f_corona) * lambda_factor;
}

} // namespace solarlens
