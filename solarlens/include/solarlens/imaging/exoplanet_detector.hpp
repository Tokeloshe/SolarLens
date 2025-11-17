#pragma once

#include "solarlens/physics/gravitational_lens.hpp"
#include <array>
#include <cstdint>

namespace solarlens {

// ============================================================================
// EXOPLANET DETECTION & IMAGING
// ============================================================================
class ExoplanetDetector {
private:
    static constexpr uint32_t IMAGE_SIZE = 1024;
    static constexpr uint32_t SPECTRUM_BINS = 2048;

    // Pre-allocated buffers (no dynamic allocation)
    std::array<std::array<float, IMAGE_SIZE>, IMAGE_SIZE> raw_image;
    std::array<std::array<float, IMAGE_SIZE>, IMAGE_SIZE> processed_image;
    std::array<float, SPECTRUM_BINS> spectrum;

    GravitationalLensPhysics physics;

public:
    struct PlanetData {
        bool detected;
        double radius_earth;
        double orbital_radius_au;
        double temperature_kelvin;
        double albedo;
        bool in_habitable_zone;
        float confidence;

        // Atmospheric composition (volume fractions)
        struct Atmosphere {
            float oxygen;      // O2
            float methane;     // CH4
            float water;       // H2O
            float co2;         // CO2
            float nitrogen;    // N2
            float biosignature_score;  // 0-1 probability of life
        } atmosphere;
    };

    // Main detection algorithm
    [[nodiscard]] PlanetData detect_exoplanet(
        const uint16_t* sensor_data,
        uint32_t integration_time_seconds,
        double target_distance_ly,
        double wavelength_nm
    );

private:
    struct Detection {
        bool found;
        double flux;
        double snr;
        double doppler_shift;
        std::array<float, SPECTRUM_BINS> spectrum;
    };

    void accumulate_photons(const uint16_t* data, uint32_t integration_seconds);
    void subtract_corona_model(double distance_ly);
    void richardson_lucy_deconvolution(const GravitationalLensPhysics::PSF& psf, uint32_t iterations);
    Detection detect_point_source();
    double estimate_radius_from_flux(double flux);
    double estimate_temperature(const std::array<float, SPECTRUM_BINS>& spectrum);
    double estimate_albedo(double flux, double temperature);
    double estimate_orbit_from_doppler(double doppler_shift);
    typename PlanetData::Atmosphere analyze_atmosphere(const std::array<float, SPECTRUM_BINS>& spectrum);

    double target_luminosity = 1.0;  // Solar units
};

} // namespace solarlens
