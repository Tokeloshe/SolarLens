#pragma once

#include "solarlens/physics/gravitational_lens.hpp"
#include "solarlens/physics/interferometry.hpp"
#include <vector>
#include <cstdint>
#include <array>

namespace solarlens {

// ============================================================================
// EXOPLANET DETECTION & CHARACTERIZATION PIPELINE
//
// Full imaging and spectroscopy pipeline for exoplanet detection using
// the SGL + distributed interferometer. Steps:
// 1. Collect photons across the satellite swarm
// 2. Cross-correlate pairs for visibility measurements
// 3. Reconstruct image using CLEAN or MEM deconvolution
// 4. Subtract stellar PSF and corona model
// 5. Detect point/extended sources above noise
// 6. Spectroscopic analysis for temperature, composition
// 7. Atmospheric biosignature assessment
// ============================================================================

class ExoplanetDetector {
public:
    static constexpr uint32_t SPECTRUM_BINS = 2048;
    static constexpr double SPECTRUM_MIN_NM = 300.0;   // UV
    static constexpr double SPECTRUM_MAX_NM = 2500.0;  // Near-IR

    // Detected planet characterization
    struct PlanetData {
        bool detected;
        float confidence;           // 0-1

        // Physical properties
        double radius_earth;        // Planet radius in Earth radii
        double mass_earth;          // Estimated mass (mass-radius relation)
        double temperature_k;       // Equilibrium temperature
        double albedo;              // Bond albedo
        double surface_gravity;     // m/s² (estimated)

        // Orbital properties
        double orbital_radius_au;   // Semi-major axis
        double orbital_period_days; // Kepler's third law
        bool in_habitable_zone;

        // Host star properties used
        double star_luminosity_solar;
        double star_distance_ly;

        // Atmospheric composition (volume mixing ratios)
        struct Atmosphere {
            double n2;              // Nitrogen
            double o2;              // Oxygen
            double co2;             // Carbon dioxide
            double h2o;             // Water vapor
            double ch4;             // Methane
            double o3;              // Ozone
            double nh3;             // Ammonia
            double so2;             // Sulfur dioxide

            double biosignature_score;  // 0-1 probability
            const char* biosignature_rationale;
        } atmosphere;

        // Resolved surface features (if SNR sufficient)
        struct SurfaceMap {
            std::vector<float> albedo_map;  // 2D albedo map
            uint32_t map_size;              // Pixels per side
            double pixel_scale_km;          // km per pixel on surface
            bool has_continents;
            bool has_oceans;
            double ocean_fraction;
            double cloud_fraction;
        } surface;
    };

    // ---------- Full Detection Pipeline ----------

    // Run complete detection pipeline on collected swarm data
    // interferometric_image: reconstructed image from CLEAN deconvolution
    // image_size: pixels per side
    // pixel_scale_rad: angular size per pixel
    // spectrum: integrated spectrum (SPECTRUM_BINS bins)
    // star_luminosity: host star luminosity in solar units
    // star_distance_ly: distance to target system
    // observer_distance_au: swarm distance from Sun
    PlanetData analyze(
        const std::vector<float>& interferometric_image,
        uint32_t image_size,
        double pixel_scale_rad,
        const std::vector<float>& spectrum,
        double star_luminosity,
        double star_distance_ly,
        double observer_distance_au
    );

    // ---------- Individual Pipeline Steps ----------

    // Richardson-Lucy deconvolution with proper PSF
    // Uses FFT-based convolution for efficiency
    static std::vector<float> richardson_lucy(
        const std::vector<float>& image,
        const std::vector<float>& psf,
        uint32_t size,
        uint32_t iterations,
        double regularization = 1e-10
    );

    // Detect point sources above threshold using matched filtering
    struct SourceDetection {
        double x_pixel, y_pixel;    // Sub-pixel position
        double flux;                // Integrated flux
        double snr;                 // Detection SNR
        double fwhm_pixels;         // Measured FWHM
        bool is_extended;           // Resolved source?
    };

    static std::vector<SourceDetection> detect_sources(
        const std::vector<float>& image,
        uint32_t size,
        double detection_threshold_sigma = 5.0
    );

    // Wien's law temperature from spectrum peak
    static double temperature_from_spectrum(
        const std::vector<float>& spectrum,
        double lambda_min_nm,
        double lambda_max_nm
    );

    // Analyze atmospheric absorption features
    static PlanetData::Atmosphere analyze_atmosphere(
        const std::vector<float>& spectrum,
        double lambda_min_nm,
        double lambda_max_nm
    );

    // Mass-radius relation (Chen & Kipping 2017 forecaster)
    static double estimate_mass(double radius_earth);

    // Habitable zone boundaries (Kopparapu et al. 2013)
    struct HabitableZone {
        double inner_au;    // Runaway greenhouse
        double outer_au;    // Maximum greenhouse
    };
    static HabitableZone habitable_zone(double star_luminosity_solar, double star_teff_k = 5772.0);

    // ---------- Synthetic Test Data ----------

    // Generate a synthetic Earth-like planet spectrum for testing
    static std::vector<float> generate_earth_spectrum(uint32_t bins, double distance_ly);

    // Generate a synthetic planet image (disk with features)
    static std::vector<float> generate_planet_image(
        uint32_t size,
        double planet_radius_pixels,
        double albedo,
        bool add_continents
    );
};

} // namespace solarlens
