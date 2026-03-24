#pragma once

#include "constants.hpp"
#include "../physics/gravitational_lens.hpp"
#include <vector>
#include <complex>
#include <cstdint>

namespace solarlens {

// ============================================================================
// APERTURE SYNTHESIS INTERFEROMETRY
//
// Each satellite in the Starlink-scale swarm acts as an element of a
// sparse-aperture interferometer. Pairs of satellites form baselines
// that sample the UV plane. The Van Cittert-Zernike theorem relates
// the visibility function (Fourier transform of sky brightness) to
// the cross-correlations between satellite pairs.
//
// Combined with the SGL's enormous gain, this enables imaging of
// exoplanet surfaces at ~10 km resolution from hundreds of light-years.
// ============================================================================

struct Baseline {
    uint32_t sat_i;           // First satellite index
    uint32_t sat_j;           // Second satellite index
    double u;                 // UV coordinate u (wavelengths)
    double v;                 // UV coordinate v (wavelengths)
    double length_m;          // Physical baseline length (m)
    std::complex<double> visibility;  // Measured complex visibility
};

struct UVPoint {
    double u, v;              // In units of wavelengths
    double weight;            // Sampling weight (natural/uniform/Briggs)
};

class Interferometer {
public:
    // ---------- UV Coverage ----------

    // Compute all baselines from satellite positions (in local plane perpendicular to target)
    // positions: Nx2 array of (x,y) offsets from swarm center in meters
    // wavelength_m: observing wavelength for converting to UV coordinates
    // Returns: vector of baselines with UV coordinates
    static std::vector<Baseline> compute_baselines(
        const std::vector<std::array<double, 2>>& positions,
        double wavelength_m
    );

    // Compute UV coverage metric (fraction of UV plane filled)
    // Higher coverage = better image fidelity
    // uv_max: maximum UV radius to consider
    // grid_size: resolution of UV grid for coverage calculation
    static double uv_coverage_fraction(
        const std::vector<Baseline>& baselines,
        double uv_max,
        uint32_t grid_size = 256
    );

    // ---------- Visibility Measurement ----------

    // Simulate visibility measurement for a model image
    // image: source brightness distribution (grid_size × grid_size, row-major)
    // grid_size: image dimension
    // pixel_scale_rad: angular size per pixel
    // baselines: UV sampling points (u,v coordinates filled in)
    // Output: fills in visibility field of each baseline
    static void measure_visibilities(
        const std::vector<float>& image,
        uint32_t grid_size,
        double pixel_scale_rad,
        std::vector<Baseline>& baselines
    );

    // ---------- Image Reconstruction ----------

    // Dirty image from visibility data (inverse Fourier transform with sampling)
    // output_size: reconstructed image dimension
    // pixel_scale_rad: desired angular pixel scale
    static std::vector<float> make_dirty_image(
        const std::vector<Baseline>& baselines,
        uint32_t output_size,
        double pixel_scale_rad
    );

    // Dirty beam (PSF of the interferometer) from UV coverage
    static std::vector<float> make_dirty_beam(
        const std::vector<Baseline>& baselines,
        uint32_t output_size,
        double pixel_scale_rad
    );

    // CLEAN deconvolution (Högbom 1974)
    // Iteratively subtracts the dirty beam from brightest point in residual
    // dirty_image: input dirty image
    // dirty_beam: interferometer PSF
    // output_size: image dimension
    // gain: loop gain (typically 0.1)
    // threshold: stopping flux threshold
    // max_iterations: iteration limit
    static std::vector<float> clean_deconvolution(
        const std::vector<float>& dirty_image,
        const std::vector<float>& dirty_beam,
        uint32_t output_size,
        double gain = 0.1,
        double threshold = 1e-4,
        uint32_t max_iterations = 10000
    );

    // ---------- Resolution & Sensitivity ----------

    // Angular resolution from maximum baseline: θ = λ / B_max
    static double angular_resolution_rad(double wavelength_m, double max_baseline_m) {
        return wavelength_m / max_baseline_m;
    }

    // Linear resolution on target at given distance
    static double linear_resolution_m(double wavelength_m, double max_baseline_m, double distance_m) {
        return distance_m * wavelength_m / max_baseline_m;
    }

    // Combined collecting area of swarm
    static double total_collecting_area(uint32_t num_satellites, double aperture_m) {
        double r = aperture_m / 2.0;
        return num_satellites * constants::PI * r * r;
    }

    // Sensitivity: minimum detectable flux density
    // SEFD: system equivalent flux density per element
    // N: number of elements
    // t_int: integration time (s)
    // bw: bandwidth (Hz)
    static double point_source_sensitivity(double sefd, uint32_t N, double t_int, double bw) {
        double n_baselines = static_cast<double>(N) * (N - 1) / 2.0;
        return sefd / (std::sqrt(2.0 * bw * t_int * n_baselines));
    }

    // ---------- SGL-Enhanced Interferometry ----------

    // With the SGL, the exoplanet's light is amplified and projected onto
    // the Einstein ring. The swarm samples this ring pattern.
    // This function computes the effective resolution combining
    // SGL gain with interferometric baseline.
    struct SGLInterferometryResult {
        double effective_aperture_m;     // SGL-equivalent aperture
        double angular_resolution_rad;   // Combined resolution
        double linear_resolution_km;     // On-target resolution
        double total_gain;               // SGL gain × interferometric gain
        double required_integration_s;   // For SNR=10 detection
        double image_pixels;             // Resolvable pixels across planet
    };

    static SGLInterferometryResult sgl_interferometry_performance(
        double wavelength_nm,
        double observer_distance_au,
        double target_distance_ly,
        double target_radius_earth,
        uint32_t num_satellites,
        double max_baseline_m,
        double aperture_per_sat_m
    );
};

} // namespace solarlens
