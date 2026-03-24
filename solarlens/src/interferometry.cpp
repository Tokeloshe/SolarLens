#include "solarlens/physics/interferometry.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace solarlens {

// ============================================================================
// UV Coverage Computation
// ============================================================================

std::vector<Baseline> Interferometer::compute_baselines(
    const std::vector<std::array<double, 2>>& positions,
    double wavelength_m
) {
    const uint32_t N = static_cast<uint32_t>(positions.size());
    std::vector<Baseline> baselines;
    baselines.reserve(static_cast<size_t>(N) * (N - 1) / 2);

    for (uint32_t i = 0; i < N; ++i) {
        for (uint32_t j = i + 1; j < N; ++j) {
            Baseline b;
            b.sat_i = i;
            b.sat_j = j;

            double dx = positions[j][0] - positions[i][0];
            double dy = positions[j][1] - positions[i][1];
            b.length_m = std::sqrt(dx * dx + dy * dy);

            // UV coordinates in units of wavelengths
            b.u = dx / wavelength_m;
            b.v = dy / wavelength_m;
            b.visibility = std::complex<double>(0.0, 0.0);

            baselines.push_back(b);

            // Also add conjugate baseline (-u, -v)
            Baseline b_conj = b;
            b_conj.u = -b.u;
            b_conj.v = -b.v;
            baselines.push_back(b_conj);
        }
    }

    return baselines;
}

double Interferometer::uv_coverage_fraction(
    const std::vector<Baseline>& baselines,
    double uv_max,
    uint32_t grid_size
) {
    // Grid the UV plane and count filled cells
    std::vector<uint8_t> grid(static_cast<size_t>(grid_size) * grid_size, 0);
    uint32_t filled = 0;

    const double cell_size = 2.0 * uv_max / grid_size;

    for (const auto& b : baselines) {
        double u_norm = b.u / uv_max;
        double v_norm = b.v / uv_max;

        if (std::abs(u_norm) > 1.0 || std::abs(v_norm) > 1.0) continue;

        int ix = static_cast<int>((u_norm + 1.0) * 0.5 * grid_size);
        int iy = static_cast<int>((v_norm + 1.0) * 0.5 * grid_size);

        ix = std::clamp(ix, 0, static_cast<int>(grid_size - 1));
        iy = std::clamp(iy, 0, static_cast<int>(grid_size - 1));

        size_t idx = static_cast<size_t>(iy) * grid_size + ix;
        if (!grid[idx]) {
            grid[idx] = 1;
            filled++;
        }
    }

    // UV plane is a disk, so max possible filled cells ≈ π/4 × grid_size²
    double max_cells = constants::PI / 4.0 * grid_size * grid_size;
    (void)cell_size;
    return static_cast<double>(filled) / max_cells;
}

// ============================================================================
// Visibility Measurement
// ============================================================================

void Interferometer::measure_visibilities(
    const std::vector<float>& image,
    uint32_t grid_size,
    double pixel_scale_rad,
    std::vector<Baseline>& baselines
) {
    // Direct Fourier Transform at each baseline's (u,v) point
    // V(u,v) = Σ I(l,m) exp(-2πi(ul + vm))
    // where l,m are direction cosines (≈ angular coordinates for small fields)

    const double center = grid_size / 2.0;

    for (auto& b : baselines) {
        std::complex<double> vis(0.0, 0.0);

        for (uint32_t iy = 0; iy < grid_size; ++iy) {
            double m = (iy - center) * pixel_scale_rad;  // Direction cosine
            for (uint32_t ix = 0; ix < grid_size; ++ix) {
                double l = (ix - center) * pixel_scale_rad;
                double I = image[iy * grid_size + ix];

                if (I == 0.0f) continue;

                // Phase: -2π(ul + vm)
                double phase = -constants::TWO_PI * (b.u * l + b.v * m);
                vis += std::complex<double>(I * std::cos(phase), I * std::sin(phase));
            }
        }

        b.visibility = vis;
    }
}

// ============================================================================
// Dirty Image (Inverse DFT with sampling)
// ============================================================================

std::vector<float> Interferometer::make_dirty_image(
    const std::vector<Baseline>& baselines,
    uint32_t output_size,
    double pixel_scale_rad
) {
    std::vector<float> image(static_cast<size_t>(output_size) * output_size, 0.0f);
    const double center = output_size / 2.0;

    for (uint32_t iy = 0; iy < output_size; ++iy) {
        double m = (iy - center) * pixel_scale_rad;
        for (uint32_t ix = 0; ix < output_size; ++ix) {
            double l = (ix - center) * pixel_scale_rad;

            std::complex<double> sum(0.0, 0.0);

            for (const auto& b : baselines) {
                double phase = constants::TWO_PI * (b.u * l + b.v * m);
                sum += b.visibility * std::complex<double>(std::cos(phase), std::sin(phase));
            }

            image[iy * output_size + ix] = static_cast<float>(sum.real() / baselines.size());
        }
    }

    return image;
}

std::vector<float> Interferometer::make_dirty_beam(
    const std::vector<Baseline>& baselines,
    uint32_t output_size,
    double pixel_scale_rad
) {
    // Dirty beam = inverse FT of the sampling function (all visibilities = 1)
    std::vector<float> beam(static_cast<size_t>(output_size) * output_size, 0.0f);
    const double center = output_size / 2.0;

    for (uint32_t iy = 0; iy < output_size; ++iy) {
        double m = (iy - center) * pixel_scale_rad;
        for (uint32_t ix = 0; ix < output_size; ++ix) {
            double l = (ix - center) * pixel_scale_rad;

            double sum = 0.0;
            for (const auto& b : baselines) {
                double phase = constants::TWO_PI * (b.u * l + b.v * m);
                sum += std::cos(phase);
            }

            beam[iy * output_size + ix] = static_cast<float>(sum / baselines.size());
        }
    }

    return beam;
}

// ============================================================================
// CLEAN Deconvolution (Högbom 1974)
// ============================================================================

std::vector<float> Interferometer::clean_deconvolution(
    const std::vector<float>& dirty_image,
    const std::vector<float>& dirty_beam,
    uint32_t output_size,
    double gain,
    double threshold,
    uint32_t max_iterations
) {
    const size_t npix = static_cast<size_t>(output_size) * output_size;

    // Residual starts as the dirty image
    std::vector<float> residual = dirty_image;
    // Clean components (delta functions)
    std::vector<float> components(npix, 0.0f);
    // Final clean image
    std::vector<float> clean_image(npix, 0.0f);

    const int half = output_size / 2;

    for (uint32_t iter = 0; iter < max_iterations; ++iter) {
        // Find peak in residual
        float max_val = 0;
        uint32_t max_idx = 0;
        for (uint32_t i = 0; i < npix; ++i) {
            if (std::abs(residual[i]) > std::abs(max_val)) {
                max_val = residual[i];
                max_idx = i;
            }
        }

        if (std::abs(max_val) < threshold) break;

        // Record clean component
        float component = static_cast<float>(gain) * max_val;
        components[max_idx] += component;

        // Subtract shifted and scaled dirty beam from residual
        int peak_y = max_idx / output_size;
        int peak_x = max_idx % output_size;

        for (uint32_t iy = 0; iy < output_size; ++iy) {
            int beam_y = static_cast<int>(iy) - peak_y + half;
            if (beam_y < 0 || beam_y >= static_cast<int>(output_size)) continue;

            for (uint32_t ix = 0; ix < output_size; ++ix) {
                int beam_x = static_cast<int>(ix) - peak_x + half;
                if (beam_x < 0 || beam_x >= static_cast<int>(output_size)) continue;

                residual[iy * output_size + ix] -=
                    component * dirty_beam[beam_y * output_size + beam_x];
            }
        }
    }

    // Restore: convolve clean components with a Gaussian "clean beam"
    // Use a simple Gaussian approximation of the main lobe
    double beam_sigma = 2.0;  // pixels
    for (uint32_t iy = 0; iy < output_size; ++iy) {
        for (uint32_t ix = 0; ix < output_size; ++ix) {
            double sum = 0;
            // Convolve with clean beam
            for (int dy = -6; dy <= 6; ++dy) {
                int sy = static_cast<int>(iy) + dy;
                if (sy < 0 || sy >= static_cast<int>(output_size)) continue;
                for (int dx = -6; dx <= 6; ++dx) {
                    int sx = static_cast<int>(ix) + dx;
                    if (sx < 0 || sx >= static_cast<int>(output_size)) continue;

                    double g = std::exp(-(dx * dx + dy * dy) / (2.0 * beam_sigma * beam_sigma));
                    sum += components[sy * output_size + sx] * g;
                }
            }
            // Add residual
            clean_image[iy * output_size + ix] = static_cast<float>(sum) + residual[iy * output_size + ix];
        }
    }

    return clean_image;
}

// ============================================================================
// SGL-Enhanced Interferometry Performance
// ============================================================================

Interferometer::SGLInterferometryResult Interferometer::sgl_interferometry_performance(
    double wavelength_nm,
    double observer_distance_au,
    double target_distance_ly,
    double target_radius_earth,
    uint32_t num_satellites,
    double max_baseline_m,
    double aperture_per_sat_m
) {
    SGLInterferometryResult result{};

    const double lambda = wavelength_nm * 1e-9;
    const double target_distance_m = target_distance_ly * constants::LY;
    const double planet_radius_m = target_radius_earth * 6.371e6;

    // SGL gain: μ = 4π² r_g / λ
    const double sgl_gain = 4.0 * constants::PI * constants::PI *
                           constants::R_SCHWARZSCHILD / lambda;

    // Effective aperture from SGL:
    // The SGL acts as a lens with effective aperture ~ sqrt(μ) × D_telescope
    // More precisely: the SGL collects light from an annulus of width ~λ
    // around the Einstein ring, so effective collecting area ≈ 2π b_E λ
    // where b_E = sqrt(r_g × z)
    double b_E = std::sqrt(constants::R_SCHWARZSCHILD * observer_distance_au * constants::AU);
    result.effective_aperture_m = std::sqrt(2.0 * constants::PI * b_E * lambda);

    // Angular resolution from the interferometric baseline
    // (The SGL amplifies the signal but the resolution is set by the baseline)
    result.angular_resolution_rad = lambda / max_baseline_m;

    // Linear resolution on the target
    result.linear_resolution_km = target_distance_m * result.angular_resolution_rad / 1000.0;

    // Total gain: SGL gain × interferometric collecting area gain
    double collecting_area = total_collecting_area(num_satellites, aperture_per_sat_m);
    result.total_gain = sgl_gain * collecting_area / (lambda * lambda);

    // How many resolvable pixels across the planet?
    double planet_angular_size = 2.0 * planet_radius_m / target_distance_m;
    result.image_pixels = planet_angular_size / result.angular_resolution_rad;

    // Required integration time for SNR=10 on each resolution element
    // Assuming Earth-like planet at target distance:
    // Planet flux ≈ albedo × L_sun / (4π d²) × (R_p / d_target)²
    // Then amplified by SGL gain
    double planet_flux = 0.3 * constants::L_SUN /
                        (4.0 * constants::PI * target_distance_m * target_distance_m) *
                        std::pow(planet_radius_m / (constants::AU), 2.0);
    // Convert to photons: E_photon = hc/λ
    double photon_energy = constants::H * constants::C / lambda;
    double photon_flux = planet_flux / photon_energy;  // photons/m²/s at Earth

    // Detected photons per second per resolution element
    double det_rate = photon_flux * collecting_area * sgl_gain * constants::CCD_QUANTUM_EFF;
    // Target SNR=10 → need 100/det_rate seconds (photon-noise dominated)
    if (det_rate > 0) {
        result.required_integration_s = 100.0 / det_rate;
    } else {
        result.required_integration_s = 1e20;
    }

    return result;
}

} // namespace solarlens
