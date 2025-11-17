#include "solarlens/imaging/exoplanet_detector.hpp"
#include <cmath>

namespace solarlens {

ExoplanetDetector::PlanetData ExoplanetDetector::detect_exoplanet(
    const uint16_t* sensor_data,
    uint32_t integration_time_seconds,
    double target_distance_ly,
    double wavelength_nm
) {
    PlanetData planet{};

    // Step 1: Accumulate photons with Poisson noise statistics
    accumulate_photons(sensor_data, integration_time_seconds);

    // Step 2: Subtract solar corona using model
    subtract_corona_model(target_distance_ly);

    // Step 3: Deconvolve PSF using Richardson-Lucy algorithm
    const auto psf = physics.calculate_psf(wavelength_nm, constants::FOCAL_OPTIMAL_AU);
    richardson_lucy_deconvolution(psf, 50);  // 50 iterations

    // Step 4: Detect point sources above noise threshold
    const auto detection = detect_point_source();

    if (detection.found) {
        planet.detected = true;
        planet.confidence = detection.snr / 10.0f;  // SNR>10 = high confidence

        // Step 5: Estimate physical parameters
        planet.radius_earth = estimate_radius_from_flux(detection.flux);
        planet.temperature_kelvin = estimate_temperature(detection.spectrum);
        planet.albedo = estimate_albedo(detection.flux, planet.temperature_kelvin);

        // Step 6: Determine orbital parameters from spectral shift
        planet.orbital_radius_au = estimate_orbit_from_doppler(detection.doppler_shift);

        // Step 7: Check habitable zone (liquid water)
        const double hz_inner = 0.95 * std::sqrt(target_luminosity);
        const double hz_outer = 1.37 * std::sqrt(target_luminosity);
        planet.in_habitable_zone = (planet.orbital_radius_au > hz_inner &&
                                   planet.orbital_radius_au < hz_outer);

        // Step 8: Atmospheric analysis via spectroscopy
        planet.atmosphere = analyze_atmosphere(detection.spectrum);
    }

    return planet;
}

void ExoplanetDetector::accumulate_photons(const uint16_t* data, uint32_t integration_seconds) {
    // Photon accumulation with shot noise
    const double dark_current = 0.01;        // e-/pixel/s at -80°C

    for (uint32_t y = 0; y < IMAGE_SIZE; ++y) {
        for (uint32_t x = 0; x < IMAGE_SIZE; ++x) {
            const uint32_t idx = y * IMAGE_SIZE + x;
            const double signal = data[idx] * integration_seconds;
            const double noise = std::sqrt(signal + dark_current * integration_seconds);
            raw_image[y][x] = signal + noise;  // Poisson noise
        }
    }
}

void ExoplanetDetector::subtract_corona_model(double distance_ly) {
    // Model-based corona subtraction
    (void)distance_ly; // Reserved for future distance-dependent corona model
    const double center_x = IMAGE_SIZE / 2.0;
    const double center_y = IMAGE_SIZE / 2.0;

    for (uint32_t y = 0; y < IMAGE_SIZE; ++y) {
        for (uint32_t x = 0; x < IMAGE_SIZE; ++x) {
            const double r = std::sqrt((x - center_x) * (x - center_x) +
                                      (y - center_y) * (y - center_y));
            const double r_solar_radii = r / 100.0;  // Pixel scale
            const double corona = physics.calculate_corona_brightness(r_solar_radii, 550.0);
            processed_image[y][x] = raw_image[y][x] - corona;
        }
    }
}

void ExoplanetDetector::richardson_lucy_deconvolution(const GravitationalLensPhysics::PSF& psf, uint32_t iterations) {
    // Richardson-Lucy deconvolution for PSF removal
    (void)psf; // Simplified version uses fixed kernel; full version would use PSF
    std::array<std::array<float, IMAGE_SIZE>, IMAGE_SIZE> estimate = processed_image;
    std::array<std::array<float, IMAGE_SIZE>, IMAGE_SIZE> ratio{};

    for (uint32_t iter = 0; iter < iterations; ++iter) {
        // Forward convolution
        for (uint32_t y = 0; y < IMAGE_SIZE; ++y) {
            for (uint32_t x = 0; x < IMAGE_SIZE; ++x) {
                float sum = 0;
                // Simplified convolution (full would use FFT)
                for (int dy = -2; dy <= 2; ++dy) {
                    for (int dx = -2; dx <= 2; ++dx) {
                        int ny = static_cast<int>(y) + dy;
                        int nx = static_cast<int>(x) + dx;
                        if (ny >= 0 && ny < static_cast<int>(IMAGE_SIZE) &&
                            nx >= 0 && nx < static_cast<int>(IMAGE_SIZE)) {
                            sum += estimate[ny][nx] * 0.04f;  // 5x5 kernel
                        }
                    }
                }
                ratio[y][x] = processed_image[y][x] / (sum + 1e-10f);
            }
        }

        // Backward convolution and update
        for (uint32_t y = 0; y < IMAGE_SIZE; ++y) {
            for (uint32_t x = 0; x < IMAGE_SIZE; ++x) {
                float sum = 0;
                for (int dy = -2; dy <= 2; ++dy) {
                    for (int dx = -2; dx <= 2; ++dx) {
                        int ny = static_cast<int>(y) + dy;
                        int nx = static_cast<int>(x) + dx;
                        if (ny >= 0 && ny < static_cast<int>(IMAGE_SIZE) &&
                            nx >= 0 && nx < static_cast<int>(IMAGE_SIZE)) {
                            sum += ratio[ny][nx] * 0.04f;
                        }
                    }
                }
                estimate[y][x] *= sum;
            }
        }
    }

    processed_image = estimate;
}

ExoplanetDetector::Detection ExoplanetDetector::detect_point_source() {
    Detection det{};

    // Find brightest pixel (simplified - real would use matched filter)
    float max_pixel = 0;
    uint32_t max_x = 0, max_y = 0;

    for (uint32_t y = 10; y < IMAGE_SIZE - 10; ++y) {
        for (uint32_t x = 10; x < IMAGE_SIZE - 10; ++x) {
            if (processed_image[y][x] > max_pixel) {
                max_pixel = processed_image[y][x];
                max_x = x;
                max_y = y;
            }
        }
    }

    // Calculate SNR in 5x5 aperture
    float signal = 0;
    float noise = 0;
    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            signal += processed_image[max_y + dy][max_x + dx];
        }
    }

    // Estimate noise from annulus
    for (int dy = -10; dy <= 10; ++dy) {
        for (int dx = -10; dx <= 10; ++dx) {
            if (std::abs(dy) > 5 || std::abs(dx) > 5) {
                noise += processed_image[max_y + dy][max_x + dx] * processed_image[max_y + dy][max_x + dx];
            }
        }
    }
    noise = std::sqrt(noise / 300.0);  // RMS noise

    det.snr = signal / (noise + 1e-10);
    det.found = det.snr > 5.0;  // 5-sigma detection threshold
    det.flux = signal;

    return det;
}

double ExoplanetDetector::estimate_radius_from_flux(double flux) {
    // Estimate radius from photon flux
    // Assumes Earth-like albedo of 0.3
    const double albedo = 0.3;
    const double stellar_luminosity = constants::L_SUN;  // Assume sun-like
    const double distance = 10.0 * constants::LY;  // 10 light years

    // Planet flux = (R_p/d)² * albedo * L_star / (4π * a²)
    // Solving for R_p
    const double radius_m = std::sqrt(flux * 4.0 * M_PI * distance * distance /
                                     (albedo * stellar_luminosity));
    const double radius_earth = radius_m / 6.371e6;

    return radius_earth;
}

double ExoplanetDetector::estimate_temperature(const std::array<float, SPECTRUM_BINS>& spectrum) {
    // Estimate temperature from blackbody peak (Wien's law)
    uint32_t peak_bin = 0;
    float max_intensity = 0;

    for (uint32_t i = 0; i < SPECTRUM_BINS; ++i) {
        if (spectrum[i] > max_intensity) {
            max_intensity = spectrum[i];
            peak_bin = i;
        }
    }

    // Convert bin to wavelength (400-2400 nm range)
    const double wavelength_nm = 400.0 + (peak_bin * 2000.0 / SPECTRUM_BINS);
    const double wavelength_m = wavelength_nm * 1e-9;

    // Wien's displacement law: λ_max * T = 2.897e-3 m·K
    const double temperature = 2.897e-3 / wavelength_m;

    return temperature;
}

double ExoplanetDetector::estimate_albedo(double flux, double temperature) {
    // Estimate Bond albedo from flux and temperature
    (void)flux; // Simplified version uses temperature; full version would use flux
    const double stefan_boltzmann = 5.67e-8;  // W/m²/K⁴
    const double emitted_power = stefan_boltzmann * std::pow(temperature, 4);
    const double incident_power = constants::L_SUN / (4.0 * M_PI * constants::AU * constants::AU);

    return 1.0 - (emitted_power / incident_power);
}

double ExoplanetDetector::estimate_orbit_from_doppler(double doppler_shift) {
    // Estimate orbital radius from Doppler shift
    // Assumes circular orbit, edge-on
    const double orbital_velocity = doppler_shift * constants::C;

    // Kepler's third law: v = sqrt(GM/r)
    const double orbital_radius = constants::G * constants::M_SUN /
                                 (orbital_velocity * orbital_velocity);

    return orbital_radius / constants::AU;
}

typename ExoplanetDetector::PlanetData::Atmosphere ExoplanetDetector::analyze_atmosphere(const std::array<float, SPECTRUM_BINS>& spectrum) {
    typename PlanetData::Atmosphere atm{};

    // Detect molecular absorption lines
    // Wavelengths in nm for key molecules
    constexpr struct {
        uint32_t wavelength_nm;
        float PlanetData::Atmosphere::* target;
        const char* molecule;
    } absorption_lines[] = {
        {760, &PlanetData::Atmosphere::oxygen, "O2"},      // O2 A-band
        {1640, &PlanetData::Atmosphere::methane, "CH4"},   // Methane
        {940, &PlanetData::Atmosphere::water, "H2O"},      // Water vapor
        {2013, &PlanetData::Atmosphere::co2, "CO2"},       // Carbon dioxide
        {2300, &PlanetData::Atmosphere::nitrogen, "N2"}    // Nitrogen
    };

    for (const auto& line : absorption_lines) {
        const uint32_t bin = (line.wavelength_nm - 400) * SPECTRUM_BINS / 2000;
        if (bin < SPECTRUM_BINS && bin >= 10 && bin + 10 < SPECTRUM_BINS) {
            // Check for absorption dip
            const float continuum = (spectrum[bin-10] + spectrum[bin+10]) / 2.0f;
            const float depth = (continuum - spectrum[bin]) / continuum;
            atm.*line.target = depth * 100.0f;  // Convert to percentage
        }
    }

    // Calculate biosignature score
    // Based on oxygen-methane disequilibrium
    const bool oxygen_present = atm.oxygen > 1.0f;
    const bool methane_present = atm.methane > 0.01f;
    const bool water_present = atm.water > 0.1f;

    if (oxygen_present && methane_present) {
        atm.biosignature_score = 0.9f;  // Strong biosignature
    } else if (oxygen_present && water_present) {
        atm.biosignature_score = 0.6f;  // Moderate biosignature
    } else if (water_present) {
        atm.biosignature_score = 0.3f;  // Weak biosignature
    } else {
        atm.biosignature_score = 0.0f;
    }

    return atm;
}

} // namespace solarlens
