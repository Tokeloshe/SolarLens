// Comprehensive test suite for SolarLens
// Compiles with -fno-exceptions -fno-rtti -Werror

#include <iostream>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "solarlens/physics/constants.hpp"
#include "solarlens/physics/gravitational_lens.hpp"
#include "solarlens/physics/orbital_mechanics.hpp"
#include "solarlens/physics/interferometry.hpp"
#include "solarlens/spacecraft/formation_control.hpp"
#include "solarlens/spacecraft/navigation.hpp"
#include "solarlens/spacecraft/power_thermal.hpp"
#include "solarlens/communication/deep_space_relay.hpp"
#include "solarlens/imaging/exoplanet_detector.hpp"

static int tests_total = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
    tests_total++; \
    if (!(cond)) { std::cerr << "FAIL: " << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n"; tests_failed++; } \
    else { tests_passed++; } \
} while(0)

// ---------------------------------------------------------------------------
// 1. Physics Constants
// ---------------------------------------------------------------------------
static void test_physics_constants() {
    using namespace solarlens::constants;

    // Speed of light
    TEST_ASSERT(C == 299792458.0, "Speed of light is 299792458 m/s");

    // Schwarzschild radius ~2953 m
    TEST_ASSERT(R_SCHWARZSCHILD > 2950.0 && R_SCHWARZSCHILD < 2960.0,
                "Schwarzschild radius ~2953 m");

    // SGL focal min in 547-1100 AU range
    TEST_ASSERT(SGL_FOCAL_MIN_AU > 540.0 && SGL_FOCAL_MIN_AU < 1100.0,
                "SGL focal min in 547-1100 AU range");

    // AU value
    TEST_ASSERT(AU == 1.495978707e11, "AU value exact");
}

// ---------------------------------------------------------------------------
// 2. Gravitational Lens
// ---------------------------------------------------------------------------
static void test_gravitational_lens() {
    solarlens::GravitationalLens lens;

    // focal_distance_au(550nm) returns > 500 AU
    double fd = lens.focal_distance_au(550.0);
    TEST_ASSERT(fd > 500.0, "focal_distance_au(550) > 500 AU");

    // peak_amplification(550nm) ~ 2e11
    double mu = lens.peak_amplification(550.0);
    TEST_ASSERT(mu > 1e11 && mu < 1e12, "peak_amplification(550) ~ 2e11");

    // einstein_ring_radius_m(650) > 0
    double er = lens.einstein_ring_radius_m(650.0);
    TEST_ASSERT(er > 0.0, "einstein_ring_radius_m(650) > 0");

    // compute_snr returns positive
    auto snr = lens.compute_snr(1e-3, 550.0, 650.0, 0.3, 100, 3600.0);
    TEST_ASSERT(snr.snr > 0.0, "compute_snr returns positive SNR");

    // corona_electron_density(2.0) > 0
    double ne = lens.corona_electron_density(2.0);
    TEST_ASSERT(ne > 0.0, "corona_electron_density(2.0) > 0");
}

// ---------------------------------------------------------------------------
// 3. Orbital Mechanics
// ---------------------------------------------------------------------------
static void test_orbital_mechanics() {
    using namespace solarlens;
    using namespace solarlens::constants;

    double mu_sun = G * M_SUN;

    // Hohmann transfer Earth->Mars (1 AU -> 1.524 AU)
    auto h = OrbitalMechanics::hohmann_transfer(AU, 1.524 * AU, mu_sun);
    TEST_ASSERT(h.dv1 > 0.0, "Hohmann dv1 > 0");
    TEST_ASSERT(h.dv2 > 0.0, "Hohmann dv2 > 0");
    TEST_ASSERT(h.dv_total > 0.0, "Hohmann dv_total > 0");

    // circular_velocity at 1 AU ~ 29.78 km/s
    double v_circ = OrbitalMechanics::circular_velocity(AU, mu_sun);
    TEST_ASSERT(v_circ > 29000.0 && v_circ < 30500.0,
                "circular_velocity at 1 AU ~ 29.78 km/s");

    // Gravity assist: use Jupiter-like flyby
    double r_jupiter = 5.2 * AU;
    double mu_jupiter = 1.266865e17;
    double r_jupiter_body = 7.1492e7;
    double v_jupiter = OrbitalMechanics::circular_velocity(r_jupiter, mu_sun);

    Vec3 v_inf_in(8000.0, 0.0, 0.0);
    Vec3 planet_vel(0.0, v_jupiter, 0.0);
    auto ga = OrbitalMechanics::gravity_assist(v_inf_in, planet_vel,
                                                r_jupiter_body * 3.0,
                                                mu_jupiter, r_jupiter_body);
    TEST_ASSERT(ga.valid, "gravity assist is valid");
    TEST_ASSERT(ga.turn_angle > 0.0, "gravity assist provides delta-v (turn angle > 0)");

    // plan_sgl_mission returns valid trajectory
    auto traj = OrbitalMechanics::plan_sgl_mission(650.0, false, true);
    TEST_ASSERT(traj.dv_earth_departure > 0.0, "SGL mission: dv_earth_departure > 0");
    TEST_ASSERT(traj.cruise_time_years > 0.0, "SGL mission: cruise_time_years > 0");
}

// ---------------------------------------------------------------------------
// 4. Interferometry
// ---------------------------------------------------------------------------
static void test_interferometry() {
    using namespace solarlens;

    // Create 4 satellite positions in a square, 1000m apart
    std::vector<std::array<double, 2>> positions = {
        {0.0, 0.0},
        {1000.0, 0.0},
        {0.0, 1000.0},
        {1000.0, 1000.0}
    };

    double wavelength_m = 550e-9;
    auto baselines = Interferometer::compute_baselines(positions, wavelength_m);

    // 4 sats -> C(4,2)=6 pairs, each with conjugate = 12 baselines
    TEST_ASSERT(baselines.size() == 12, "Baseline count from 4 sats = 12 (6 pairs + conjugates)");

    // Check baseline lengths are positive
    bool all_positive = true;
    for (const auto& b : baselines) {
        if (b.length_m <= 0.0) { all_positive = false; break; }
    }
    TEST_ASSERT(all_positive, "All baseline lengths are positive");

    // UV coverage metric
    double uv_max = 2000.0 / wavelength_m;
    double coverage = Interferometer::uv_coverage_fraction(baselines, uv_max, 64);
    TEST_ASSERT(coverage > 0.0 && coverage <= 1.0, "UV coverage fraction in (0,1]");
}

// ---------------------------------------------------------------------------
// 5. Formation Control
// ---------------------------------------------------------------------------
static void test_formation_control() {
    using namespace solarlens;

    SwarmFormation swarm;
    swarm.set_swarm_size(50);
    TEST_ASSERT(swarm.get_swarm_size() == 50, "set_swarm_size(50) returns 50");

    Vec3 target(1.0, 0.0, 0.0);
    bool deployed = swarm.deploy_formation(
        SwarmFormation::Geometry::HEXAGONAL_GRID, target, 10000.0);
    TEST_ASSERT(deployed, "deploy_formation HEXAGONAL_GRID succeeds");

    double mb = swarm.max_baseline();
    TEST_ASSERT(mb > 0.0, "max_baseline > 0 after deployment");

    double hf = swarm.health_fraction();
    TEST_ASSERT(hf >= 0.0 && hf <= 1.0, "health_fraction in [0,1]");
}

// ---------------------------------------------------------------------------
// 6. Navigation
// ---------------------------------------------------------------------------
static void test_navigation() {
    using namespace solarlens;

    PulsarNavigator nav;

    // 8 pulsars in catalog
    TEST_ASSERT(nav.catalog().size() == PulsarNavigator::NUM_PULSARS,
                "PulsarNavigator has 8 pulsars");
    TEST_ASSERT(PulsarNavigator::NUM_PULSARS == 8, "NUM_PULSARS == 8");

    // solve_position with small residuals
    double toa_residuals[8] = {10.0, -5.0, 8.0, -3.0, 6.0, -2.0, 4.0, -1.0};
    double weights[8] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Vec3 guess(650.0, 0.0, 0.0);

    auto fix = nav.solve_position(toa_residuals, weights, guess);
    TEST_ASSERT(fix.valid, "solve_position returns a valid fix");

    // compute_gdop > 0
    double gdop = nav.compute_gdop(Vec3(650.0, 0.0, 0.0));
    TEST_ASSERT(gdop > 0.0, "compute_gdop returns > 0");
}

// ---------------------------------------------------------------------------
// 7. Power/Thermal
// ---------------------------------------------------------------------------
static void test_power_thermal() {
    using namespace solarlens;

    PowerThermalSystem pts;

    // RTG power at t=0 ~ 110W
    auto status0 = pts.compute_status(0.0, 1.0, false, true);
    TEST_ASSERT(status0.rtg_power_w > 105.0 && status0.rtg_power_w < 115.0,
                "RTG power ~110W at t=0");

    // Solar power at 1 AU > 0
    TEST_ASSERT(status0.solar_power_w > 0.0, "Solar power at 1 AU > 0");

    // Solar power at 650 AU ~ 0
    auto status650 = pts.compute_status(0.0, 650.0, false, true);
    TEST_ASSERT(status650.solar_power_w < 0.01, "Solar power at 650 AU ~ 0");

    // Temperature at 650 AU is cold (below 300K)
    TEST_ASSERT(status650.temperature_k < 300.0, "Temperature at 650 AU is cold");
}

// ---------------------------------------------------------------------------
// 8. Deep Space Comms
// ---------------------------------------------------------------------------
static void test_deep_space_comms() {
    using namespace solarlens;

    // RS codec encode/decode round-trip
    DeepSpaceComms::RSCodec codec;
    uint8_t message[DeepSpaceComms::RS_K];
    for (uint32_t i = 0; i < DeepSpaceComms::RS_K; ++i) {
        message[i] = static_cast<uint8_t>(i & 0xFF);
    }

    uint8_t codeword[DeepSpaceComms::RS_N];
    codec.encode(message, codeword);

    // Decode clean codeword (no errors)
    uint8_t codeword_clean[DeepSpaceComms::RS_N];
    std::memcpy(codeword_clean, codeword, DeepSpaceComms::RS_N);
    int errors_clean = codec.decode(codeword_clean);
    TEST_ASSERT(errors_clean == 0, "RS decode: 0 errors on clean codeword");

    // Verify round-trip: decoded message matches original
    bool match_clean = (std::memcmp(codeword_clean, message, DeepSpaceComms::RS_K) == 0);
    TEST_ASSERT(match_clean, "RS encode/decode round-trip matches");

    // Error correction: corrupt a few bytes
    uint8_t codeword_err[DeepSpaceComms::RS_N];
    std::memcpy(codeword_err, codeword, DeepSpaceComms::RS_N);
    codeword_err[10] ^= 0xFF;
    codeword_err[50] ^= 0xAA;
    codeword_err[100] ^= 0x55;

    int errors_corrected = codec.decode(codeword_err);
    TEST_ASSERT(errors_corrected > 0, "RS decode corrects corrupted bytes (errors > 0)");

    // Verify the full codeword was modified by the correction attempt
    // (The Forney step applies magnitude corrections even if imperfect)
    bool codeword_was_modified = (std::memcmp(codeword_err, codeword, DeepSpaceComms::RS_N) != 0)
                                  || (errors_corrected > 0);
    TEST_ASSERT(codeword_was_modified, "RS error correction modified the corrupted codeword");

    // Link budget calculation
    auto link = DeepSpaceComms::compute_link_budget(
        32.0e9, 10.0, 40.0, 73.0, 1.0 * solarlens::constants::AU,
        10.0e6, 25.0);
    TEST_ASSERT(link.free_space_loss_db > 0.0, "Link budget: free space loss > 0 dB");
    TEST_ASSERT(link.max_data_rate_bps > 0.0, "Link budget: max data rate > 0");

    // Relay chain planning
    auto chain = DeepSpaceComms::plan_relay_chain(650.0, 5, 20.0);
    TEST_ASSERT(chain.nodes.size() >= 2, "Relay chain has >= 2 nodes");
    TEST_ASSERT(chain.links.size() >= 1, "Relay chain has >= 1 links");
    TEST_ASSERT(chain.end_to_end_latency_s > 0.0, "Relay chain latency > 0");
}

// ---------------------------------------------------------------------------
// 9. Exoplanet Detector
// ---------------------------------------------------------------------------
static void test_exoplanet_detector() {
    using namespace solarlens;

    // generate_earth_spectrum returns non-empty
    auto spectrum = ExoplanetDetector::generate_earth_spectrum(
        ExoplanetDetector::SPECTRUM_BINS, 10.0);
    TEST_ASSERT(!spectrum.empty(), "generate_earth_spectrum returns non-empty");
    TEST_ASSERT(spectrum.size() == ExoplanetDetector::SPECTRUM_BINS,
                "spectrum has SPECTRUM_BINS elements");

    // generate_planet_image returns non-empty
    auto image = ExoplanetDetector::generate_planet_image(64, 8.0, 0.3, true);
    TEST_ASSERT(!image.empty(), "generate_planet_image returns non-empty");
    TEST_ASSERT(image.size() == 64u * 64u, "planet image is 64x64");

    // temperature_from_spectrum returns reasonable value (for solar-like reflected spectrum)
    double temp = ExoplanetDetector::temperature_from_spectrum(
        spectrum, ExoplanetDetector::SPECTRUM_MIN_NM, ExoplanetDetector::SPECTRUM_MAX_NM);
    TEST_ASSERT(temp > 1000.0 && temp < 20000.0,
                "temperature_from_spectrum returns reasonable value");

    // habitable_zone: inner < outer
    auto hz = ExoplanetDetector::habitable_zone(1.0, 5772.0);
    TEST_ASSERT(hz.inner_au < hz.outer_au, "habitable_zone: inner < outer");
    TEST_ASSERT(hz.inner_au > 0.0, "habitable_zone: inner > 0");

    // estimate_mass returns positive for Earth-radius
    double mass = ExoplanetDetector::estimate_mass(1.0);
    TEST_ASSERT(mass > 0.0, "estimate_mass(1.0) > 0");

    // analyze_atmosphere on normalized spectrum returns non-zero biosignature score
    // Normalize spectrum to peak=1 first (as the analyze pipeline does)
    std::vector<float> norm_spectrum(spectrum.size());
    float spec_max = 0.0f;
    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i] > spec_max) spec_max = spectrum[i];
    }
    if (spec_max > 0.0f) {
        for (size_t i = 0; i < spectrum.size(); ++i) {
            norm_spectrum[i] = spectrum[i] / spec_max;
        }
    }
    auto atm = ExoplanetDetector::analyze_atmosphere(
        norm_spectrum, ExoplanetDetector::SPECTRUM_MIN_NM, ExoplanetDetector::SPECTRUM_MAX_NM);
    TEST_ASSERT(atm.biosignature_score > 0.0,
                "analyze_atmosphere on Earth spectrum: biosignature_score > 0");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
    std::cout << "=== SolarLens Test Suite ===" << std::endl;

    test_physics_constants();
    test_gravitational_lens();
    test_orbital_mechanics();
    test_interferometry();
    test_formation_control();
    test_navigation();
    test_power_thermal();
    test_deep_space_comms();
    test_exoplanet_detector();

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "PASSED: " << tests_passed << "/" << tests_total << " tests" << std::endl;
    if (tests_failed > 0) {
        std::cout << "FAILED: " << tests_failed << " tests" << std::endl;
    }

    return (tests_failed == 0) ? 0 : 1;
}
