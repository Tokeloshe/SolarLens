#include "solarlens/physics/gravitational_lens.hpp"
#include "solarlens/physics/orbital_mechanics.hpp"
#include "solarlens/physics/interferometry.hpp"
#include "solarlens/imaging/exoplanet_detector.hpp"
#include "solarlens/spacecraft/formation_control.hpp"
#include "solarlens/spacecraft/navigation.hpp"
#include "solarlens/spacecraft/power_thermal.hpp"
#include "solarlens/communication/deep_space_relay.hpp"

#include <iostream>
#include <iomanip>
#include <cstring>

namespace solarlens {

// ============================================================================
// MAIN MISSION CONTROLLER
//
// 6-phase state machine controlling a Starlink-scale CubeSat swarm
// at the solar gravitational lens focal region (~650 AU).
// ============================================================================
class SolarLensMission {
private:
    GravitationalLens sgl_;
    SwarmFormation swarm_;
    PulsarNavigator navigator_;
    PowerThermalSystem power_;
    DeepSpaceComms::RSCodec rs_codec_;
    ExoplanetDetector detector_;

    enum class Phase : uint8_t {
        LAUNCH,
        CRUISE,
        ARRIVAL,
        FORMATION,
        OBSERVATION,
        TRANSMISSION
    } current_phase_;

    uint64_t mission_time_ms_;
    uint32_t cycle_count_;

public:
    SolarLensMission()
        : current_phase_(Phase::LAUNCH)
        , mission_time_ms_(0)
        , cycle_count_(0)
    {
        // Initialize swarm with 16 spacecraft for demo
        swarm_.set_swarm_size(16);
    }

    void execute() {
        mission_time_ms_ += 100000;  // Accelerated: 100s per cycle for demo
        cycle_count_++;

        // Print status periodically
        if (cycle_count_ % 50 == 0 && cycle_count_ > 0) {
            print_status();
        }

        // Use nominal distance (full XNAV only in observation phase)
        double distance_au = constants::SGL_FOCAL_NOMINAL_AU;

        switch (current_phase_) {
            case Phase::LAUNCH:
                if (cycle_count_ > 20) {  // After 20 cycles
                    std::cout << "\n=== TRANSITION TO CRUISE PHASE ===\n";
                    current_phase_ = Phase::CRUISE;
                }
                break;

            case Phase::CRUISE:
                // Demonstrate trajectory planning
                if (cycle_count_ == 200) {
                    std::cout << "\n--- Trajectory Analysis ---\n";
                    auto plan = OrbitalMechanics::plan_sgl_mission(
                        constants::SGL_FOCAL_NOMINAL_AU, true, true
                    );
                    std::cout << "Earth departure Δv: "
                              << std::fixed << std::setprecision(2)
                              << plan.dv_earth_departure / 1000.0 << " km/s\n";
                    std::cout << "Jupiter assist Δv gained: "
                              << plan.dv_jupiter_assist / 1000.0 << " km/s\n";
                    std::cout << "Total Δv: " << plan.dv_total / 1000.0 << " km/s\n";
                    std::cout << "Cruise time: " << plan.cruise_time_years << " years\n";
                    std::cout << "Arrival speed: " << plan.arrival_speed_km_s << " km/s\n";
                }
                if (cycle_count_ > 80) {
                    std::cout << "\n=== ARRIVAL AT FOCAL REGION (" << distance_au << " AU) ===\n";
                    current_phase_ = Phase::ARRIVAL;
                }
                break;

            case Phase::ARRIVAL:
                std::cout << "\n=== DEPLOYING FORMATION ===\n";
                current_phase_ = Phase::FORMATION;
                break;

            case Phase::FORMATION: {
                Vec3 target_dir{1.0, 0.0, 0.0};  // Toward Alpha Centauri
                bool deployed = swarm_.deploy_formation(
                    SwarmFormation::Geometry::HEXAGONAL_GRID,
                    target_dir,
                    50000.0  // 50 km baseline
                );
                if (deployed) {
                    std::cout << "Formation deployed: " << swarm_.get_swarm_size()
                              << " satellites\n";
                    std::cout << "Max baseline: " << std::scientific
                              << swarm_.max_baseline() << " m\n";
                    std::cout << "Formation error RMS: " << std::fixed << std::setprecision(3)
                              << swarm_.formation_error_rms() << " m\n";
                    std::cout << "Health fraction: " << std::setprecision(1)
                              << swarm_.health_fraction() * 100.0 << "%\n";

                    std::cout << "\n=== FORMATION COMPLETE - STARTING OBSERVATION ===\n";
                    current_phase_ = Phase::OBSERVATION;
                }
                break;
            }

            case Phase::OBSERVATION: {
                if (cycle_count_ % 50 == 0 && cycle_count_ > 100) {
                    run_observation_cycle(distance_au);
                }
                break;
            }

            case Phase::TRANSMISSION: {
                run_transmission_cycle(distance_au);
                break;
            }
        }
    }

private:
    void run_observation_cycle(double distance_au) {
        std::cout << "\n--- EXOPLANET OBSERVATION CYCLE ---\n";
        std::cout << "Target: Alpha Centauri system (4.37 ly)\n";
        std::cout << "Observer distance: " << std::fixed << std::setprecision(1)
                  << distance_au << " AU\n";

        // SGL physics
        const double wavelength_nm = 550.0;
        double gain = sgl_.peak_amplification(wavelength_nm);
        double focal_dist = sgl_.focal_distance_au(wavelength_nm);
        double einstein_r = sgl_.einstein_ring_radius_m(distance_au);

        std::cout << "Wavelength: " << wavelength_nm << " nm\n";
        std::cout << "SGL peak amplification: " << std::scientific << std::setprecision(3)
                  << gain << "\n";
        std::cout << "Focal distance: " << std::fixed << std::setprecision(1)
                  << focal_dist << " AU\n";
        std::cout << "Einstein ring radius: " << std::setprecision(1) << einstein_r << " m\n";

        // SNR calculation
        double planet_flux = 1e-4;  // photons/m²/s at 1 AU for Earth twin at 4.37 ly
        auto snr = sgl_.compute_snr(
            planet_flux, wavelength_nm, distance_au,
            constants::SAT_APERTURE_M,
            swarm_.get_swarm_size(),
            3600.0  // 1 hour integration
        );

        std::cout << "\nSNR Analysis (1 hour integration, " << swarm_.get_swarm_size() << " sats):\n";
        std::cout << "  Signal photons: " << std::scientific << std::setprecision(2)
                  << snr.signal_photons << "\n";
        std::cout << "  Corona photons: " << snr.corona_photons << "\n";
        std::cout << "  SNR: " << std::fixed << std::setprecision(1) << snr.snr << "\n";

        // Generate synthetic planet data for demonstration
        if (snr.snr > 5.0 || cycle_count_ > 100) {
            std::cout << "\n*** EXOPLANET CANDIDATE DETECTED ***\n";

            // Generate synthetic Earth spectrum and image for pipeline demo
            auto spectrum = ExoplanetDetector::generate_earth_spectrum(
                ExoplanetDetector::SPECTRUM_BINS, 4.37);
            auto planet_img = ExoplanetDetector::generate_planet_image(64, 8.0, 0.3, true);

            // Run full analysis pipeline
            auto planet = detector_.analyze(
                planet_img, 64, 1e-9,
                spectrum, 1.0, 4.37, distance_au
            );

            if (planet.detected) {
                print_planet_data(planet);
                std::cout << "\nPreparing to transmit discovery...\n";
                current_phase_ = Phase::TRANSMISSION;
            }
        }
    }

    void run_transmission_cycle(double distance_au) {
        std::cout << "\n--- TRANSMITTING DISCOVERY ---\n";

        // Encode message with Reed-Solomon
        const char* msg = "EXOPLANET BIOSIGNATURE DETECTED - ALPHA CENTAURI";
        uint8_t message[DeepSpaceComms::RS_K] = {};
        std::memcpy(message, msg, std::min(std::strlen(msg), static_cast<size_t>(DeepSpaceComms::RS_K)));

        uint8_t codeword[DeepSpaceComms::RS_N] = {};
        rs_codec_.encode(message, codeword);

        std::cout << "RS(" << DeepSpaceComms::RS_N << "," << DeepSpaceComms::RS_K
                  << ") encoded: " << std::strlen(msg) << " bytes -> "
                  << DeepSpaceComms::RS_N << " bytes\n";

        // Simulate errors and correction
        codeword[10] ^= 0xFF;
        codeword[50] ^= 0xAB;
        codeword[100] ^= 0x55;
        int corrected = rs_codec_.decode(codeword);
        std::cout << "Simulated 3 symbol errors, corrected: " << corrected << "\n";

        // Link budget: satellite to Earth via DSN
        auto link_dsn = DeepSpaceComms::satellite_to_earth_link(distance_au);
        print_link_budget("Direct Ka-band to DSN", link_dsn);

        // Link budget: SGL-amplified (using Sun as antenna for return signal)
        auto link_sgl = DeepSpaceComms::sgl_amplified_link(4.37, 550.0);
        print_link_budget("SGL-amplified optical", link_sgl);

        // Relay chain
        auto relay = DeepSpaceComms::plan_relay_chain(distance_au, 5, 10.0);
        std::cout << "\nRelay Chain (" << relay.nodes.size() << " nodes):\n";
        std::cout << "  End-to-end rate: " << std::scientific << std::setprecision(2)
                  << relay.end_to_end_rate_bps << " bps\n";
        std::cout << "  End-to-end latency: " << std::fixed << std::setprecision(1)
                  << relay.end_to_end_latency_s / 3600.0 << " hours\n";
        std::cout << "  Total power: " << relay.total_power_w << " W\n";

        if (link_dsn.link_margin_db > 0 || relay.end_to_end_rate_bps > 0) {
            std::cout << "\n=== TRANSMISSION SUCCESSFUL ===\n";
            std::cout << "Returning to observation mode...\n\n";
            current_phase_ = Phase::OBSERVATION;
        }
    }

    void print_status() {
        std::cout << "\n--- Mission Time: " << std::fixed << std::setprecision(1)
                  << (mission_time_ms_ / 1000.0) << "s | Phase: "
                  << get_phase_name() << " ---\n";

        // SGL physics summary
        const double focal_dist = sgl_.focal_distance_au(550.0);
        const double gain = sgl_.peak_amplification(550.0);
        std::cout << "SGL focal distance (550nm): " << std::fixed << std::setprecision(1)
                  << focal_dist << " AU | Peak gain: " << std::scientific << std::setprecision(2)
                  << gain << "\n";

        // Navigation
        std::cout << "XNAV pulsars: " << PulsarNavigator::NUM_PULSARS
                  << " | GDOP: " << std::fixed << std::setprecision(2)
                  << navigator_.compute_gdop(Vec3{constants::SGL_FOCAL_NOMINAL_AU, 0, 0}) << "\n";

        // Power at 650 AU
        auto pwr = power_.compute_status(mission_time_ms_ / 1000.0,
                                         constants::SGL_FOCAL_NOMINAL_AU,
                                         false, true);
        std::cout << "Power: RTG " << std::fixed << std::setprecision(1)
                  << pwr.rtg_power_w << "W + Solar " << pwr.solar_power_w
                  << "W | Temp: " << pwr.temperature_k << "K"
                  << (pwr.heaters_active ? " [HEATERS ON]" : "") << "\n";
    }

    void print_planet_data(const ExoplanetDetector::PlanetData& planet) {
        std::cout << "\nPlanet Characterization:\n";
        std::cout << "  Radius: " << std::fixed << std::setprecision(2)
                  << planet.radius_earth << " Earth radii\n";
        std::cout << "  Mass: " << planet.mass_earth << " Earth masses\n";
        std::cout << "  Temperature: " << std::setprecision(0)
                  << planet.temperature_k << " K\n";
        std::cout << "  Albedo: " << std::setprecision(2) << planet.albedo << "\n";
        std::cout << "  Surface gravity: " << std::setprecision(1)
                  << planet.surface_gravity << " m/s²\n";
        std::cout << "  Orbital radius: " << std::setprecision(2)
                  << planet.orbital_radius_au << " AU\n";
        std::cout << "  Orbital period: " << std::setprecision(1)
                  << planet.orbital_period_days << " days\n";
        std::cout << "  Habitable zone: " << (planet.in_habitable_zone ? "YES" : "NO") << "\n";
        std::cout << "  Confidence: " << std::setprecision(1)
                  << planet.confidence * 100.0f << "%\n";

        std::cout << "\nAtmospheric Composition (volume mixing ratios):\n";
        std::cout << "  N2:  " << std::setprecision(3) << planet.atmosphere.n2 << "\n";
        std::cout << "  O2:  " << planet.atmosphere.o2 << "\n";
        std::cout << "  CO2: " << planet.atmosphere.co2 << "\n";
        std::cout << "  H2O: " << planet.atmosphere.h2o << "\n";
        std::cout << "  CH4: " << planet.atmosphere.ch4 << "\n";
        std::cout << "  O3:  " << planet.atmosphere.o3 << "\n";

        std::cout << "\nBiosignature Score: " << std::setprecision(1)
                  << planet.atmosphere.biosignature_score * 100.0 << "%\n";
        std::cout << "Assessment: " << planet.atmosphere.biosignature_rationale << "\n";
    }

    void print_link_budget(const char* label, const DeepSpaceComms::LinkBudget& link) {
        std::cout << "\n" << label << " Link Budget:\n";
        std::cout << "  Frequency: " << std::scientific << std::setprecision(2)
                  << link.frequency_hz << " Hz\n";
        std::cout << "  TX Power: " << std::fixed << std::setprecision(1)
                  << link.tx_power_w << " W\n";
        std::cout << "  Path Loss: " << link.free_space_loss_db << " dB\n";
        std::cout << "  SNR: " << link.snr_db << " dB\n";
        std::cout << "  Achievable Rate: " << std::scientific << std::setprecision(2)
                  << link.achievable_rate_bps << " bps\n";
        std::cout << "  Link Margin: " << std::fixed << std::setprecision(1)
                  << link.link_margin_db << " dB\n";
    }

    const char* get_phase_name() const {
        switch (current_phase_) {
            case Phase::LAUNCH:       return "LAUNCH";
            case Phase::CRUISE:       return "CRUISE";
            case Phase::ARRIVAL:      return "ARRIVAL";
            case Phase::FORMATION:    return "FORMATION";
            case Phase::OBSERVATION:  return "OBSERVATION";
            case Phase::TRANSMISSION: return "TRANSMISSION";
            default:                  return "UNKNOWN";
        }
    }
};

} // namespace solarlens

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================
int main() {
    std::cout << "==================================================\n";
    std::cout << "  SOLARLENS — Solar Gravitational Lens Mission\n";
    std::cout << "  Starlink-Scale CubeSat Swarm Control System\n";
    std::cout << "==================================================\n\n";

    // Print key SGL parameters
    std::cout << "Physical Constants:\n";
    std::cout << "  Schwarzschild radius: " << std::fixed << std::setprecision(2)
              << solarlens::constants::R_SCHWARZSCHILD << " m\n";
    std::cout << "  SGL focal minimum: " << std::setprecision(1)
              << solarlens::constants::SGL_FOCAL_MIN_AU << " AU\n";
    std::cout << "  SGL nominal distance: "
              << solarlens::constants::SGL_FOCAL_NOMINAL_AU << " AU\n\n";

    solarlens::SolarLensMission mission;

    std::cout << "Mission initialized. Running control loop (500 cycles)...\n\n";
    std::cout.flush();

    for (int i = 0; i < 500; ++i) {
        mission.execute();
    }

    std::cout << "\n==================================================\n";
    std::cout << "  Mission demonstration complete\n";
    std::cout << "==================================================\n";

    return 0;
}
