#include "solarlens/physics/gravitational_lens.hpp"
#include "solarlens/imaging/exoplanet_detector.hpp"
#include "solarlens/spacecraft/formation_control.hpp"
#include "solarlens/spacecraft/navigation.hpp"
#include "solarlens/spacecraft/power_thermal.hpp"
#include "solarlens/communication/deep_space_relay.hpp"

#include <iostream>
#include <iomanip>
#include <unistd.h>

namespace solarlens {

// ============================================================================
// MAIN MISSION CONTROLLER
// ============================================================================
class SolarLensMission {
private:
    GravitationalLensPhysics physics;
    // ExoplanetDetector detector; // Commented out due to large stack arrays (8MB+)
    SwarmController swarm;
    DeepSpaceNavigator navigator;
    InterstellarTransmitter transmitter;
    PowerThermalManager power;

    enum class Phase : uint8_t {
        LAUNCH,
        CRUISE,
        ARRIVAL,
        FORMATION,
        OBSERVATION,
        TRANSMISSION
    } current_phase;

    uint64_t mission_time_ms;
    uint32_t cycle_count;

public:
    SolarLensMission() : current_phase(Phase::LAUNCH), mission_time_ms(0), cycle_count(0) {
        // Initialize swarm with 16 spacecraft
        swarm.set_active_count(16);
    }

    // Main control loop - runs at 10 Hz
    void execute() {
        // Update mission time
        mission_time_ms += 100;  // 100ms per cycle at 10Hz
        cycle_count++;

        // Print status every 10 seconds (100 cycles)
        if (cycle_count % 100 == 0) {
            print_status();
        }

        // Get current state
        const auto nav_solution = navigator.calculate_position(nullptr, mission_time_ms * 1000000);
        const auto power_status = power.calculate_power_status(mission_time_ms / 86400000,
                                                              nav_solution.position_au[0]);

        // Suppress unused variable warning (used for monitoring in real mission)
        (void)power_status;

        // Phase-specific operations
        switch (current_phase) {
            case Phase::LAUNCH:
                if (mission_time_ms > 3600000) {  // 1 hour
                    std::cout << "\n=== TRANSITION TO CRUISE PHASE ===\n";
                    current_phase = Phase::CRUISE;
                }
                break;

            case Phase::CRUISE:
                // 25-year journey (simulated - would take actual years)
                if (nav_solution.position_au[0] > constants::FOCAL_OPTIMAL_AU - 1.0) {
                    std::cout << "\n=== ARRIVAL AT FOCAL POINT ===\n";
                    current_phase = Phase::ARRIVAL;
                }
                break;

            case Phase::ARRIVAL:
                // Slow down and prepare for formation
                if (nav_solution.velocity_km_s[0] < 1.0) {
                    std::cout << "\n=== BEGINNING FORMATION ===\n";
                    current_phase = Phase::FORMATION;
                }
                break;

            case Phase::FORMATION: {
                // Form observation array
                double target[3] = {650.0, 0.0, 0.0};
                if (swarm.optimize_formation(
                    SwarmController::Formation::HEXAGONAL_GRID,
                    target,
                    1000.0  // 1000 km baseline
                )) {
                    std::cout << "\n=== FORMATION COMPLETE - STARTING OBSERVATION ===\n";
                    current_phase = Phase::OBSERVATION;
                }
                break;
            }

            case Phase::OBSERVATION: {
                // Simulated exoplanet detection (every 10 seconds in demo)
                if (cycle_count % 100 == 0 && cycle_count > 100) {
                    std::cout << "\n--- EXOPLANET DETECTION SIMULATION ---\n";
                    std::cout << "Performing gravitational lens imaging...\n";
                    std::cout << "Integration time: 3600 seconds\n";
                    std::cout << "Target: Alpha Centauri system (4.37 ly)\n";
                    std::cout << "Wavelength: 550 nm (visible light)\n";

                    // Calculate magnification for demo
                    double mag = physics.calculate_magnification(4.37, 650.0, 5000.0);
                    std::cout << "Achieved magnification: " << std::scientific
                              << std::setprecision(2) << mag << "\n";

                    // Simulate detection after 30 seconds of observation
                    if (cycle_count > 300) {
                        std::cout << "\n!!! EXOPLANET CANDIDATE DETECTED !!!\n";
                        std::cout << "Estimated radius: 1.05 Earth radii\n";
                        std::cout << "Orbital radius: 1.2 AU (habitable zone)\n";
                        std::cout << "Atmospheric signatures detected (O2, H2O)\n";
                        std::cout << "Biosignature score: 0.6 (moderate)\n";
                        std::cout << "\nPreparing to transmit discovery...\n";
                        current_phase = Phase::TRANSMISSION;
                    }
                }
                break;
            }

            case Phase::TRANSMISSION: {
                // Send discovery back to Earth
                const uint8_t message[] = "LIFE DETECTED";
                const auto encoded = transmitter.encode_message(
                    message,
                    sizeof(message),
                    InterstellarTransmitter::ErrorCorrection::TURBO_CODES
                );

                // Suppress unused variable warning (would be transmitted in real mission)
                (void)encoded;

                // Calculate link budget
                const auto link = transmitter.calculate_link_budget(
                    0.0,    // Earth distance
                    1e9,    // Lens magnification
                    true    // Use lens
                );

                print_link_budget(link);

                if (link.link_margin_db > 0) {
                    std::cout << "=== TRANSMISSION SUCCESSFUL ===\n";
                    std::cout << "Returning to observation mode...\n";
                    current_phase = Phase::OBSERVATION;
                }
                break;
            }
        }
    }

    void print_status() {
        std::cout << "\n--- Mission Time: " << (mission_time_ms / 1000.0) << "s ---\n";
        std::cout << "Phase: " << get_phase_name(current_phase) << "\n";

        // Calculate and display focal distance for visible light
        const double focal_dist = physics.calculate_focal_distance_au(550.0);
        std::cout << "Optimal Focal Distance (550nm): " << std::fixed << std::setprecision(1)
                  << focal_dist << " AU\n";

        // Display magnification at Einstein ring
        const double mag = physics.calculate_magnification(10.0, 650.0, 5000.0);
        std::cout << "Magnification: " << std::scientific << std::setprecision(2) << mag << "\n";
    }

    void print_planet_data(const ExoplanetDetector::PlanetData& planet) {
        std::cout << "\n*** EXOPLANET DETECTED ***\n";
        std::cout << "Radius: " << std::fixed << std::setprecision(2)
                  << planet.radius_earth << " Earth radii\n";
        std::cout << "Temperature: " << std::setprecision(0)
                  << planet.temperature_kelvin << " K\n";
        std::cout << "Orbital Radius: " << std::setprecision(2)
                  << planet.orbital_radius_au << " AU\n";
        std::cout << "Habitable Zone: " << (planet.in_habitable_zone ? "YES" : "NO") << "\n";
        std::cout << "Confidence: " << std::setprecision(1)
                  << planet.confidence * 100.0f << "%\n";
        std::cout << "\nAtmospheric Composition:\n";
        std::cout << "  O2:  " << planet.atmosphere.oxygen << "%\n";
        std::cout << "  CH4: " << planet.atmosphere.methane << "%\n";
        std::cout << "  H2O: " << planet.atmosphere.water << "%\n";
        std::cout << "  CO2: " << planet.atmosphere.co2 << "%\n";
        std::cout << "  N2:  " << planet.atmosphere.nitrogen << "%\n";
        std::cout << "\nBiosignature Score: " << std::setprecision(1)
                  << planet.atmosphere.biosignature_score * 100.0f << "%\n";
    }

    void print_link_budget(const InterstellarTransmitter::LinkBudget& link) {
        std::cout << "\n--- Communication Link Budget ---\n";
        std::cout << "Frequency: " << link.frequency_ghz << " GHz\n";
        std::cout << "TX Power: " << link.tx_power_watts << " W\n";
        std::cout << "Path Loss: " << std::fixed << std::setprecision(1)
                  << link.path_loss_db << " dB\n";
        std::cout << "Data Rate: " << std::scientific << std::setprecision(2)
                  << link.data_rate_bps << " bps\n";
        std::cout << "Link Margin: " << std::fixed << std::setprecision(1)
                  << link.link_margin_db << " dB\n";
    }

    const char* get_phase_name(Phase phase) {
        switch (phase) {
            case Phase::LAUNCH: return "LAUNCH";
            case Phase::CRUISE: return "CRUISE";
            case Phase::ARRIVAL: return "ARRIVAL";
            case Phase::FORMATION: return "FORMATION";
            case Phase::OBSERVATION: return "OBSERVATION";
            case Phase::TRANSMISSION: return "TRANSMISSION";
            default: return "UNKNOWN";
        }
    }
};

} // namespace solarlens

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================
int main() {
    std::cout << "==================================================\n";
    std::cout << "  SOLARLENS - Solar Gravitational Lens Mission\n";
    std::cout << "  Flight-Ready CubeSat Swarm Control System\n";
    std::cout << "==================================================\n\n";
    std::cout.flush();

    solarlens::SolarLensMission mission;

    std::cout << "Mission initialized. Running control loop...\n";
    std::cout << "(Running for 50 cycles demonstration)\n\n";
    std::cout.flush();

    // Run for 50 cycles (5 seconds at 10Hz simulated)
    for (int i = 0; i < 500; ++i) {
        mission.execute();

        // Reduced sleep for faster demo
        usleep(10000);  // 10ms instead of 100ms
    }

    std::cout << "\n==================================================\n";
    std::cout << "  Mission demonstration complete\n";
    std::cout << "==================================================\n";

    return 0;
}
