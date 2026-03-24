#pragma once

#include "solarlens/physics/constants.hpp"
#include <cstdint>
#include <cmath>

namespace solarlens {

// ============================================================================
// POWER & THERMAL MANAGEMENT FOR DEEP SPACE SWARM
//
// Models power generation, storage, and thermal balance for satellites
// operating from 1 AU to 900 AU. Power sources include:
// - Solar panels (useful to ~5 AU, negligible beyond Jupiter)
// - Radioisotope Thermoelectric Generator (RTG) — primary beyond 5 AU
// - Next-gen: betavoltaic or Am-241 based sources
//
// Thermal model uses radiative equilibrium with internal heat sources.
// ============================================================================

class PowerThermalSystem {
public:
    // Configuration for a single satellite's power system
    struct Config {
        // RTG parameters (based on MMRTG / Next-Gen RTG)
        double rtg_initial_power_w  = 110.0;    // Beginning-of-life electrical power
        double rtg_thermal_w        = 2000.0;    // Total thermal output (Pu-238)
        double rtg_efficiency       = 0.065;     // Thermoelectric conversion efficiency
        double rtg_half_life_years  = 87.7;      // Pu-238 half-life
        double rtg_mass_kg          = 45.0;      // RTG mass

        // Solar panel parameters
        double solar_panel_area_m2  = 12.0;      // Panel area (Starlink-scale)
        double solar_cell_eff       = 0.32;       // Triple-junction GaAs efficiency
        double panel_degradation_per_year = 0.01; // 1% per year radiation damage

        // Battery (Li-ion space-grade)
        double battery_capacity_wh  = 2400.0;    // 100 Ah × 24V
        double battery_mass_kg      = 15.0;
        double charge_efficiency    = 0.95;
        double discharge_efficiency = 0.97;
        double min_soc              = 0.10;       // Minimum state of charge
        double max_soc              = 0.95;       // Maximum to preserve lifetime

        // Thermal
        double spacecraft_surface_area_m2 = 6.0;  // Total radiating surface
        double emissivity           = 0.85;        // Thermal coating emissivity
        double absorptivity         = 0.20;        // Solar absorptivity (white paint)
        double internal_dissipation_w = 15.0;      // Electronics waste heat
        double heater_power_w       = 20.0;        // Survival heaters max power

        // Temperature limits
        double t_min_operational_k  = 233.0;      // -40°C
        double t_max_operational_k  = 343.0;      // +70°C
        double t_min_survival_k     = 213.0;      // -60°C
        double t_max_survival_k     = 358.0;      // +85°C
    };

    // Power consumption budget for one satellite
    struct PowerBudget {
        double comms_w        = 25.0;    // Laser inter-sat link + downlink
        double computer_w     = 15.0;    // On-board processing
        double camera_w       = 8.0;     // Telescope + detector cooling
        double attitude_w     = 5.0;     // Reaction wheels + star tracker
        double propulsion_w   = 0.0;     // Ion engine (when thrusting)
        double heater_w       = 0.0;     // Computed dynamically
        double margin_w       = 5.0;     // 10% margin

        double total() const {
            return comms_w + computer_w + camera_w + attitude_w +
                   propulsion_w + heater_w + margin_w;
        }
    };

    struct SystemStatus {
        // Power
        double rtg_power_w;
        double solar_power_w;
        double total_generation_w;
        double total_consumption_w;
        double power_margin_w;
        double battery_soc;               // State of charge (0-1)
        double battery_runtime_hours;     // Time at current draw

        // Thermal
        double temperature_k;
        double heater_power_needed_w;
        bool heaters_active;

        // Status flags
        bool nominal;                     // All systems within limits
        bool low_power_warning;           // SOC < 30%
        bool critical_power;              // SOC < 15%
        bool thermal_warning;             // Outside operational range
        bool thermal_critical;            // Outside survival range
    };

    PowerThermalSystem();
    explicit PowerThermalSystem(const Config& cfg);

    // Compute full power and thermal status
    // mission_elapsed_s: seconds since launch
    // distance_au: current distance from Sun
    // is_thrusting: whether ion engine is firing
    // camera_active: whether detector is powered
    SystemStatus compute_status(
        double mission_elapsed_s,
        double distance_au,
        bool is_thrusting,
        bool camera_active
    );

    // Step the battery state forward by dt seconds
    void update_battery(double dt_s, double power_surplus_w);

    // Get current battery SOC
    double battery_soc() const { return battery_soc_; }

    const Config& config() const { return cfg_; }

private:
    Config cfg_;
    double battery_soc_;          // Current state of charge (0-1)

    // Compute equilibrium temperature from heat balance
    double compute_temperature(double distance_au, double internal_heat_w) const;

    // Solar power at given distance
    double solar_power_at(double distance_au, double mission_elapsed_s) const;

    // RTG power at given mission time
    double rtg_power_at(double mission_elapsed_s) const;
};

} // namespace solarlens
