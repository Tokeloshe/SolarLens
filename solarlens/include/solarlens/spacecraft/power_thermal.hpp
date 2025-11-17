#pragma once

#include "solarlens/physics/constants.hpp"
#include <cstdint>
#include <cmath>

namespace solarlens {

// ============================================================================
// POWER & THERMAL MANAGEMENT
// ============================================================================
class PowerThermalManager {
private:
    // Radioisotope Thermoelectric Generator model
    struct RTG {
        double initial_power_w = 10.0;
        double decay_rate_per_year = 0.02;
        double efficiency = 0.07;
        double heat_output_w = 140.0;
    } rtg;

    // Battery model
    struct Battery {
        double capacity_wh = 100.0;
        double charge_efficiency = 0.95;
        double discharge_efficiency = 0.98;
        double current_soc = 0.8;  // State of charge
        double temperature_k = 273.0;
    } battery;

public:
    struct PowerStatus {
        double available_power_w;
        double power_consumption_w;
        double battery_soc;
        double battery_runtime_hours;
        double rtg_output_w;
        double temperature_k;
        bool low_power_mode;
        bool critical_power;
    };

    [[nodiscard]] PowerStatus calculate_power_status(
        uint32_t mission_days,
        double distance_au
    );
};

} // namespace solarlens
