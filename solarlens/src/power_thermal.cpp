#include "solarlens/spacecraft/power_thermal.hpp"
#include <cmath>
#include <algorithm>

namespace solarlens {

PowerThermalManager::PowerStatus PowerThermalManager::calculate_power_status(
    uint32_t mission_days,
    double distance_au
) {
    PowerStatus status{};

    // RTG output with decay
    const double years_elapsed = mission_days / 365.0;
    status.rtg_output_w = rtg.initial_power_w *
                         std::exp(-rtg.decay_rate_per_year * years_elapsed);

    // Solar panels (negligible beyond Jupiter)
    const double solar_power_w = (distance_au < 10.0) ?
                                100.0 / (distance_au * distance_au) : 0.0;

    // Total available power
    status.available_power_w = status.rtg_output_w + solar_power_w;

    // Power consumption budget
    const double comm_power = 3.0;      // Communications
    const double compute_power = 2.0;   // Processing
    const double sensors_power = 2.0;   // Instruments
    const double thermal_power = 1.0;   // Heaters
    const double attitude_power = 1.0;  // Reaction wheels

    status.power_consumption_w = comm_power + compute_power + sensors_power +
                                thermal_power + attitude_power;

    // Battery state
    if (status.available_power_w > status.power_consumption_w) {
        // Charging
        const double charge_rate = status.available_power_w - status.power_consumption_w;
        battery.current_soc += (charge_rate * battery.charge_efficiency) / battery.capacity_wh;
        battery.current_soc = std::min(battery.current_soc, 1.0);
    } else {
        // Discharging
        const double discharge_rate = status.power_consumption_w - status.available_power_w;
        battery.current_soc -= discharge_rate / (battery.capacity_wh * battery.discharge_efficiency);
        battery.current_soc = std::max(battery.current_soc, 0.0);
    }

    status.battery_soc = battery.current_soc;
    status.battery_runtime_hours = (battery.current_soc * battery.capacity_wh) /
                                  status.power_consumption_w;

    // Temperature calculation (simplified)
    const double solar_heating_w = 1361.0 / (distance_au * distance_au);  // Solar constant
    const double rtg_heating_w = rtg.heat_output_w * (1.0 - rtg.efficiency);
    const double total_heating_w = solar_heating_w + rtg_heating_w;

    // Radiative cooling (Stefan-Boltzmann)
    const double surface_area_m2 = 0.1;  // CubeSat surface
    const double emissivity = 0.9;
    const double sigma = 5.67e-8;

    // Equilibrium temperature
    status.temperature_k = std::pow(total_heating_w /
                                   (emissivity * sigma * surface_area_m2), 0.25);

    // Power modes
    status.low_power_mode = status.battery_soc < 0.3;
    status.critical_power = status.battery_soc < 0.1;

    return status;
}

} // namespace solarlens
