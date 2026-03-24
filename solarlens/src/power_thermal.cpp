#include "solarlens/spacecraft/power_thermal.hpp"
#include <algorithm>

namespace solarlens {

PowerThermalSystem::PowerThermalSystem()
    : cfg_{}, battery_soc_(0.8) {}

PowerThermalSystem::PowerThermalSystem(const Config& cfg)
    : cfg_(cfg), battery_soc_(0.8) {}

double PowerThermalSystem::rtg_power_at(double mission_elapsed_s) const {
    // Pu-238 exponential decay: P(t) = P_0 × exp(-λt)
    // λ = ln(2) / t_half
    const double t_half_s = cfg_.rtg_half_life_years * 365.25 * 86400.0;
    const double lambda = 0.693147 / t_half_s;
    return cfg_.rtg_initial_power_w * std::exp(-lambda * mission_elapsed_s);
}

double PowerThermalSystem::solar_power_at(double distance_au, double mission_elapsed_s) const {
    // Solar irradiance: 1361 W/m² at 1 AU, scales as 1/r²
    // Negligible beyond ~5 AU
    if (distance_au > 10.0) return 0.0;

    const double irradiance = 1361.0 / (distance_au * distance_au);
    const double years = mission_elapsed_s / (365.25 * 86400.0);
    const double degradation = std::pow(1.0 - cfg_.panel_degradation_per_year, years);

    return irradiance * cfg_.solar_panel_area_m2 * cfg_.solar_cell_eff * degradation;
}

double PowerThermalSystem::compute_temperature(double distance_au, double internal_heat_w) const {
    // Radiative equilibrium: absorbed power = emitted power
    //
    // Absorbed: Q_solar + Q_internal
    // Q_solar = α × S(r) × A_cross (solar absorptivity × irradiance × cross-section)
    //   where A_cross ≈ A_surface / 6 for a cube
    // Q_internal = RTG waste heat + electronics dissipation
    //
    // Emitted: ε σ A T⁴

    const double solar_irradiance = 1361.0 / (distance_au * distance_au);  // W/m²
    const double cross_section = cfg_.spacecraft_surface_area_m2 / 6.0;

    const double q_solar = cfg_.absorptivity * solar_irradiance * cross_section;
    const double q_total = q_solar + internal_heat_w;

    // T = (Q / (ε σ A))^(1/4)
    const double t4 = q_total / (cfg_.emissivity * constants::SIGMA_SB * cfg_.spacecraft_surface_area_m2);
    if (t4 <= 0) return 2.7;  // CMB temperature minimum
    return std::pow(t4, 0.25);
}

void PowerThermalSystem::update_battery(double dt_s, double power_surplus_w) {
    // power_surplus_w > 0 means charging, < 0 means discharging
    if (power_surplus_w > 0) {
        double energy_in = power_surplus_w * dt_s / 3600.0;  // Wh
        battery_soc_ += (energy_in * cfg_.charge_efficiency) / cfg_.battery_capacity_wh;
        battery_soc_ = std::min(battery_soc_, cfg_.max_soc);
    } else {
        double energy_out = -power_surplus_w * dt_s / 3600.0;
        battery_soc_ -= energy_out / (cfg_.battery_capacity_wh * cfg_.discharge_efficiency);
        battery_soc_ = std::max(battery_soc_, 0.0);
    }
}

PowerThermalSystem::SystemStatus PowerThermalSystem::compute_status(
    double mission_elapsed_s,
    double distance_au,
    bool is_thrusting,
    bool camera_active
) {
    SystemStatus status{};

    // Power generation
    status.rtg_power_w = rtg_power_at(mission_elapsed_s);
    status.solar_power_w = solar_power_at(distance_au, mission_elapsed_s);
    status.total_generation_w = status.rtg_power_w + status.solar_power_w;

    // Power consumption
    PowerBudget budget;
    if (is_thrusting) {
        budget.propulsion_w = 150.0;  // Ion engine power
    }
    if (!camera_active) {
        budget.camera_w = 0.5;  // Standby
    }

    // RTG waste heat (thermal power minus electrical)
    double rtg_thermal = cfg_.rtg_thermal_w * std::exp(
        -0.693147 * mission_elapsed_s / (cfg_.rtg_half_life_years * 365.25 * 86400.0));
    double rtg_waste = rtg_thermal - status.rtg_power_w;

    // Temperature calculation
    double internal_heat = rtg_waste + cfg_.internal_dissipation_w;
    status.temperature_k = compute_temperature(distance_au, internal_heat);

    // Heater power needed if below operational minimum
    status.heaters_active = false;
    status.heater_power_needed_w = 0;
    if (status.temperature_k < cfg_.t_min_operational_k) {
        // Need to add heat to reach minimum operational temperature
        // ΔQ = ε σ A (T_target⁴ - T_current⁴)
        double t_target4 = std::pow(cfg_.t_min_operational_k, 4.0);
        double t_current4 = std::pow(status.temperature_k, 4.0);
        status.heater_power_needed_w = cfg_.emissivity * constants::SIGMA_SB *
                                       cfg_.spacecraft_surface_area_m2 * (t_target4 - t_current4);
        status.heater_power_needed_w = std::min(status.heater_power_needed_w, cfg_.heater_power_w);
        status.heaters_active = true;
        status.temperature_k = cfg_.t_min_operational_k;  // Heaters maintain minimum
    }
    budget.heater_w = status.heater_power_needed_w;

    status.total_consumption_w = budget.total();
    status.power_margin_w = status.total_generation_w - status.total_consumption_w;

    // Battery state
    status.battery_soc = battery_soc_;
    if (status.total_consumption_w > 0) {
        status.battery_runtime_hours = (battery_soc_ * cfg_.battery_capacity_wh) /
                                       status.total_consumption_w;
    }

    // Status flags
    status.low_power_warning = battery_soc_ < 0.30;
    status.critical_power = battery_soc_ < 0.15;
    status.thermal_warning = (status.temperature_k < cfg_.t_min_operational_k ||
                              status.temperature_k > cfg_.t_max_operational_k);
    status.thermal_critical = (status.temperature_k < cfg_.t_min_survival_k ||
                               status.temperature_k > cfg_.t_max_survival_k);
    status.nominal = !status.low_power_warning && !status.thermal_warning &&
                     status.power_margin_w > 0;

    return status;
}

} // namespace solarlens
