#include "solarlens/communication/deep_space_relay.hpp"
#include <cmath>

namespace solarlens {

InterstellarTransmitter::LinkBudget InterstellarTransmitter::calculate_link_budget(
    double distance_ly,
    double lens_magnification,
    bool use_lens
) {
    LinkBudget link{};

    // Transmission parameters
    link.frequency_ghz = 32.0;        // Ka-band
    link.tx_power_watts = 10.0;       // 10W for CubeSat
    link.tx_gain_dbi = 30.0;          // 1m dish

    // Convert to distance in meters
    const double distance_m = distance_ly * constants::LY;

    // Free space path loss (Friis equation)
    const double wavelength_m = constants::C / (link.frequency_ghz * 1e9);
    link.path_loss_db = 20.0 * std::log10(4.0 * M_PI * distance_m / wavelength_m);

    // Receiver parameters (Earth DSN)
    link.rx_gain_dbi = 73.0;          // 70m DSN antenna
    link.system_noise_k = 20.0;       // System temperature

    // Add gravitational lens gain if enabled
    if (use_lens) {
        link.tx_gain_dbi += 10.0 * std::log10(lens_magnification);
    }

    // Calculate received power
    const double rx_power_dbm = 10.0 * std::log10(link.tx_power_watts * 1000.0) +
                               link.tx_gain_dbi - link.path_loss_db + link.rx_gain_dbi;

    // Noise power
    const double noise_power_dbm = 10.0 * std::log10(constants::K_B * link.system_noise_k * 1000.0);

    // SNR and data rate (Shannon limit)
    const double snr_db = rx_power_dbm - noise_power_dbm;
    const double bandwidth_hz = 10e6;  // 10 MHz
    link.data_rate_bps = bandwidth_hz * std::log2(1.0 + std::pow(10.0, snr_db / 10.0));

    // BER for BPSK
    const double eb_n0 = std::pow(10.0, snr_db / 10.0) / (link.data_rate_bps / bandwidth_hz);
    link.bit_error_rate = 0.5 * std::erfc(std::sqrt(eb_n0));

    // Link margin
    link.link_margin_db = snr_db - 10.0;  // 10 dB required SNR

    return link;
}

std::array<uint8_t, 4096> InterstellarTransmitter::encode_message(
    const uint8_t* data,
    uint32_t length,
    ErrorCorrection ecc
) {
    std::array<uint8_t, 4096> encoded{};

    // Suppress unused parameter warning (simplified version uses fixed encoding)
    (void)ecc;

    // Add sync header (alternating 0xAA 0x55 pattern)
    encoded[0] = 0xAA;
    encoded[1] = 0x55;
    encoded[2] = 0xAA;
    encoded[3] = 0x55;

    // Message length
    encoded[4] = (length >> 8) & 0xFF;
    encoded[5] = length & 0xFF;

    // Copy message
    for (uint32_t i = 0; i < length && i < 2048; ++i) {
        encoded[6 + i] = data[i];
    }

    // Add Reed-Solomon parity (simplified)
    // Real implementation would use full RS(255,223) encoding
    uint32_t parity_start = 6 + length;
    for (uint32_t i = 0; i < 32; ++i) {
        uint8_t parity = 0;
        for (uint32_t j = 0; j < length; ++j) {
            parity ^= data[j];  // Simplified parity
        }
        encoded[parity_start + i] = parity;
    }

    return encoded;
}

} // namespace solarlens
