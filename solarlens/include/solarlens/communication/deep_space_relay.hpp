#pragma once

#include "solarlens/physics/constants.hpp"
#include <array>
#include <cstdint>
#include <cmath>

namespace solarlens {

// ============================================================================
// INTERSTELLAR COMMUNICATION
// ============================================================================
class InterstellarTransmitter {
public:
    // Error correction codes
    enum class ErrorCorrection : uint8_t {
        REED_SOLOMON,   // Classic deep space
        TURBO_CODES,    // Modern, high performance
        LDPC,           // Low density parity check
        POLAR_CODES     // Quantum-ready
    };

private:
    // Modulation schemes for deep space
    enum class Modulation : uint8_t {
        BPSK,      // Binary phase shift keying (robust)
        QPSK,      // Quadrature PSK (2x data rate)
        QAM16,     // 16-QAM (4x data rate)
        CHIRP,     // Frequency chirp (doppler resistant)
        MFSK       // Multi-frequency shift keying
    };

public:
    struct LinkBudget {
        double frequency_ghz;
        double tx_power_watts;
        double tx_gain_dbi;
        double path_loss_db;
        double rx_gain_dbi;
        double system_noise_k;
        double data_rate_bps;
        double bit_error_rate;
        double link_margin_db;
    };

    // Calculate link budget with gravitational lens gain
    [[nodiscard]] LinkBudget calculate_link_budget(
        double distance_ly,
        double lens_magnification,
        bool use_lens
    );

    // Encode message for interstellar transmission
    [[nodiscard]] std::array<uint8_t, 4096> encode_message(
        const uint8_t* data,
        uint32_t length,
        ErrorCorrection ecc
    );
};

} // namespace solarlens
