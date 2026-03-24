#pragma once

#include "solarlens/physics/constants.hpp"
#include <array>
#include <vector>
#include <cstdint>
#include <cmath>

namespace solarlens {

// ============================================================================
// DEEP SPACE COMMUNICATION & RELAY NETWORK
//
// Multi-layer communication architecture:
// 1. Inter-satellite laser links (1 Gbps, within swarm)
// 2. Swarm-to-relay laser links (10 Mbps, to relay chain)
// 3. Relay chain back to inner solar system
// 4. Final relay to Earth DSN (Ka-band or optical)
//
// Includes real Reed-Solomon error correction, link budget calculations
// with Friis equation, and optional SGL-amplified communication.
// ============================================================================

class DeepSpaceComms {
public:
    // ---------- Link Budget ----------

    struct LinkBudget {
        // Input parameters
        double frequency_hz;
        double tx_power_w;
        double tx_antenna_gain_dbi;
        double rx_antenna_gain_dbi;
        double distance_m;
        double bandwidth_hz;
        double system_noise_temp_k;

        // Computed results
        double eirp_dbw;              // Effective isotropic radiated power
        double free_space_loss_db;    // Path loss (Friis)
        double atmospheric_loss_db;   // Additional losses
        double received_power_dbw;    // Power at receiver
        double noise_power_dbw;       // Noise floor
        double snr_db;                // Signal-to-noise ratio
        double eb_n0_db;              // Energy per bit to noise
        double max_data_rate_bps;     // Shannon capacity
        double achievable_rate_bps;   // Practical rate (with coding overhead)
        double ber;                   // Bit error rate
        double link_margin_db;        // Margin above minimum required Eb/N0
    };

    // Compute link budget between two points
    // If sgl_gain > 1, applies gravitational lens amplification
    static LinkBudget compute_link_budget(
        double frequency_hz,
        double tx_power_w,
        double tx_gain_dbi,
        double rx_gain_dbi,
        double distance_m,
        double bandwidth_hz,
        double system_noise_k,
        double sgl_gain = 1.0
    );

    // Preset: inter-satellite laser link within swarm
    static LinkBudget inter_satellite_link(double distance_m);

    // Preset: satellite to Earth via DSN
    static LinkBudget satellite_to_earth_link(double distance_au);

    // Preset: SGL-amplified link (using Sun as antenna)
    static LinkBudget sgl_amplified_link(double distance_ly, double wavelength_nm);

    // ---------- Reed-Solomon Error Correction ----------

    // GF(2^8) Reed-Solomon codec
    // RS(255, 223) — corrects up to 16 symbol errors per block
    static constexpr uint32_t RS_N = 255;       // Codeword length
    static constexpr uint32_t RS_K = 223;       // Message length
    static constexpr uint32_t RS_T = 16;        // Error correction capability
    static constexpr uint32_t RS_PARITY = 32;   // Parity symbols (2*T)

    struct RSCodec {
        uint8_t gf_exp[512];       // GF(2^8) exponential table
        uint8_t gf_log[256];       // GF(2^8) logarithm table
        uint8_t generator[RS_PARITY + 1];  // Generator polynomial

        RSCodec();

        // Encode: input RS_K bytes, output RS_N bytes (message + parity)
        void encode(const uint8_t* message, uint8_t* codeword) const;

        // Decode: input RS_N bytes, correct errors in-place
        // Returns number of errors corrected, or -1 if uncorrectable
        int decode(uint8_t* codeword) const;

    private:
        uint8_t gf_mul(uint8_t a, uint8_t b) const;
        uint8_t gf_div(uint8_t a, uint8_t b) const;
        uint8_t gf_pow(uint8_t a, int n) const;
        uint8_t gf_inv(uint8_t a) const;
        void gf_poly_mul(const uint8_t* p, uint32_t pn,
                         const uint8_t* q, uint32_t qn,
                         uint8_t* result) const;
    };

    // ---------- Message Framing ----------

    struct Frame {
        static constexpr uint32_t SYNC_SIZE = 4;
        static constexpr uint32_t HEADER_SIZE = 16;
        static constexpr uint32_t MAX_PAYLOAD = RS_K - HEADER_SIZE;
        static constexpr uint32_t MAX_FRAME = RS_N + SYNC_SIZE;

        uint8_t sync[SYNC_SIZE];      // Sync pattern
        uint8_t version;              // Protocol version
        uint8_t frame_type;           // 0=telemetry, 1=science, 2=command
        uint16_t sequence;            // Frame counter
        uint32_t timestamp;           // Mission time (seconds)
        uint16_t source_id;           // Satellite ID
        uint16_t dest_id;             // Destination (0xFFFF = broadcast)
        uint16_t payload_length;      // Payload bytes
        uint8_t payload[MAX_PAYLOAD]; // Actual data
    };

    // Encode a frame with RS error correction and sync pattern
    static std::vector<uint8_t> encode_frame(const Frame& frame, const RSCodec& codec);

    // Decode a received frame, correcting errors
    // Returns true if frame decoded successfully
    static bool decode_frame(const uint8_t* data, uint32_t length,
                            Frame& frame, const RSCodec& codec);

    // ---------- Relay Network ----------

    struct RelayNode {
        uint32_t id;
        double distance_au;          // Distance from Sun
        double longitude_rad;        // Heliocentric longitude
        double tx_power_w;
        double antenna_gain_dbi;
        double data_buffer_bytes;    // On-board storage
    };

    // Plan relay chain from swarm at SGL focal distance to Earth
    // Returns optimal relay node positions and link budgets
    struct RelayChain {
        std::vector<RelayNode> nodes;
        std::vector<LinkBudget> links;
        double end_to_end_rate_bps;
        double end_to_end_latency_s;
        double total_power_w;
    };

    static RelayChain plan_relay_chain(
        double swarm_distance_au,
        uint32_t num_relays,
        double relay_tx_power_w
    );
};

} // namespace solarlens
