#include "solarlens/communication/deep_space_relay.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>

namespace solarlens {

// ============================================================================
// Link Budget Calculations
// ============================================================================

DeepSpaceComms::LinkBudget DeepSpaceComms::compute_link_budget(
    double frequency_hz,
    double tx_power_w,
    double tx_gain_dbi,
    double rx_gain_dbi,
    double distance_m,
    double bandwidth_hz,
    double system_noise_k,
    double sgl_gain
) {
    LinkBudget link{};
    link.frequency_hz = frequency_hz;
    link.tx_power_w = tx_power_w;
    link.tx_antenna_gain_dbi = tx_gain_dbi;
    link.rx_antenna_gain_dbi = rx_gain_dbi;
    link.distance_m = distance_m;
    link.bandwidth_hz = bandwidth_hz;
    link.system_noise_temp_k = system_noise_k;

    // EIRP (dBW)
    link.eirp_dbw = 10.0 * std::log10(tx_power_w) + tx_gain_dbi;

    // Free space path loss: L = (4πd/λ)²
    double wavelength_m = constants::C / frequency_hz;
    if (distance_m > 0 && wavelength_m > 0) {
        link.free_space_loss_db = 20.0 * std::log10(4.0 * constants::PI * distance_m / wavelength_m);
    }

    // Atmospheric/implementation losses
    link.atmospheric_loss_db = 1.0;  // 1 dB margin for pointing, polarization, etc.

    // Add SGL gain to receiver side
    double sgl_gain_db = (sgl_gain > 1.0) ? 10.0 * std::log10(sgl_gain) : 0.0;

    // Received power (dBW)
    link.received_power_dbw = link.eirp_dbw - link.free_space_loss_db -
                              link.atmospheric_loss_db + rx_gain_dbi + sgl_gain_db;

    // Noise power: N = kTB (dBW)
    link.noise_power_dbw = 10.0 * std::log10(constants::K_B * system_noise_k * bandwidth_hz);

    // SNR
    link.snr_db = link.received_power_dbw - link.noise_power_dbw;

    // Shannon capacity: C = B × log2(1 + SNR)
    double snr_linear = std::pow(10.0, link.snr_db / 10.0);
    link.max_data_rate_bps = bandwidth_hz * std::log2(1.0 + snr_linear);

    // Practical rate with coding overhead (~50% of Shannon for RS + turbo)
    link.achievable_rate_bps = link.max_data_rate_bps * 0.5;

    // Eb/N0 for achieved data rate
    if (link.achievable_rate_bps > 0) {
        link.eb_n0_db = link.snr_db - 10.0 * std::log10(link.achievable_rate_bps / bandwidth_hz);
    }

    // BER for BPSK with Eb/N0
    double eb_n0_linear = std::pow(10.0, link.eb_n0_db / 10.0);
    link.ber = 0.5 * std::erfc(std::sqrt(eb_n0_linear));

    // Link margin (relative to required Eb/N0 = 2.5 dB for RS-coded BPSK at BER=10^-6)
    link.link_margin_db = link.eb_n0_db - 2.5;

    return link;
}

DeepSpaceComms::LinkBudget DeepSpaceComms::inter_satellite_link(double distance_m) {
    // 1550nm laser link (telecom wavelength)
    return compute_link_budget(
        constants::C / 1.55e-6,  // ~193 THz (1550nm laser)
        constants::LASER_LINK_POWER_W,
        90.0,          // 10cm laser aperture ~ 90 dBi at 1550nm
        90.0,          // Same receiver aperture
        distance_m,
        1.0e9,         // 1 GHz bandwidth
        300.0          // ~300K noise temperature for photodetector
    );
}

DeepSpaceComms::LinkBudget DeepSpaceComms::satellite_to_earth_link(double distance_au) {
    // Ka-band (32 GHz) to DSN 70m antenna
    return compute_link_budget(
        32.0e9,        // 32 GHz Ka-band
        10.0,          // 10W transmitter
        40.0,          // 0.5m dish at Ka-band
        73.0,          // DSN 70m antenna gain
        distance_au * constants::AU,
        10.0e6,        // 10 MHz bandwidth
        20.0           // 20K system noise (cooled receiver)
    );
}

DeepSpaceComms::LinkBudget DeepSpaceComms::sgl_amplified_link(double distance_ly, double wavelength_nm) {
    // Using the SGL to communicate with a distant star system
    double lambda = wavelength_nm * 1e-9;
    double sgl_gain = 4.0 * constants::PI * constants::PI * constants::R_SCHWARZSCHILD / lambda;

    return compute_link_budget(
        constants::C / lambda,
        10.0,          // 10W laser
        90.0,          // 10cm aperture
        90.0,          // Receiver at other end
        distance_ly * constants::LY,
        1.0e9,
        300.0,
        sgl_gain
    );
}

// ============================================================================
// Reed-Solomon GF(2^8) Codec
// ============================================================================

DeepSpaceComms::RSCodec::RSCodec() {
    // Initialize GF(2^8) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
    uint32_t x = 1;
    for (uint32_t i = 0; i < 255; ++i) {
        gf_exp[i] = static_cast<uint8_t>(x);
        gf_log[x] = static_cast<uint8_t>(i);
        x <<= 1;
        if (x & 0x100) x ^= 0x11D;
    }
    // Extend exp table for convenience
    for (uint32_t i = 255; i < 512; ++i) {
        gf_exp[i] = gf_exp[i - 255];
    }
    gf_log[0] = 0;  // Undefined, but needed for safety

    // Build generator polynomial: g(x) = Π(x - α^i) for i = 0..2T-1
    // Start with g(x) = 1
    std::memset(generator, 0, sizeof(generator));
    generator[0] = 1;
    uint32_t gen_len = 1;

    for (uint32_t i = 0; i < RS_PARITY; ++i) {
        // Multiply g(x) by (x - α^i)
        uint8_t factor[2] = {1, gf_exp[i]};
        uint8_t result[RS_PARITY + 1] = {};
        gf_poly_mul(generator, gen_len, factor, 2, result);
        gen_len++;
        std::memcpy(generator, result, gen_len);
    }
}

uint8_t DeepSpaceComms::RSCodec::gf_mul(uint8_t a, uint8_t b) const {
    if (a == 0 || b == 0) return 0;
    return gf_exp[gf_log[a] + gf_log[b]];
}

uint8_t DeepSpaceComms::RSCodec::gf_div(uint8_t a, uint8_t b) const {
    if (b == 0) return 0;  // Division by zero
    if (a == 0) return 0;
    return gf_exp[(gf_log[a] + 255 - gf_log[b]) % 255];
}

uint8_t DeepSpaceComms::RSCodec::gf_pow(uint8_t a, int n) const {
    if (n == 0) return 1;
    if (a == 0) return 0;
    return gf_exp[(gf_log[a] * n) % 255];
}

uint8_t DeepSpaceComms::RSCodec::gf_inv(uint8_t a) const {
    if (a == 0) return 0;
    return gf_exp[255 - gf_log[a]];
}

void DeepSpaceComms::RSCodec::gf_poly_mul(
    const uint8_t* p, uint32_t pn,
    const uint8_t* q, uint32_t qn,
    uint8_t* result
) const {
    for (uint32_t i = 0; i < pn; ++i) {
        for (uint32_t j = 0; j < qn; ++j) {
            result[i + j] ^= gf_mul(p[i], q[j]);
        }
    }
}

void DeepSpaceComms::RSCodec::encode(const uint8_t* message, uint8_t* codeword) const {
    // Systematic RS encoding:
    // 1. Copy message to first K positions of codeword
    // 2. Compute remainder of message polynomial / generator polynomial
    // 3. Append remainder as parity symbols

    std::memcpy(codeword, message, RS_K);
    std::memset(codeword + RS_K, 0, RS_PARITY);

    // Polynomial long division
    uint8_t feedback[RS_PARITY + 1];
    std::memset(feedback, 0, sizeof(feedback));

    for (uint32_t i = 0; i < RS_K; ++i) {
        uint8_t coef = codeword[i] ^ feedback[0];

        // Shift feedback register
        std::memmove(feedback, feedback + 1, RS_PARITY);
        feedback[RS_PARITY] = 0;

        if (coef != 0) {
            for (uint32_t j = 0; j < RS_PARITY; ++j) {
                feedback[j] ^= gf_mul(generator[j + 1], coef);
            }
        }
    }

    // Copy parity
    for (uint32_t i = 0; i < RS_PARITY; ++i) {
        codeword[RS_K + i] = feedback[i];
    }
}

int DeepSpaceComms::RSCodec::decode(uint8_t* codeword) const {
    // Step 1: Compute syndromes S_i = C(α^i) for i = 0..2T-1
    uint8_t syndromes[RS_PARITY];
    bool has_errors = false;

    for (uint32_t i = 0; i < RS_PARITY; ++i) {
        uint8_t s = 0;
        for (uint32_t j = 0; j < RS_N; ++j) {
            s = gf_mul(s, gf_exp[i]) ^ codeword[j];
        }
        syndromes[i] = s;
        if (s != 0) has_errors = true;
    }

    if (!has_errors) return 0;  // No errors

    // Step 2: Berlekamp-Massey algorithm to find error locator polynomial
    uint8_t sigma[RS_PARITY + 1] = {};  // Error locator
    uint8_t old_sigma[RS_PARITY + 1] = {};
    sigma[0] = 1;
    old_sigma[0] = 1;
    uint32_t sigma_len = 1;
    uint32_t old_len = 1;

    for (uint32_t n = 0; n < RS_PARITY; ++n) {
        // Compute discrepancy
        uint8_t delta = syndromes[n];
        for (uint32_t j = 1; j < sigma_len; ++j) {
            delta ^= gf_mul(sigma[j], syndromes[n - j]);
        }

        // Shift old_sigma
        uint8_t shifted[RS_PARITY + 1] = {};
        for (uint32_t j = 0; j < old_len; ++j) {
            shifted[j + 1] = old_sigma[j];
        }

        if (delta != 0) {
            uint8_t new_sigma[RS_PARITY + 1];
            uint32_t new_len = std::max(sigma_len, old_len + 1);
            for (uint32_t j = 0; j < new_len; ++j) {
                new_sigma[j] = sigma[j] ^ gf_mul(delta, shifted[j]);
            }

            if (2 * (sigma_len - 1) <= n) {
                uint8_t inv_delta = gf_inv(delta);
                for (uint32_t j = 0; j < sigma_len; ++j) {
                    old_sigma[j] = gf_mul(sigma[j], inv_delta);
                }
                old_len = sigma_len;
                sigma_len = new_len;
            }

            std::memcpy(sigma, new_sigma, new_len);
        } else {
            // Just shift old
        }

        // Update old_sigma for next iteration if not updated above
        uint8_t temp[RS_PARITY + 1];
        std::memcpy(temp, old_sigma, sizeof(old_sigma));
        (void)temp;
    }

    // Step 3: Chien search — find error locations
    uint32_t num_errors = sigma_len - 1;
    if (num_errors > RS_T) return -1;  // Too many errors

    uint8_t error_locs[RS_T];
    uint32_t found = 0;

    for (uint32_t i = 0; i < RS_N; ++i) {
        uint8_t sum = 0;
        for (uint32_t j = 0; j < sigma_len; ++j) {
            sum ^= gf_mul(sigma[j], gf_pow(gf_exp[i], j));
        }
        if (sum == 0) {
            // α^i is a root → error at position N-1-i
            uint32_t pos = RS_N - 1 - i;
            if (pos < RS_N) {
                error_locs[found++] = static_cast<uint8_t>(pos);
            }
        }
    }

    if (found != num_errors) return -1;  // Couldn't find all roots

    // Step 4: Forney algorithm — compute error magnitudes
    for (uint32_t k = 0; k < found; ++k) {
        uint8_t Xi = gf_exp[RS_N - 1 - error_locs[k]];  // Error locator value
        uint8_t Xi_inv = gf_inv(Xi);

        // Evaluate syndrome polynomial at Xi_inv
        uint8_t omega = 0;
        uint8_t power = 1;
        for (uint32_t j = 0; j < RS_PARITY; ++j) {
            omega ^= gf_mul(syndromes[j], power);
            power = gf_mul(power, Xi_inv);
        }

        // Evaluate derivative of sigma at Xi_inv
        uint8_t sigma_prime = 0;
        for (uint32_t j = 1; j < sigma_len; j += 2) {
            sigma_prime ^= gf_mul(sigma[j], gf_pow(Xi_inv, j - 1));
        }

        if (sigma_prime == 0) return -1;  // Degenerate

        // Error magnitude
        uint8_t magnitude = gf_mul(omega, gf_inv(sigma_prime));
        codeword[error_locs[k]] ^= magnitude;
    }

    return static_cast<int>(found);
}

// ============================================================================
// Frame Encoding/Decoding
// ============================================================================

std::vector<uint8_t> DeepSpaceComms::encode_frame(const Frame& frame, const RSCodec& codec) {
    // Pack frame into message bytes
    uint8_t message[RS_K] = {};
    uint32_t offset = 0;

    message[offset++] = frame.version;
    message[offset++] = frame.frame_type;
    message[offset++] = (frame.sequence >> 8) & 0xFF;
    message[offset++] = frame.sequence & 0xFF;
    message[offset++] = (frame.timestamp >> 24) & 0xFF;
    message[offset++] = (frame.timestamp >> 16) & 0xFF;
    message[offset++] = (frame.timestamp >> 8) & 0xFF;
    message[offset++] = frame.timestamp & 0xFF;
    message[offset++] = (frame.source_id >> 8) & 0xFF;
    message[offset++] = frame.source_id & 0xFF;
    message[offset++] = (frame.dest_id >> 8) & 0xFF;
    message[offset++] = frame.dest_id & 0xFF;
    message[offset++] = (frame.payload_length >> 8) & 0xFF;
    message[offset++] = frame.payload_length & 0xFF;

    // Pad header to HEADER_SIZE
    while (offset < Frame::HEADER_SIZE) message[offset++] = 0;

    // Copy payload
    uint32_t copy_len = std::min(static_cast<uint32_t>(frame.payload_length),
                                  Frame::MAX_PAYLOAD);
    std::memcpy(message + Frame::HEADER_SIZE, frame.payload, copy_len);

    // RS encode
    uint8_t codeword[RS_N];
    codec.encode(message, codeword);

    // Build output: sync + codeword
    std::vector<uint8_t> output(Frame::SYNC_SIZE + RS_N);
    output[0] = 0xEB;  // CCSDS-like sync marker
    output[1] = 0x90;
    output[2] = 0xEB;
    output[3] = 0x90;
    std::memcpy(output.data() + Frame::SYNC_SIZE, codeword, RS_N);

    return output;
}

bool DeepSpaceComms::decode_frame(
    const uint8_t* data, uint32_t length,
    Frame& frame, const RSCodec& codec
) {
    if (length < Frame::SYNC_SIZE + RS_N) return false;

    // Verify sync pattern
    if (data[0] != 0xEB || data[1] != 0x90 || data[2] != 0xEB || data[3] != 0x90) {
        return false;
    }

    // RS decode
    uint8_t codeword[RS_N];
    std::memcpy(codeword, data + Frame::SYNC_SIZE, RS_N);
    int errors = codec.decode(codeword);
    if (errors < 0) return false;  // Uncorrectable

    // Unpack frame header
    uint32_t offset = 0;
    frame.sync[0] = data[0]; frame.sync[1] = data[1];
    frame.sync[2] = data[2]; frame.sync[3] = data[3];
    frame.version = codeword[offset++];
    frame.frame_type = codeword[offset++];
    frame.sequence = (static_cast<uint16_t>(codeword[offset]) << 8) | codeword[offset + 1]; offset += 2;
    frame.timestamp = (static_cast<uint32_t>(codeword[offset]) << 24) |
                      (static_cast<uint32_t>(codeword[offset + 1]) << 16) |
                      (static_cast<uint32_t>(codeword[offset + 2]) << 8) |
                      codeword[offset + 3];
    offset += 4;
    frame.source_id = (static_cast<uint16_t>(codeword[offset]) << 8) | codeword[offset + 1]; offset += 2;
    frame.dest_id = (static_cast<uint16_t>(codeword[offset]) << 8) | codeword[offset + 1]; offset += 2;
    frame.payload_length = (static_cast<uint16_t>(codeword[offset]) << 8) | codeword[offset + 1]; offset += 2;

    offset = Frame::HEADER_SIZE;
    uint32_t copy_len = std::min(static_cast<uint32_t>(frame.payload_length), Frame::MAX_PAYLOAD);
    std::memcpy(frame.payload, codeword + offset, copy_len);

    return true;
}

// ============================================================================
// Relay Network Planning
// ============================================================================

DeepSpaceComms::RelayChain DeepSpaceComms::plan_relay_chain(
    double swarm_distance_au,
    uint32_t num_relays,
    double relay_tx_power_w
) {
    RelayChain chain;

    // Place relays at geometrically-spaced distances between Earth (1 AU) and swarm
    // Geometric spacing optimizes the worst-case link (equalizes path losses)
    double ratio = std::pow(swarm_distance_au, 1.0 / (num_relays + 1));

    // First node is the swarm
    RelayNode swarm_node;
    swarm_node.id = 0;
    swarm_node.distance_au = swarm_distance_au;
    swarm_node.longitude_rad = 0;
    swarm_node.tx_power_w = relay_tx_power_w;
    swarm_node.antenna_gain_dbi = 40.0;
    swarm_node.data_buffer_bytes = 1e12;  // 1 TB
    chain.nodes.push_back(swarm_node);

    // Relay nodes
    for (uint32_t i = 0; i < num_relays; ++i) {
        RelayNode node;
        node.id = i + 1;
        node.distance_au = swarm_distance_au / std::pow(ratio, i + 1);
        node.longitude_rad = 0;  // Along radial line (simplified)
        node.tx_power_w = relay_tx_power_w;
        node.antenna_gain_dbi = 45.0;  // Larger relay antennas
        node.data_buffer_bytes = 1e12;
        chain.nodes.push_back(node);
    }

    // Earth station
    RelayNode earth;
    earth.id = num_relays + 1;
    earth.distance_au = 1.0;
    earth.longitude_rad = 0;
    earth.tx_power_w = 0;  // Receiver only
    earth.antenna_gain_dbi = 73.0;  // DSN 70m
    earth.data_buffer_bytes = 0;
    chain.nodes.push_back(earth);

    // Compute link budgets between consecutive nodes
    chain.end_to_end_rate_bps = 1e15;  // Will take minimum
    chain.end_to_end_latency_s = 0;
    chain.total_power_w = 0;

    for (uint32_t i = 0; i < chain.nodes.size() - 1; ++i) {
        double dist_au = std::abs(chain.nodes[i].distance_au - chain.nodes[i + 1].distance_au);
        double dist_m = dist_au * constants::AU;

        auto link = compute_link_budget(
            32.0e9,  // Ka-band
            chain.nodes[i].tx_power_w,
            chain.nodes[i].antenna_gain_dbi,
            chain.nodes[i + 1].antenna_gain_dbi,
            dist_m,
            10.0e6,  // 10 MHz bandwidth
            25.0     // System noise
        );

        chain.links.push_back(link);

        // End-to-end rate is bottleneck (minimum link)
        if (link.achievable_rate_bps < chain.end_to_end_rate_bps) {
            chain.end_to_end_rate_bps = link.achievable_rate_bps;
        }

        // Latency: light travel time
        chain.end_to_end_latency_s += dist_m / constants::C;
        chain.total_power_w += chain.nodes[i].tx_power_w;
    }

    return chain;
}

} // namespace solarlens
