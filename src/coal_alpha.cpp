//
// Created by mafu on 12/15/2023.
//
#include "../include/coal_alpha.h"
#include <algorithm>
#include <filesystem>
#include <queue>
#include <random>

ParticleData calculate_info_alpha_twobody(const ParticleData &d1, const ParticleData &d2,
                                          const config_in &config_input, std::vector<double> &alpha_mix_spv) {
    double rms   = config_input.alpha_rms;
    double hbar2 = 0.19733 * 0.19733;
    ParticleData alpha{};
    //Rapidity check
    double rapidity_deutron1 = calculate_rap(d1.p0, d1.pz);
    double rapidity_deutron2 = calculate_rap(d2.p0, d2.pz);
    if (std::abs(rapidity_deutron1) > config_input.rap_cut_nucl || std::abs(rapidity_deutron2) > config_input.rap_cut_nucl) {
        return ParticleData{};
    }
    //alpha rapidity check
    double px_total       = d1.px + d2.px;
    double py_total       = d1.py + d2.py;
    double pz_total       = d1.pz + d2.pz;
    double p0_total       = d1.p0 + d2.p0;
    double rapidity_alpha = calculate_rap(p0_total, pz_total);
    if (std::abs(rapidity_alpha) > config_input.rap_cut_coal) {
        return ParticleData{};
    }
    //lorentz boost calculate
    const double beta_x   = px_total / p0_total;
    const double beta_y   = py_total / p0_total;
    const double beta_z   = pz_total / p0_total;
    ParticleData boost_d1 = d1.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boost_d2 = d2.lorentz_boost(beta_x, beta_y, beta_z);
    //find t_max in d1 d2
    double t_max = std::max(boost_d1.freeze_out_time, boost_d2.freeze_out_time);
    //update position
    boost_d1.update_position(t_max);
    boost_d2.update_position(t_max);
    //Calculate the position difference and momentum difference between two particles
    auto [diff_dx, diff_dy, diff_dz,
          diff_dpx, diff_dpy, diff_dpz] = calculate_differences(boost_d1, boost_d2);

    double diff_dr = sqrt(diff_dx * diff_dx + diff_dy * diff_dy + diff_dz * diff_dz);
    double diff_dp = sqrt(diff_dpx * diff_dpx + diff_dpy * diff_dpy + diff_dpz * diff_dpz);

    if (diff_dr > config_input.cut_dr * rms || diff_dp > config_input.cut_dp * 0.19733 / rms) {
        return ParticleData{};
    }
    // gc = 2j+1/2^N
    // deutron : j = 1, N = 2
    // helium3 : j = 1/2, N = 3
    // alpha : j = 0, N = 4

    alpha.probability = 1.0 / 4 * 8 * exp(-diff_dr * diff_dr / rms / rms - diff_dp * diff_dp * rms * rms / hbar2);
    alpha.calculate_deutron_data(boost_d1, boost_d2);
    ParticleData alpha_boost = alpha.lorentz_boost(-beta_x, -beta_y, -beta_z);
    const double pt          = sqrt(px_total * px_total + py_total * py_total);
    update_momentum_array(pt, alpha_boost.probability, config_input, alpha_mix_spv);
    return alpha_boost;
}

ParticleData calculate_info_alpha_fourbody(const ParticleData &p1, const ParticleData &p2,
                                           const ParticleData &n1, const ParticleData &n2,
                                           const config_in &config_input, std::vector<double> &alpha_mix_spv) {
    double sig1                  = config_input.alpha_rms;
    double sig2                  = config_input.alpha_rms;
    double sig3                  = config_input.alpha_rms;
    constexpr const double hbar2 = 0.19733 * 0.19733;

    ParticleData alpha_particle{};
    //Rapidity check
    double rapidity_proton1  = calculate_rap(p1.p0, p1.pz);
    double rapidity_proton2  = calculate_rap(p2.p0, p2.pz);
    double rapidity_neutron1 = calculate_rap(n1.p0, n1.pz);
    double rapidity_neutron2 = calculate_rap(n2.p0, n2.pz);
    if (std::abs(rapidity_proton1) > config_input.rap_cut_nucl || std::abs(rapidity_proton2) > config_input.rap_cut_nucl ||
        std::abs(rapidity_neutron1) > config_input.rap_cut_nucl || std::abs(rapidity_neutron2) > config_input.rap_cut_nucl) {
        return ParticleData{};
    }
    //alpha rapidity check
    double px_total       = p1.px + p2.px + n1.px + n2.px;
    double py_total       = p1.py + p2.py + n1.py + n2.py;
    double pz_total       = p1.pz + p2.pz + n1.pz + n2.pz;
    double p0_total       = p1.p0 + p2.p0 + n1.p0 + n2.p0;
    double rapidity_alpha = calculate_rap(p0_total, pz_total);
    if (std::abs(rapidity_alpha) > config_input.rap_cut_coal) {
        return ParticleData{};
    }
    //lorentz boost calculate
    const double beta_x   = px_total / p0_total;
    const double beta_y   = py_total / p0_total;
    const double beta_z   = pz_total / p0_total;
    ParticleData boost_p1 = p1.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boost_p2 = p2.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boost_n1 = n1.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boost_n2 = n2.lorentz_boost(beta_x, beta_y, beta_z);
    //find t_max in p1 p2 n1 n2
    double t_max = std::max(std::max(boost_p1.freeze_out_time, boost_p2.freeze_out_time),
                            std::max(boost_n1.freeze_out_time, boost_n2.freeze_out_time));
    //update position
    boost_p1.update_position(t_max);
    boost_p2.update_position(t_max);
    boost_n1.update_position(t_max);
    boost_n2.update_position(t_max);
    //Calculate the position difference and momentum difference between two particles two by two
    auto [diff_dr_p1p2, diff_dp_p1p2, diff_dr_p1p2n1,
          diff_dp_p1p2n1, diff_dr_p1p2n1n2, diff_dp_p1p2n1n2] =
            calculate_diff_fourbody(boost_p1, boost_p2, boost_n1, boost_n2);

    if (diff_dr_p1p2 > config_input.cut_dr * sig1 || diff_dp_p1p2 > config_input.cut_dp * 0.19733 / sig1 ||
        diff_dr_p1p2n1 > config_input.cut_dr * sig2 || diff_dp_p1p2n1 > config_input.cut_dp * 0.19733 / sig2 ||
        diff_dr_p1p2n1n2 > config_input.cut_dr * sig3 || diff_dp_p1p2n1n2 > config_input.cut_dp * 0.19733 / sig3) {
        return ParticleData{};
    }

    double diff_dr = diff_dr_p1p2 * diff_dr_p1p2 / sig1 / sig1 + diff_dp_p1p2 * diff_dp_p1p2 * sig1 * sig1 / hbar2 +
                     diff_dr_p1p2n1 * diff_dr_p1p2n1 / sig2 / sig2 + diff_dp_p1p2n1 * diff_dp_p1p2n1 * sig2 * sig2 / hbar2 +
                     diff_dr_p1p2n1n2 * diff_dr_p1p2n1n2 / sig3 / sig3 + diff_dp_p1p2n1n2 * diff_dp_p1p2n1n2 * sig3 * sig3 / hbar2;
    if (diff_dr < 40) {
        alpha_particle.probability = 1.0 / 16 * 8 * 8 * 8 * exp(-diff_dr);
    } else {
        alpha_particle.probability = 0;
    }
    alpha_particle.calculate_fourbody_data(boost_p1, boost_p2, boost_n1, boost_n2);
    ParticleData alpha_particle_boost = alpha_particle.lorentz_boost(-beta_x, -beta_y, -beta_z);
    const double pt                   = sqrt(px_total * px_total + py_total * py_total);
    update_momentum_array(pt, alpha_particle_boost.probability, config_input, alpha_mix_spv);
    return alpha_particle_boost;
}
std::tuple<double, double, double, double, double, double>
calculate_diff_fourbody(const ParticleData &p1, const ParticleData &p2,
                        const ParticleData &n1, const ParticleData &n2) {
    double dx_p1p2      = p1.x - p2.x;
    double dy_p1p2      = p1.y - p2.y;
    double dz_p1p2      = p1.z - p2.z;
    double dpx_p1p2     = (p2.mass * p1.px - p1.mass * p2.px) / (p1.mass + p2.mass) * sqrt(2.0);
    double dpy_p1p2     = (p2.mass * p1.py - p1.mass * p2.py) / (p1.mass + p2.mass) * sqrt(2.0);
    double dpz_p1p2     = (p2.mass * p1.pz - p1.mass * p2.pz) / (p1.mass + p2.mass) * sqrt(2.0);
    double diff_dr_p1p2 = sqrt(dx_p1p2 * dx_p1p2 + dy_p1p2 * dy_p1p2 + dz_p1p2 * dz_p1p2);
    double diff_dp_p1p2 = sqrt(dpx_p1p2 * dpx_p1p2 + dpy_p1p2 * dpy_p1p2 + dpz_p1p2 * dpz_p1p2);

    double dx_p1p2n1 = (p1.mass * p1.x + p2.mass * p2.x - (p1.mass + p2.mass) * n1.x) /
                       (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dy_p1p2n1 = (p1.mass * p1.y + p2.mass * p2.y - (p1.mass + p2.mass) * n1.y) /
                       (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dz_p1p2n1 = (p1.mass * p1.z + p2.mass * p2.z - (p1.mass + p2.mass) * n1.z) /
                       (p1.mass + p2.mass) * sqrt(2.0 / 3.0);
    double dpx_p1p2n1 = (n1.mass * p1.px + n1.mass * p2.px - (p1.mass + p2.mass) * n1.px) /
                        (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double dpy_p1p2n1 = (n1.mass * p1.py + n1.mass * p2.py - (p1.mass + p2.mass) * n1.py) /
                        (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double dpz_p1p2n1 = (n1.mass * p1.pz + n1.mass * p2.pz - (p1.mass + p2.mass) * n1.pz) /
                        (p1.mass + p2.mass + n1.mass) * sqrt(6) / 2;
    double diff_dr_p1p2n1 = sqrt(dx_p1p2n1 * dx_p1p2n1 + dy_p1p2n1 * dy_p1p2n1 + dz_p1p2n1 * dz_p1p2n1);
    double diff_dp_p1p2n1 =
            sqrt(dpx_p1p2n1 * dpx_p1p2n1 + dpy_p1p2n1 * dpy_p1p2n1 + dpz_p1p2n1 * dpz_p1p2n1);

    double dx_p1p2n1n2 = (p1.mass * p1.x + p2.mass * p2.x +
                          n1.mass * n1.x - (p1.mass + p2.mass + n1.mass) * n2.x) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dy_p1p2n1n2 = (p1.mass * p1.y + p2.mass * p2.y +
                          n1.mass * n1.y - (p1.mass + p2.mass + n1.mass) * n2.y) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dz_p1p2n1n2 = (p1.mass * p1.z + p2.mass * p2.z +
                          n1.mass * n1.z - (p1.mass + p2.mass + n1.mass) * n2.z) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dpx_p1p2n1n2 = (n2.mass * (p1.px + p2.px + n1.px) -
                           (p1.mass + p2.mass + n1.mass) * n2.px) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double dpy_p1p2n1n2 = (n2.mass * (p1.py + p2.py + n1.py) -
                           (p1.mass + p2.mass + n1.mass) * n2.py) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double dpz_p1p2n1n2 = (n2.mass * (p1.pz + p2.pz + n1.pz) -
                           (p1.mass + p2.mass + n1.mass) * n2.pz) /
                          (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double diff_dr_p1p2n1n2 = sqrt(dx_p1p2n1n2 * dx_p1p2n1n2 + dy_p1p2n1n2 * dy_p1p2n1n2 +
                                   dz_p1p2n1n2 * dz_p1p2n1n2);
    double diff_dp_p1p2n1n2 = sqrt(dpx_p1p2n1n2 * dpx_p1p2n1n2 + dpy_p1p2n1n2 * dpy_p1p2n1n2 +
                                   dpz_p1p2n1n2 * dpz_p1p2n1n2);
    return std::make_tuple(diff_dr_p1p2, diff_dp_p1p2, diff_dr_p1p2n1, diff_dp_p1p2n1, diff_dr_p1p2n1n2, diff_dp_p1p2n1n2);
}

void calculate_alpha_fourbody_batch(const std::vector<ParticleData> &protons, const std::vector<ParticleData> &neutrons,
                                    const config_in &config_input, std::vector<double> &alpha_mix_spv,
                                    const std::vector<double> &alpha_mix_ptv, double &batch_alpha, int eventsInBatch,
                                    std::vector<ParticleData> &alpha) {
    alpha.clear();
    batch_alpha = 0.0;

    std::vector<ParticleData> protons_fraction, neutrons_fraction;
    double factor = calculate_fraction_batch_factor(protons, neutrons,
                                                    protons_fraction, neutrons_fraction, config_input);

    const int mixEvents = eventsInBatch * eventsInBatch * eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_alpha;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < protons_fraction.size(); i++) {
        for (size_t j = i + 1; j < protons_fraction.size(); j++) {
            for (size_t k = 0; k < neutrons_fraction.size(); k++) {
                for (size_t l = k + 1; l < neutrons_fraction.size(); l++) {
                    ParticleData alpha_particle = calculate_info_alpha_fourbody(
                            protons_fraction[i], protons_fraction[j],
                            neutrons_fraction[k], neutrons_fraction[l],
                            config_input, alpha_mix_spv);
                    if (alpha_particle.probability > 0) {
                        potential_alpha.emplace_back(alpha_particle, alpha_particle.probability);
                        batch_alpha += alpha_particle.probability * factor;
                        cumulated_probabilities.push_back(batch_alpha);
                    }
                }
            }
        }
    }
    weighted_sampling(alpha, potential_alpha, cumulated_probabilities, batch_alpha);

    for (int k = 0; k < alpha_mix_spv.size(); ++k) {
        if (alpha_mix_ptv[k] > 0) {
            alpha_mix_spv[k] = alpha_mix_spv[k] / (2 * M_PI * alpha_mix_ptv[k] * config_input.d_mix_dpt * mixEvents) * factor;
        }
    }
}

void calculate_alpha_fourbody(const std::string &protonFile, const std::string &neutronFile,
                              const std::string &alphaFile, const config_in &configInput, int batchSize,
                              std::vector<double> &alphaMixSpv, const std::vector<double> &alphaMixPtv,
                              std::vector<ParticleData> &alpha) {
    BatchMap protonBatches, neutronBatches;
    double total_alpha = 0.0;
    int total_batches  = 0;

    read_batch_nuclei(protonFile, batchSize, protonBatches);
    read_batch_nuclei(neutronFile, batchSize, neutronBatches);
    std::ofstream alphaFileOut(alphaFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons       = batchData.particles;
        const auto &neutrons      = neutronBatches[batchNumber].particles;
        int eventsInBatch         = batchData.eventCount;
        double batch_number_alpha = 0.0;
        std::vector<ParticleData> batch_alpha;
        calculate_alpha_fourbody_batch(protons, neutrons, configInput,
                                       alphaMixSpv, alphaMixPtv, batch_number_alpha,
                                       eventsInBatch, batch_alpha);
        total_alpha += batch_number_alpha / eventsInBatch / eventsInBatch / eventsInBatch / eventsInBatch;
        total_batches++;
        alpha.insert(alpha.end(), batch_alpha.begin(), batch_alpha.end());
        if (alphaFileOut.is_open()) {
            output_deutrons(batch_alpha, alphaFileOut);
        }
        std::cout << "average number of alpha per batch:" << batch_number_alpha / eventsInBatch / eventsInBatch / eventsInBatch / eventsInBatch
                  << std::endl;
    }
    alphaFileOut.close();
    const double average_alpha = total_alpha > 0 ? total_alpha / total_batches : 0.0;
    std::cout << "average number of alpha:" << average_alpha << std::endl;
    output_d_mix_spv(alphaMixSpv, alphaMixPtv, "tem/alpha_mix_spv_four.dat", total_batches);
}
void calculate_alpha_twobody_batch(const std::vector<ParticleData> &deutrons, const config_in &config_input,
                                   std::vector<double> &alpha_mix_spv, const std::vector<double> &alpha_mix_ptv,
                                   double &batch_alpha, int eventsInBatch, std::vector<ParticleData> &alpha) {
    alpha.clear();
    batch_alpha         = 0.0;
    const int mixEvents = eventsInBatch * eventsInBatch;

    std::vector<std::pair<ParticleData, double>> potential_alpha;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < deutrons.size(); i++) {
        for (size_t j = i + 1; j < deutrons.size(); j++) {
            ParticleData alpha_particle = calculate_info_alpha_twobody(
                    deutrons[i], deutrons[j],
                    config_input, alpha_mix_spv);
            if (alpha_particle.probability > 0) {
                potential_alpha.emplace_back(alpha_particle, alpha_particle.probability);
                batch_alpha += alpha_particle.probability;
                cumulated_probabilities.push_back(batch_alpha);
            }
        }
    }

    weighted_sampling(alpha, potential_alpha, cumulated_probabilities, batch_alpha);

    for (int k = 0; k < alpha_mix_spv.size(); ++k) {
        if (alpha_mix_ptv[k] > 0) {
            alpha_mix_spv[k] /= (2 * M_PI * alpha_mix_ptv[k] * config_input.d_mix_dpt * mixEvents);
        }
    }
}
void calculate_alpha_twobody(const std::string &deuteronFile, const std::string &alphaFile,
                             const config_in &configInput, int batchSize, std::vector<double> &alphaMixSpv,
                             const std::vector<double> &alphaMixPtv, std::vector<ParticleData> &alpha) {
    double total_alpha = 0.0;
    int total_batches  = 0;
    int mixEvents      = batchSize * batchSize;
    std::vector<BatchData> deuteronBatches;
    read_batch_deutrons(deuteronFile, deuteronBatches);

    std::ofstream alphaFileOut(alphaFile, std::ios::out);

    for (const auto &batchData: deuteronBatches) {
        const auto &deutrons = batchData.particles;
        int batchNumber      = batchData.eventCount;

        double batch_number_alpha = 0.0;
        std::vector<ParticleData> batch_alpha;
        calculate_alpha_twobody_batch(deutrons, configInput,
                                      alphaMixSpv, alphaMixPtv,
                                      batch_number_alpha, mixEvents, batch_alpha);

        total_alpha += batch_number_alpha / mixEvents / mixEvents;
        total_batches++;

        alpha.insert(alpha.end(), batch_alpha.begin(), batch_alpha.end());
        if (alphaFileOut.is_open()) {
            output_deutrons(batch_alpha, alphaFileOut);
        }
        std::cout << "Average number of alpha per batch (Batch " << batchNumber << "): " << batch_number_alpha / mixEvents / mixEvents << std::endl;
    }

    alphaFileOut.close();
    const double average_alpha = total_batches > 0 ? total_alpha / total_batches : 0.0;
    std::cout << "Average number of alpha: " << average_alpha << std::endl;

    output_d_mix_spv(alphaMixSpv, alphaMixPtv, "tem/alpha_mix_spv_two.dat", total_batches);
}
