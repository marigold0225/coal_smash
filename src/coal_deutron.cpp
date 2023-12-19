//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_deutron.h"
#include <algorithm>
#include <iomanip>
#include <queue>
#include <random>

void update_momentum_array(double pt, double probability, const config_in &config_input,
                           std::vector<double> &pt_array) {
    int npt = static_cast<int>(pt / config_input.d_mix_dpt);
    if (npt < pt_array.size()) {
        pt_array[npt] += probability;
    }
}
double calculate_pro_deutron(const double diff_r, const double diff_p, const double rms) {
    const double diff_r2 = diff_r * diff_r;
    const double diff_p2 = diff_p * diff_p;
    const double rms_2   = rms * rms;
    const double probability =
            3. / 4. * 8 * std::exp(-diff_r2 / rms_2 - diff_p2 * rms_2 / 0.19733 / 0.19733);
    return probability;
}
double calculate_rap(const double p0, const double pz) {
    if (p0 == pz) {
        return 0;
    }
    if (p0 == -pz) {
        return 1e10;
    }
    if (p0 < pz) {
        std::cout << "Error: p0 < pz" << std::endl;
        return 0;
    }
    const double rap = 0.5 * log((p0 + pz) / (p0 - pz));
    return rap;
}

std::tuple<double, double, double, double, double, double>
calculateTwoBodyDifferences(const ParticleData &p1, const ParticleData &p2) {
    double diff_x  = (p1.x - p2.x) / sqrt(2.0);
    double diff_y  = (p1.y - p2.y) / sqrt(2.0);
    double diff_z  = (p1.z - p2.z) / sqrt(2.0);
    double diff_px = (p1.px - p2.px) / sqrt(2.0);
    double diff_py = (p1.py - p2.py) / sqrt(2.0);
    double diff_pz = (p1.pz - p2.pz) / sqrt(2.0);
    return {diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz};
}

ParticleData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                    const config_in &config_input, std::vector<double> &pt_array) {
    ParticleData deutron_data{};
    double rap_nucl = config_input.d_rap_cut_nucl;
    double rap_coal = config_input.d_rap_cut_coal;
    // Rapidity checks
    if (std::abs(p1.get_rapidity()) > rap_nucl || std::abs(p2.get_rapidity()) > rap_nucl) {
        return ParticleData{};
    }

    const double px_total = p1.px + p2.px;
    const double py_total = p1.py + p2.py;
    const double pz_total = p1.pz + p2.pz;
    const double p0_total = p1.p0 + p2.p0;
    //Deutron rapidity check
    if (std::abs(calculate_rap(p0_total, pz_total)) > rap_coal) {
        return ParticleData{};
    }

    const double beta_x = px_total / p0_total;
    const double beta_y = py_total / p0_total;
    const double beta_z = pz_total / p0_total;

    ParticleData boosted_proton  = p1.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boosted_neutron = p2.lorentz_boost(beta_x, beta_y, beta_z);
    const double t_max = std::max(boosted_proton.freeze_out_time, boosted_neutron.freeze_out_time);

    boosted_proton.update_position(t_max);
    boosted_neutron.update_position(t_max);
    const auto [diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz] =
            calculateTwoBodyDifferences(boosted_proton, boosted_neutron);
    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);

    if (diff_p > config_input.cut_dp * 0.19733 / config_input.rms ||
        diff_r > config_input.cut_dr * config_input.rms) {
        return ParticleData{};
    }
    deutron_data.probability = calculate_pro_deutron(diff_r, diff_p, config_input.rms);

    deutron_data.get_twobody_data(boosted_proton, boosted_neutron);
    ParticleData boosted_deutron = deutron_data.lorentz_boost(-beta_x, -beta_y, -beta_z);
    const double pt              = sqrt(px_total * px_total + py_total * py_total);
    update_momentum_array(pt, deutron_data.probability, config_input, pt_array);
    return boosted_deutron;
}

void calculate_deutron_batch(const std::vector<ParticleData> &protons,
                             const std::vector<ParticleData> &neutrons,
                             const config_in &config_input, std::vector<double> &pt_array,
                             const std::vector<double> &d_pt, double &batch_deutrons,
                             int eventsInBatch, std::vector<ParticleData> &deutrons) {
    deutrons.clear();
    batch_deutrons      = 0.0;
    const int mixEvents = eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_deutrons;
    std::vector<double> cumulative_probabilities;

    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            ParticleData deutron = calculate_info_deutron(proton, neutron, config_input, pt_array);
            if (deutron.probability > 0) {
                potential_deutrons.emplace_back(deutron, deutron.probability);
                batch_deutrons += deutron.probability;
                cumulative_probabilities.push_back(batch_deutrons);
            }
        }
    }

    weighted_sampling(deutrons, potential_deutrons, cumulative_probabilities, batch_deutrons);

    for (int k = 0; k < 100; ++k) {
        if (d_pt[k] > 0) {
            pt_array[k] /= (2 * M_PI * d_pt[k] * config_input.d_mix_dpt * mixEvents);
        }
    }
}

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const std::string &deutronFile, const std::string &ptFile,
                        const config_in &configInput, std::vector<double> &pt_array) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    int batchSize         = configInput.mix_events;
    read_batch_nuclei(protonFile, batchSize, protonBatches);
    read_batch_nuclei(neutronFile, batchSize, neutronBatches);

    std::vector<double> deutron_d_pt(100);
    for (int i = 0; i < 100; ++i) {
        deutron_d_pt[i] =
                configInput.d_mix_dpt / 2 + static_cast<double>(i) * configInput.d_mix_dpt;
    }

    std::ofstream output(deutronFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons       = batchData.particles;
        const auto &neutrons      = neutronBatches[batchNumber].particles;
        double num_batch_deutrons = 0.0;
        const int eventsInBatch   = batchData.eventCount;
        std::vector<ParticleData> batchDeutrons;
        calculate_deutron_batch(protons, neutrons, configInput, pt_array, deutron_d_pt,
                                num_batch_deutrons, eventsInBatch, batchDeutrons);
        total_deutrons += num_batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
        if (output.is_open()) {
            output_cluster(batchDeutrons, output);
        }
        std::cout << "Average number of deuteron per batch: "
                  << num_batch_deutrons / eventsInBatch / eventsInBatch << std::endl;
    }
    output.close();
    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per event: " << average_deutrons << std::endl;
    output_spv(pt_array, deutron_d_pt, ptFile, total_batches);
}
void weighted_sampling(std::vector<ParticleData> &sample_particles,
                       std::vector<std::pair<ParticleData, double>> &potential_particles,
                       std::vector<double> &cumulative_probabilities, double number_of_particles) {
    if (potential_particles.empty() || cumulative_probabilities.empty()) {
        return;
    }

    int expected_particles = static_cast<int>(number_of_particles);
    double fraction        = number_of_particles - expected_particles;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (dis(gen) < fraction) {
        expected_particles++;
    }

    while (expected_particles > 0 && !cumulative_probabilities.empty()) {
        double upper_bound   = cumulative_probabilities.back();
        double random_number = dis(gen) * upper_bound;
        auto it = std::lower_bound(cumulative_probabilities.begin(), cumulative_probabilities.end(),
                                   random_number);

        if (it != cumulative_probabilities.end()) {
            size_t index = std::distance(cumulative_probabilities.begin(), it);
            sample_particles.push_back(potential_particles[index].first);
            cumulative_probabilities.erase(it);
            potential_particles.erase(potential_particles.begin() +
                                      static_cast<std::ptrdiff_t>(index));
        }
        expected_particles--;
    }
}
