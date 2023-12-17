//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_deutron.h"
#include <algorithm>
#include <iomanip>
#include <queue>
#include <random>

double calculate_fraction_batch_factor(
        const std::vector<ParticleData> &protons,
        const std::vector<ParticleData> &neutrons,
        std::vector<ParticleData> &protons_fraction,
        std::vector<ParticleData> &neutrons_fraction,
        const config_in &config_input) {
    std::random_device rd;
    std::mt19937 g(rd());

    protons_fraction = protons;
    std::shuffle(protons_fraction.begin(), protons_fraction.end(), g);
    if (protons_fraction.size() > static_cast<size_t>(config_input.proton_limit)) {
        protons_fraction.resize(config_input.proton_limit);
    }

    neutrons_fraction = neutrons;
    std::shuffle(neutrons_fraction.begin(), neutrons_fraction.end(), g);
    if (neutrons_fraction.size() > static_cast<size_t>(config_input.neutron_limit)) {
        neutrons_fraction.resize(config_input.neutron_limit);
    }

    auto kp    = static_cast<double>(protons.size());
    auto kn    = static_cast<double>(neutrons.size());
    auto kpmix = static_cast<double>(protons_fraction.size());
    auto knmix = static_cast<double>(neutrons_fraction.size());

    double sfactor = (kp * (kp - 1) / (kpmix * (kpmix - 1)) * 2.0) *
                     (kn * (kn - 1) / (knmix * (knmix - 1)) * 2.0);
    //    double sfactor = (kp * kn) / (kpmix * knmix);

    return sfactor;
}
void update_momentum_array(double pt, double probability, const config_in &config_input,
                           std::vector<double> &d_mix_spv) {
    int npt = static_cast<int>(pt / config_input.d_mix_dpt);
    if (npt < d_mix_spv.size()) {
        d_mix_spv[npt] += probability;
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
calculate_differences(const ParticleData &p1, const ParticleData &p2) {
    double diff_x  = (p1.x - p2.x) / sqrt(2.0);
    double diff_y  = (p1.y - p2.y) / sqrt(2.0);
    double diff_z  = (p1.z - p2.z) / sqrt(2.0);
    double diff_px = (p1.px - p2.px) / sqrt(2.0);
    double diff_py = (p1.py - p2.py) / sqrt(2.0);
    double diff_pz = (p1.pz - p2.pz) / sqrt(2.0);
    return {diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz};
}

ParticleData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                    const config_in &config_input, std::vector<double> &d_mix_spv) {
    ParticleData deutron_data{};
    // Rapidity checks
    if (std::abs(calculate_rap(p1.p0, p1.pz)) > config_input.rap_cut_nucl ||
        std::abs(calculate_rap(p2.p0, p2.pz)) > config_input.rap_cut_nucl) {
        return ParticleData{};
    }

    const double px_total = p1.px + p2.px;
    const double py_total = p1.py + p2.py;
    const double pz_total = p1.pz + p2.pz;
    const double p0_total = p1.p0 + p2.p0;
    //Deutron rapidity check
    if (std::abs(calculate_rap(p0_total, pz_total)) > config_input.rap_cut_coal) {
        return ParticleData{};
    }

    const double beta_x = px_total / p0_total;
    const double beta_y = py_total / p0_total;
    const double beta_z = pz_total / p0_total;

    ParticleData boosted_proton  = p1.lorentz_boost(beta_x, beta_y, beta_z);
    ParticleData boosted_neutron = p2.lorentz_boost(beta_x, beta_y, beta_z);
    const double t_max           = std::max(boosted_proton.freeze_out_time, boosted_neutron.freeze_out_time);

    boosted_proton.update_position(t_max);
    boosted_neutron.update_position(t_max);
    const auto [diff_x, diff_y, diff_z,
                diff_px, diff_py, diff_pz] =
            calculate_differences(boosted_proton, boosted_neutron);
    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);

    if (diff_p > config_input.cut_dp * 0.19733 / config_input.rms || diff_r > config_input.cut_dr * config_input.rms) {
        return ParticleData{};
    }
    deutron_data.probability = calculate_pro_deutron(diff_r, diff_p, config_input.rms);

    deutron_data.calculate_deutron_data(boosted_proton, boosted_neutron);
    ParticleData boosted_deutron = deutron_data.lorentz_boost(-beta_x, -beta_y, -beta_z);
    const double pt              = sqrt(px_total * px_total + py_total * py_total);
    update_momentum_array(pt, deutron_data.probability, config_input, d_mix_spv);
    return boosted_deutron;
}

//void calculate_deutron_batch(const std::vector<ParticleData> &protons,
//                             const std::vector<ParticleData> &neutrons,
//                             const config_in &config_input, std::vector<double> &d_mix_spv,
//                             const std::vector<double> &d_mix_ptv, double &batch_deutrons,
//                             int eventsInBatch, std::vector<ParticleData> &deutrons) {
//    deutrons.clear();
//    batch_deutrons      = 0.0;
//    const int mixEvents = eventsInBatch * eventsInBatch;
//
//    auto compare = [](const ParticleData &a, const ParticleData &b) {
//        return a.probability < b.probability;
//    };
//    std::priority_queue<ParticleData, std::vector<ParticleData>, decltype(compare)> deuteron_queue(
//            compare);
//
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> dis(0.0, 1.0);
//
//    if (config_input.set_on) {
//
//        std::vector<ParticleData> protons_fraction, neutrons_fraction;
//        double factor = calculate_fraction_batch_factor(protons, neutrons, protons_fraction, neutrons_fraction, config_input);
//        std::cout << "Factor: " << factor << "proton size" << protons_fraction.size() << std::endl;
//        for (const auto &proton: protons_fraction) {
//            for (const auto &neutron: neutrons_fraction) {
//                ParticleData deutron = calculate_info_deutron(proton, neutron, config_input, d_mix_spv);
//                batch_deutrons += deutron.probability * factor;
//                if (deutron.probability > 0) {
//                    deuteron_queue.push(deutron);
//                }
//            }
//        }
//
//        for (int k = 0; k < d_mix_spv.size(); ++k) {
//            if (d_mix_ptv[k] > 0) {
//                d_mix_spv[k] = d_mix_spv[k] / (2 * M_PI * d_mix_ptv[k] * config_input.d_mix_dpt * mixEvents) * factor;
//            }
//        }
//    } else {
//        for (const auto &proton: protons) {
//            for (const auto &neutron: neutrons) {
//                ParticleData deutron = calculate_info_deutron(proton, neutron, config_input, d_mix_spv);
//                batch_deutrons += deutron.probability;
//                if (deutron.probability > 0) {
//                    deuteron_queue.push(deutron);
//                }
//            }
//        }
//
//        for (int k = 0; k < d_mix_spv.size(); ++k) {
//            if (d_mix_ptv[k] > 0) {
//                d_mix_spv[k] = d_mix_spv[k] / (2 * M_PI * d_mix_ptv[k] * config_input.d_mix_dpt * mixEvents);
//            }
//        }
//    }
//
//    int expected_deutrons = static_cast<int>(batch_deutrons);
//    double fraction       = batch_deutrons - expected_deutrons;
//    if (dis(gen) < fraction) {
//        expected_deutrons++;
//    }
//    while (!deuteron_queue.empty() && expected_deutrons > 0) {
//        deutrons.push_back(deuteron_queue.top());
//        deuteron_queue.pop();
//        expected_deutrons--;
//    }
//}

void calculate_deutron_batch(const std::vector<ParticleData> &protons,
                             const std::vector<ParticleData> &neutrons,
                             const config_in &config_input, std::vector<double> &d_mix_spv,
                             const std::vector<double> &d_mix_ptv, double &batch_deutrons,
                             int eventsInBatch, std::vector<ParticleData> &deutrons) {
    deutrons.clear();
    batch_deutrons      = 0.0;
    const int mixEvents = eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_deutrons;
    std::vector<double> cumulative_probabilities;

    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            ParticleData deutron = calculate_info_deutron(proton, neutron, config_input, d_mix_spv);
            if (deutron.probability > 0) {
                potential_deutrons.emplace_back(deutron, deutron.probability);
                batch_deutrons += deutron.probability;
                cumulative_probabilities.push_back(batch_deutrons);
            }
        }
    }

    weighted_sampling(deutrons, potential_deutrons, cumulative_probabilities, batch_deutrons);

    for (int k = 0; k < 100; ++k) {
        if (d_mix_ptv[k] > 0) {
            d_mix_spv[k] /= (2 * M_PI * d_mix_ptv[k] * config_input.d_mix_dpt * mixEvents);
        }
    }
}

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const std::string &deutronFile,
                        const config_in &configInput, int batchSize,
                        std::vector<double> &dMixSpv, const std::vector<double> &dMixPtv,
                        std::vector<ParticleData> &deutrons) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    read_batch_nuclei(protonFile, batchSize, protonBatches);
    read_batch_nuclei(neutronFile, batchSize, neutronBatches);

    std::ofstream output(deutronFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons       = batchData.particles;
        const auto &neutrons      = neutronBatches[batchNumber].particles;
        double num_batch_deutrons = 0.0;
        const int eventsInBatch   = batchData.eventCount;
        std::vector<ParticleData> batchDeutrons;
        calculate_deutron_batch(protons,
                                neutrons,
                                configInput,
                                dMixSpv,
                                dMixPtv,
                                num_batch_deutrons,
                                eventsInBatch, batchDeutrons);
        total_deutrons += num_batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
        deutrons.insert(deutrons.end(), batchDeutrons.begin(), batchDeutrons.end());
        if (output.is_open()) {
            output_deutrons(batchDeutrons, output);
        }
        std::cout << "Average number of deuteron per batch: " << num_batch_deutrons / eventsInBatch / eventsInBatch << std::endl;
    }
    output.close();
    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per event: " << average_deutrons << std::endl;
    output_d_mix_spv(dMixSpv, dMixPtv, "tem/d_mix_spv.dat", total_batches);
}
void weighted_sampling(std::vector<ParticleData> &sample_particles,
                       std::vector<std::pair<ParticleData, double>> &potential_particles,
                       std::vector<double> &cumulative_probabilities,
                       double number_of_particles) {
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
        auto it              = std::lower_bound(cumulative_probabilities.begin(),
                                                cumulative_probabilities.end(), random_number);

        if (it != cumulative_probabilities.end()) {
            size_t index = std::distance(cumulative_probabilities.begin(), it);
            sample_particles.push_back(potential_particles[index].first);
            cumulative_probabilities.erase(it);
            potential_particles.erase(potential_particles.begin() + static_cast<std::ptrdiff_t>(index));
        }
        expected_particles--;
    }
}
