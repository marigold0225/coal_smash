//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_deutron.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <queue>
#include <random>
ParticleData lorentz_boost(const double beta_x, const double beta_y, const double beta_z,
                           const ParticleData &p) {
    ParticleData boost_p = p;
    if (const double beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z; beta2 > 1.0e-5) {
        const double gamma   = 1.0 / sqrt(1.0 - beta2);
        const double x_lam00 = gamma;
        const double x_lam01 = -gamma * beta_x;
        const double x_lam02 = -gamma * beta_y;
        const double x_lam03 = -gamma * beta_z;
        const double x_lam11 = 1.0 + (gamma - 1.0) * beta_x * beta_x / beta2;
        const double x_lam22 = 1.0 + (gamma - 1.0) * beta_y * beta_y / beta2;
        const double x_lam33 = 1.0 + (gamma - 1.0) * beta_z * beta_z / beta2;
        const double x_lam12 = (gamma - 1.0) * beta_x * beta_y / beta2;
        const double x_lam13 = (gamma - 1.0) * beta_x * beta_z / beta2;
        const double x_lam23 = (gamma - 1.0) * beta_y * beta_z / beta2;
        const double new_t =
                p.freeze_out_time * x_lam00 + p.x * x_lam01 + p.y * x_lam02 + p.z * x_lam03;
        const double new_x =
                p.freeze_out_time * x_lam01 + p.x * x_lam11 + p.y * x_lam12 + p.z * x_lam13;
        const double new_y =
                p.freeze_out_time * x_lam02 + p.x * x_lam12 + p.y * x_lam22 + p.z * x_lam23;
        const double new_z =
                p.freeze_out_time * x_lam03 + p.x * x_lam13 + p.y * x_lam23 + p.z * x_lam33;
        const double new_p0     = p.p0 * x_lam00 + p.px * x_lam01 + p.py * x_lam02 + p.pz * x_lam03;
        const double new_px     = p.p0 * x_lam01 + p.px * x_lam11 + p.py * x_lam12 + p.pz * x_lam13;
        const double new_py     = p.p0 * x_lam02 + p.px * x_lam12 + p.py * x_lam22 + p.pz * x_lam23;
        const double new_pz     = p.p0 * x_lam03 + p.px * x_lam13 + p.py * x_lam23 + p.pz * x_lam33;
        boost_p.freeze_out_time = new_t;
        boost_p.x               = new_x;
        boost_p.y               = new_y;
        boost_p.z               = new_z;
        boost_p.p0              = new_p0;
        boost_p.px              = new_px;
        boost_p.py              = new_py;
        boost_p.pz              = new_pz;
    }
    return boost_p;
}

void update_d_mix_spv(const double pt, const double probability, const config_in &config_input,
                      std::vector<double> &d_mix_spv) {
    int npt = static_cast<int>(pt / config_input.d_mix_dpt);
    if (npt <= d_mix_spv.size()) {
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

DeutronData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                   const config_in &config_input, std::vector<double> &d_mix_spv) {
    DeutronData deutron_data{};
    // Rapidity checks
    if (std::abs(calculate_rap(p1.p0, p1.pz)) > config_input.rap_cut_nucl ||
        std::abs(calculate_rap(p2.p0, p2.pz)) > config_input.rap_cut_nucl) {
        return DeutronData{};
    }

    const double px_total = p1.px + p2.px;
    const double py_total = p1.py + p2.py;
    const double pz_total = p1.pz + p2.pz;
    const double p0_total = p1.p0 + p2.p0;
    //Deutron rapidity check
    if (std::abs(calculate_rap(p0_total, pz_total)) > config_input.rap_cut_coal) {
        return DeutronData{};
    }

    const double beta_x          = px_total / p0_total;
    const double beta_y          = py_total / p0_total;
    const double beta_z          = pz_total / p0_total;
    ParticleData boosted_proton  = lorentz_boost(beta_x, beta_y, beta_z, p1);
    ParticleData boosted_neutron = lorentz_boost(beta_x, beta_y, beta_z, p2);

    const double t_max = std::max(boosted_proton.freeze_out_time, boosted_neutron.freeze_out_time);

    boosted_proton.update_position(t_max);
    boosted_neutron.update_position(t_max);
    const auto [diff_x, diff_y, diff_z,
                diff_px, diff_py, diff_pz] =
            calculate_differences(boosted_proton, boosted_neutron);
    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);

    if (diff_p > config_input.cut_dp * 0.19733 / config_input.rms || diff_r > config_input.cut_dr * config_input.rms) {
        return DeutronData{};
    }
    deutron_data.probability = calculate_pro_deutron(diff_r, diff_p, config_input.rms);

    deutron_data.calculate_deutron_data(boosted_proton, boosted_neutron);
    //DeutronData boosted_deutron = lorentz_boost(-beta_x, -beta_y, -beta_z, deutron_data);
    const double pt = sqrt(px_total * px_total + py_total * py_total);

    update_d_mix_spv(pt, deutron_data.probability, config_input, d_mix_spv);

    return deutron_data;
}
void calculate_one_batch(const std::vector<ParticleData> &protons,
                         const std::vector<ParticleData> &neutrons,
                         const config_in &config_input, std::vector<double> &d_mix_spv,
                         const std::vector<double> &d_mix_ptv, double &batch_deutrons,
                         int eventsInBatch, std::vector<DeutronData> &deutrons) {
    batch_deutrons      = 0.0;
    const int mixEvents = eventsInBatch * eventsInBatch;
    auto compare        = [](const DeutronData &a, const DeutronData &b) {
        return a.probability < b.probability;
    };
    std::priority_queue<DeutronData, std::vector<DeutronData>, decltype(compare)> deuteron_queue(
            compare);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            DeutronData deutron = calculate_info_deutron(proton, neutron, config_input, d_mix_spv);
            batch_deutrons += deutron.probability;
            if (deutron.probability > 0) {
                deuteron_queue.push(deutron);
            }
        }
    }
    int expected_deutrons = static_cast<int>(batch_deutrons);
    double fraction       = batch_deutrons - expected_deutrons;
    if (dis(gen) < fraction) {
        expected_deutrons++;
    }
    while (!deuteron_queue.empty() && expected_deutrons > 0) {
        deutrons.push_back(deuteron_queue.top());
        deuteron_queue.pop();
        expected_deutrons--;
    }
    for (int k = 0; k < d_mix_spv.size(); ++k) {
        if (d_mix_ptv[k] > 0) {
            d_mix_spv[k] /= (2 * M_PI * d_mix_ptv[k] * config_input.d_mix_dpt * mixEvents);
        }
    }
}
void calculate_deuteron(const std::string &protonFile,
                        const std::string &neutronFile,
                        const config_in &configInput,
                        const int batchSize, std::vector<double> &dMixSpv,
                        const std::vector<double> &dMixPtv,
                        std::vector<DeutronData> &deutrons) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    read_batch_size(protonFile, batchSize, protonBatches);
    read_batch_size(neutronFile, batchSize, neutronBatches);
    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons     = batchData.particles;
        const auto &neutrons    = neutronBatches[batchNumber].particles;
        double batch_deutrons   = 0.0;
        const int eventsInBatch = batchData.eventCount;
        calculate_one_batch(protons,
                            neutrons,
                            configInput,
                            dMixSpv,
                            dMixPtv,
                            batch_deutrons,
                            eventsInBatch, deutrons);
        total_deutrons += batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
        std::cout << "Average number of deuteron per batch: " << batch_deutrons / eventsInBatch / eventsInBatch << std::endl;
    }
    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per event: " << average_deutrons << std::endl;

    output_deutrons(deutrons, "data/deuteron.dat");
    output_d_mix_spv(dMixSpv, dMixPtv, "tem/d_mix_spv.dat", total_batches);
}
