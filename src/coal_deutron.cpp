//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_deutron.h"
#include <algorithm>
#include <iomanip>
#include <queue>
#include <random>


double calculateDeutronProbability(const double diff_r, const double diff_p, const double rms) {
    const double diff_r2 = diff_r * diff_r;
    const double diff_p2 = diff_p * diff_p;
    const double rms_2   = rms * rms / 0.75;
    const double probability =
            3. / 4. * 8 * std::exp(-diff_r2 / rms_2 - diff_p2 * rms_2 / 0.19733 / 0.19733);
    return probability;
}


ParticleData pNToDeutron(const ParticleData &p1, const ParticleData &p2,
                         const config_in &config_input,
                         std::map<std::string, std::vector<double>> &pt_array,
                         const std::map<std::string, RapidityRange> &rapidityRange) {
    ParticleData deutron_data{};
    double rap_nucl = config_input.d_rap_cut_nucl;
    double rap_coal = config_input.d_rap_cut_coal;
    double d_pt     = config_input.d_mix_dpt;
    // Rapidity checks
    if (std::abs(p1.getRapidity()) > rap_nucl || std::abs(p2.getRapidity()) > rap_nucl) {
        return ParticleData{};
    }

    const double px_total = p1.px + p2.px;
    const double py_total = p1.py + p2.py;
    const double pz_total = p1.pz + p2.pz;
    const double p0_total = p1.p0 + p2.p0;
    //Deutron rapidity check
    double deutron_rap = calculateParticleRapidity(p0_total, pz_total);
    if (std::abs(deutron_rap) > rap_coal) {
        return ParticleData{};
    }

    const double beta_x = px_total / p0_total;
    const double beta_y = py_total / p0_total;
    const double beta_z = pz_total / p0_total;

    ParticleData boosted_proton  = p1.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boosted_neutron = p2.lorentzBoost(beta_x, beta_y, beta_z);
    const double t_boost_max =
            std::max(boosted_proton.freeze_out_time, boosted_neutron.freeze_out_time);

    boosted_proton.updatePosition(t_boost_max);
    boosted_neutron.updatePosition(t_boost_max);
    const auto [diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz] =
            TwoBodyJacobi(boosted_proton, boosted_neutron);
    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);

    if (diff_p > config_input.cut_dp * 0.19733 / config_input.rms / sqrt(4. / 3.) ||
        diff_r > config_input.cut_dr * config_input.rms * sqrt(4. / 3.)) {
        return ParticleData{};
    }
    deutron_data.probability = calculateDeutronProbability(diff_r, diff_p, config_input.rms);
    deutron_data.getTwobodyData(p1, p2);
    const double pt = sqrt(px_total * px_total + py_total * py_total);
    updateMomentumArray(pt, deutron_data.probability, d_pt, deutron_rap, pt_array, rapidityRange);
    return deutron_data;
}

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const config_in &config_input,
                     std::map<std::string, std::vector<double>> &pt_array,
                     const std::map<std::string, RapidityRange> &rapidityRange,
                     double &batch_deutrons, int eventsInBatch,
                     std::vector<ParticleData> &deutrons) {
    deutrons.clear();
    batch_deutrons      = 0.0;
    double d_pt         = config_input.d_mix_dpt;
    int ptBins          = 10;
    const int mixEvents = eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_deutrons;
    std::vector<double> cumulative_probabilities;

    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            ParticleData deutron =
                    pNToDeutron(proton, neutron, config_input, pt_array, rapidityRange);
            if (deutron.probability > 0) {
                potential_deutrons.emplace_back(deutron, deutron.probability);
                batch_deutrons += deutron.probability;
                cumulative_probabilities.push_back(batch_deutrons);
            }
        }
    }

    weightedSampling(deutrons, potential_deutrons, cumulative_probabilities, batch_deutrons);

    for (auto &[label, pts]: pt_array) {
        for (size_t k = 0; k < ptBins; ++k) {
            double pt = d_pt / 2 + static_cast<double>(k) * d_pt;
            pts[k] /= (2 * M_PI * pt * d_pt * mixEvents);
        }
    }
}

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput,
                      std::map<std::string, std::vector<double>> &pt_array,
                      std::map<std::string, RapidityRange> &rapidityRanges) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    int batchSize         = configInput.mix_events;
    double d_pt           = configInput.d_mix_dpt;
    int ptBins            = 10;
    std::map<std::string, double> deutronCountByRapidity;
    for (auto &[label, _]: rapidityRanges) {
        pt_array[label]               = std::vector<double>(ptBins, 0.0);
        deutronCountByRapidity[label] = 0.0;
    }

    readBatchNuclei(protonFile, batchSize, protonBatches);
    readBatchNuclei(neutronFile, batchSize, neutronBatches);
    std::ofstream output(deutronFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons       = batchData.particles;
        const auto &neutrons      = neutronBatches[batchNumber].particles;
        double num_batch_deutrons = 0.0;
        const int eventsInBatch   = batchData.eventCount;
        std::vector<ParticleData> batchDeutrons;
        DeutronOneBatch(protons, neutrons, configInput, pt_array, rapidityRanges,
                        num_batch_deutrons, eventsInBatch, batchDeutrons);
        total_deutrons += num_batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
        if (output.is_open()) {
            outputCluster(batchDeutrons, output);
        }
        std::cout << "Average number of deuteron per batch (Batch " << batchNumber
                  << "): " << num_batch_deutrons / eventsInBatch / eventsInBatch << std::endl;

        for (const auto &deuteron: batchDeutrons) {
            double rapidity = deuteron.getRapidity();
            for (const auto &[label, range]: rapidityRanges) {
                if (rapidity >= range.min && rapidity < range.max) {
                    deutronCountByRapidity[label] +=
                            1.0 / static_cast<double>(eventsInBatch) / eventsInBatch;
                    break;
                }
            }
        }
    }
    output.close();
    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per event: " << average_deutrons << std::endl;
    outputPt(pt_array, deutronCountByRapidity, d_pt, ptBins, ptFile, total_batches);
}
