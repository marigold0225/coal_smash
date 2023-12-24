//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_deutron.h"
#include <algorithm>
#include <iomanip>
#include <queue>
#include <random>

ParticleData pNToDeutron(const ParticleData &p1, const ParticleData &p2,
                         const reactionConfig &deutronConfig, ptArray &pt_array,
                         const RapidityMap &rapidityRange,
                         std::map<std::string, double> &clusterCountByRapidity) {
    ParticleData deutron_data{};
    int cut_dr      = 5;
    int cut_dp      = 5;
    double rms      = deutronConfig.rms;
    double hbar2    = 0.19733 * 0.19733;
    double sig      = rms * sqrt(4. / 3.);
    double rap_nucl = deutronConfig.rap_cut_nucl;
    double rap_coal = deutronConfig.rap_cut_coal;
    // Rapidity checks
    if (std::abs(p1.getRapidity()) > rap_nucl || std::abs(p2.getRapidity()) > rap_nucl) {
        deutron_data.probability = 0.0;
        return deutron_data;
    }
    deutron_data.getTwobodyData(p1, p2);
    //Deutron rapidity check
    double deutron_rap = deutron_data.getRapidity();
    if (std::abs(deutron_rap) > rap_coal) {
        deutron_data.probability = 0.0;
        return deutron_data;
    }

    const double beta_x = deutron_data.px / deutron_data.p0;
    const double beta_y = deutron_data.py / deutron_data.p0;
    const double beta_z = deutron_data.pz / deutron_data.p0;

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

    if (diff_p > cut_dp * 0.19733 / sig || diff_r > cut_dr * sig) {
        deutron_data.probability = 0.0;
        return deutron_data;
    }
    deutron_data.probability =
            3. / 4. * 8 * exp(-diff_r * diff_r / sig / sig - diff_p * diff_p * sig * sig / hbar2);
    updateMomentumArray(deutron_data, deutronConfig, pt_array, rapidityRange,
                        clusterCountByRapidity);
    return deutron_data;
}

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const reactionConfig &deutronConfig,
                     ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_deutrons,
                     int eventsInBatch, std::vector<ParticleData> &deutrons,
                     std::map<std::string, double> &clusterCountByRapidity) {
    deutrons.clear();
    batch_deutrons      = 0.0;
    double d_pt         = deutronConfig.dpt;
    int ptBins          = deutronConfig.ptBins;
    const int mixEvents = eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_deutrons;
    std::vector<double> cumulative_probabilities;
    ptArray batch_pt_array;
    std::map<std::string, double> clusterCountOneBatch;
    for (auto &[label, _]: rapidityRange) {
        batch_pt_array[label]       = std::vector<double>(ptBins, 0.0);
        clusterCountOneBatch[label] = 0.0;
    }

    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            ParticleData deutron = pNToDeutron(proton, neutron, deutronConfig, batch_pt_array,
                                               rapidityRange, clusterCountOneBatch);
            if (deutron.probability > 0) {
                potential_deutrons.emplace_back(deutron, deutron.probability);
                batch_deutrons += deutron.probability;
                cumulative_probabilities.push_back(batch_deutrons);
            }
        }
    }

    weightedSampling(deutrons, potential_deutrons, cumulative_probabilities, batch_deutrons);

    for (auto &[label, pts]: batch_pt_array) {
        for (size_t k = 0; k < ptBins; ++k) {
            double pt = d_pt / 2 + static_cast<double>(k) * d_pt;
            pts[k] /= (2 * M_PI * pt * d_pt * mixEvents);
            pt_array[label][k] += pts[k];
        }
    }
    for (auto &[label, count]: clusterCountOneBatch) {
        clusterCountByRapidity[label] += count / mixEvents;
    }
}

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput, const RapidityMap &rapidityRanges) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    int batchSize         = configInput.mix_events;
    int ptBins            = configInput.deuteron.ptBins;
    ptArray pt_array;
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
        std::vector<ParticleData> batchDeutrons{};
        DeutronOneBatch(protons, neutrons, configInput.deuteron, pt_array, rapidityRanges,
                        num_batch_deutrons, eventsInBatch, batchDeutrons, deutronCountByRapidity);
        //        for (const auto &deuteron: batchDeutrons) {
        //            double rapidity = deuteron.getRapidity();
        //            for (const auto &[label, range]: rapidityRanges) {
        //                if (rapidity >= range.min && rapidity < range.max) {
        //                    deutronCountByRapidity[label] +=
        //                            1.0 / static_cast<double>(eventsInBatch) / eventsInBatch;
        //                    break;
        //                }
        //            }
        //        }
        if (output.is_open()) {
            outputCluster(batchDeutrons, eventsInBatch, output);
        }
        std::cout << "number of deutron per batch (Batch " << batchNumber
                  << "): " << batchDeutrons.size() << std::endl;
        total_deutrons += num_batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
    }
    output.close();
    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per event: " << average_deutrons << std::endl;
    outputPt(pt_array, deutronCountByRapidity, configInput.deuteron, ptFile, total_batches);
}
