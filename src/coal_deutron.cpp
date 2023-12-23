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
                         const RapidityMap &rapidityRange) {
    ParticleData deutron_data{};
    int cut_dr = 5;
    int cut_dp = 5;
    // gc = 2j+1/2^N
    // deutron : j = 1, N = 2
    // helium3 : j = 1/2, N = 3
    // alpha : j = 0, N = 4

    double rms      = deutronConfig.rms;
    double hbar2    = 0.19733 * 0.19733;
    double sig      = rms * sqrt(4. / 3.);
    double rap_nucl = deutronConfig.rap_cut_nucl;
    double rap_coal = deutronConfig.rap_cut_coal;
    double d_pt     = deutronConfig.dpt;
    int d_ptBins    = deutronConfig.ptBins;
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

    if (diff_p > cut_dp * 0.19733 / sig || diff_r > cut_dr * sig) {
        return ParticleData{};
    }
    deutron_data.probability =
            3. / 4. * 8 * exp(-diff_r * diff_r / sig / sig - diff_p * diff_p * sig * sig / hbar2);
    deutron_data.getTwobodyData(p1, p2);
    const double pt = sqrt(px_total * px_total + py_total * py_total);
    updateMomentumArray(pt, deutron_data.probability, d_pt, d_ptBins, deutron_rap, pt_array,
                        rapidityRange);
    return deutron_data;
}

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const reactionConfig &deutronConfig,
                     ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_deutrons,
                     int eventsInBatch, std::vector<ParticleData> &deutrons) {
    deutrons.clear();
    batch_deutrons      = 0.0;
    double d_pt         = deutronConfig.dpt;
    int ptBins          = deutronConfig.ptBins;
    const int mixEvents = eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_deutrons;
    std::vector<double> cumulative_probabilities;
    ptArray batch_pt_array;
    for (auto &[label, _]: rapidityRange) {
        batch_pt_array[label] = std::vector<double>(ptBins, 0.0);
    }

    for (const auto &proton: protons) {
        for (const auto &neutron: neutrons) {
            ParticleData deutron =
                    pNToDeutron(proton, neutron, deutronConfig, batch_pt_array, rapidityRange);
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
}

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput, ptArray &pt_array,
                      const RapidityMap &rapidityRanges) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches     = 0;
    int batchSize         = configInput.mix_events;
    int ptBins            = configInput.deuteron.ptBins;
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
        DeutronOneBatch(protons, neutrons, configInput.deuteron, pt_array, rapidityRanges,
                        num_batch_deutrons, eventsInBatch, batchDeutrons);
        total_deutrons += num_batch_deutrons / eventsInBatch / eventsInBatch;
        total_batches++;
        if (output.is_open()) {
            outputCluster(batchDeutrons, eventsInBatch, output);
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
    outputPt(pt_array, deutronCountByRapidity, configInput.deuteron, ptFile, total_batches);
}
