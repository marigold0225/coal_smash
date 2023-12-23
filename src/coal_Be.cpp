//
// Created by mafu on 12/23/2023.
//
#include "../include/coal_Be.h"


ParticleData he4He4ToBe(const ParticleData &he4_1, const ParticleData &he4_2,
                        const reactionConfig &BeConfig, ptArray &pt_array,
                        const RapidityMap &rapidityRange) {
    ParticleData be_data{};
    int cut_dr = 5;
    int cut_dp = 5;
    // gc = 2j+1/2^N
    // deutron : j = 1, N = 2
    // helium3 : j = 1/2, N = 3
    // alpha : j = 0, N = 4
    double rms      = BeConfig.rms;
    double hbar2    = 0.19733 * 0.19733;
    double sig      = rms * sqrt(4. / 3.);
    double rap_nucl = BeConfig.rap_cut_nucl;
    double rap_coal = BeConfig.rap_cut_coal;
    double d_pt     = BeConfig.dpt;
    // Rapidity checks
    if (std::abs(he4_1.getRapidity()) > rap_nucl || std::abs(he4_2.getRapidity()) > rap_nucl) {
        return ParticleData{};
    }

    const double px_total = he4_1.px + he4_2.px;
    const double py_total = he4_1.py + he4_2.py;
    const double pz_total = he4_1.pz + he4_2.pz;
    const double p0_total = he4_1.p0 + he4_2.p0;
    //Deutron rapidity check
    double be_rap = calculateParticleRapidity(p0_total, pz_total);
    if (std::abs(be_rap) > rap_coal) {
        return ParticleData{};
    }

    const double beta_x = px_total / p0_total;
    const double beta_y = py_total / p0_total;
    const double beta_z = pz_total / p0_total;

    ParticleData boosted_he4_1 = he4_1.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boosted_he4_2 = he4_2.lorentzBoost(beta_x, beta_y, beta_z);
    const double t_boost_max =
            std::max(boosted_he4_1.freeze_out_time, boosted_he4_2.freeze_out_time);

    boosted_he4_1.updatePosition(t_boost_max);
    boosted_he4_2.updatePosition(t_boost_max);
    const auto [diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz] =
            TwoBodyJacobi(boosted_he4_1, boosted_he4_2);
    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);

    if (diff_p > cut_dp * 0.19733 / sig || diff_r > cut_dr * sig) {
        return ParticleData{};
    }
    be_data.probability =
            1. / 4. * 8 * exp(-diff_r * diff_r / sig / sig - diff_p * diff_p * sig * sig / hbar2);
    be_data.getTwobodyData(he4_1, he4_2);
    const double pt = sqrt(diff_px * diff_px + diff_py * diff_py);
    updateMomentumArray(pt, be_data.probability, d_pt, BeConfig.ptBins, be_rap, pt_array,
                        rapidityRange);
    return be_data;
}
void processBeOneBatch(const std::vector<ParticleData> &alpha, const reactionConfig &BeConfig,
                       ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_be,
                       double mixEvents, std::vector<ParticleData> &be) {
    be.clear();
    batch_be    = 0.0;
    int ptBins  = BeConfig.ptBins;
    double d_pt = BeConfig.dpt;
    std::map<std::string, std::vector<double>> batch_pt_array;
    for (auto &[label, _]: rapidityRange) {
        batch_pt_array[label] = std::vector<double>(ptBins, 0.0);
    }

    std::vector<std::pair<ParticleData, double>> potential_be;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < alpha.size(); i++) {
        for (size_t j = i + 1; j < alpha.size(); j++) {
            ParticleData be_particle =
                    he4He4ToBe(alpha[i], alpha[j], BeConfig, batch_pt_array, rapidityRange);
            if (be_particle.probability > 0) {
                potential_be.emplace_back(be_particle, be_particle.probability);
                batch_be += be_particle.probability;
                cumulated_probabilities.push_back(batch_be);
            }
        }
    }

    weightedSampling(be, potential_be, cumulated_probabilities, batch_be);

    for (auto &[label, pts]: batch_pt_array) {
        for (size_t k = 0; k < ptBins; ++k) {
            double pt = d_pt / 2 + static_cast<double>(k) * d_pt;
            if (pt > 0) {
                pts[k] = pts[k] / (2 * M_PI * pt * d_pt * mixEvents);
                pt_array[label][k] += pts[k];
            }
        }
    }
}
void calculateBeAllBatch(const std::string &alphaFile, const std::string &beFile,
                         const std::string &ptFile, const reactionConfig &BeConfig,
                         ptArray &pt_array, const RapidityMap &rapidityRange) {
    double total_be   = 0.0;
    int total_batches = 0;
    int ptBins        = BeConfig.ptBins;
    std::map<std::string, double> clusterCountByRapidity;
    for (auto &[label, _]: rapidityRange) {
        pt_array[label]               = std::vector<double>(ptBins, 0.0);
        clusterCountByRapidity[label] = 0.0;
    }

    BatchMap alphaBatches;
    readBatchDeutrons(alphaFile, alphaBatches);

    std::ofstream beOutFile(beFile, std::ios::out);

    for (const auto &[batchNumber, alpha]: alphaBatches) {
        const auto &alphas = alpha.particles;
        int eventInBatch   = alpha.eventCount;
        double mixEvents   = std::pow(eventInBatch, 8);

        double batch_number_be = 0.0;
        std::vector<ParticleData> batch_be;
        processBeOneBatch(alphas, BeConfig, pt_array, rapidityRange, batch_number_be, mixEvents,
                          batch_be);
        total_be += batch_number_be / mixEvents;
        total_batches++;
        if (beOutFile.is_open()) {
            outputCluster(batch_be, eventInBatch, beOutFile);
        }
        std::cout << "Average number of Be per batchNumber (Batch " << batchNumber
                  << "): " << batch_number_be / mixEvents << std::endl;
        for (const auto &be: batch_be) {
            double rapidity = be.getRapidity();
            for (auto &[label, range]: rapidityRange) {
                if (rapidity >= range.min && rapidity < range.max) {
                    clusterCountByRapidity[label] += 1.0 / static_cast<double>(mixEvents);
                    break;
                }
            }
        }
    }
    beOutFile.close();
    const double average_be = total_batches > 0 ? total_be / total_batches : 0.0;
    std::cout << "Average number of Be: " << average_be << std::endl;
    outputPt(pt_array, clusterCountByRapidity, BeConfig, ptFile, total_batches);
}
