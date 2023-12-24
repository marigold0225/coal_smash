//
// Created by mafu on 12/23/2023.
//
#include "../include/coal_Be.h"


ParticleData he4He4ToBe(const ParticleData &he4_1, const ParticleData &he4_2,
                        const reactionConfig &BeConfig, ptArray &pt_array,
                        const RapidityMap &rapidityRange,
                        std::map<std::string, double> &clusterCountOneBatch) {
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
    // Rapidity checks
    ParticleData be_data{};
    if (std::abs(he4_1.getRapidity()) > rap_nucl || std::abs(he4_2.getRapidity()) > rap_nucl) {
        be_data.probability = 0.0;
        return be_data;
    }

    be_data.getTwobodyData(he4_1, he4_2);
    //Deutron rapidity check
    double be_rap = be_data.getRapidity();
    if (std::abs(be_rap) > rap_coal) {
        be_data.probability = 0.0;
        return be_data;
    }

    const double beta_x        = be_data.px / be_data.p0;
    const double beta_y        = be_data.py / be_data.p0;
    const double beta_z        = be_data.pz / be_data.p0;
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
        be_data.probability = 0.0;
        return be_data;
    }
    be_data.probability =
            1.0 / 4 * 8 * exp(-diff_r * diff_r / sig / sig - diff_p * diff_p * sig * sig / hbar2);
    updateMomentumArray(be_data, BeConfig, pt_array, rapidityRange, clusterCountOneBatch);
    return be_data;
}
void processBeOneBatch(const std::vector<ParticleData> &alpha, const reactionConfig &BeConfig,
                       ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_be,
                       int mixEvents, std::vector<ParticleData> &be,
                       std::map<std::string, double> &clusterCountByRapidity) {
    be.clear();
    batch_be    = 0.0;
    int ptBins  = BeConfig.ptBins;
    double d_pt = BeConfig.dpt;
    ptArray batch_pt_array;
    std::map<std::string, double> clusterCountOneBatch;
    for (auto &[label, _]: rapidityRange) {
        batch_pt_array[label]       = std::vector<double>(ptBins, 0.0);
        clusterCountOneBatch[label] = 0.0;
    }

    std::vector<std::pair<ParticleData, double>> potential_be;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < alpha.size(); i++) {
        for (size_t j = i + 1; j < alpha.size(); j++) {
            ParticleData be_particle = he4He4ToBe(alpha[i], alpha[j], BeConfig, batch_pt_array,
                                                  rapidityRange, clusterCountOneBatch);
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
    for (auto &[label, count]: clusterCountOneBatch) {
        clusterCountByRapidity[label] += count / mixEvents;
    }
}
void calculateBeAllBatch(const std::string &alphaFile, const std::string &beFile,
                         const std::string &ptFile, const reactionConfig &BeConfig,
                         const RapidityMap &rapidityRange) {
    double total_be   = 0.0;
    int total_batches = 0;
    int ptBins        = BeConfig.ptBins;
    ptArray pt_array;
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
        int mixEvents = eventInBatch * eventInBatch * eventInBatch * eventInBatch * eventInBatch *
                        eventInBatch * eventInBatch * eventInBatch;

        double batch_number_be = 0.0;
        std::vector<ParticleData> batch_be;
        processBeOneBatch(alphas, BeConfig, pt_array, rapidityRange, batch_number_be, mixEvents,
                          batch_be, clusterCountByRapidity);
        total_be += batch_number_be / mixEvents;
        total_batches++;
        if (beOutFile.is_open()) {
            outputCluster(batch_be, eventInBatch, beOutFile);
        }
        std::cout << "Average number of Be per batchNumber (Batch " << batchNumber
                  << "): " << batch_number_be / mixEvents << std::endl;
    }
    beOutFile.close();
    const double average_be = total_batches > 0 ? total_be / total_batches : 0.0;
    std::cout << "Average number of Be: " << average_be << std::endl;
    outputPt(pt_array, clusterCountByRapidity, BeConfig, ptFile, total_batches);
}
