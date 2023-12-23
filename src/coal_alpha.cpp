//
// Created by mafu on 12/15/2023.
//
#include "../include/coal_alpha.h"
#include <algorithm>
#include <filesystem>
#include <queue>
#include <random>

ParticleData dDToAlpha(const ParticleData &d1, const ParticleData &d2,
                       const reactionConfig &alphaConfig, ptArray &pt_array,
                       const RapidityMap &rapidityRange) {
    int cut_dr      = 5;
    int cut_dp      = 5;
    double rms      = alphaConfig.rms;
    double hbar2    = 0.19733 * 0.19733;
    double sig      = rms * sqrt(4. / 3.);
    double rap_nucl = alphaConfig.rap_cut_nucl;
    double rap_coal = alphaConfig.rap_cut_coal;
    double d_pt     = alphaConfig.dpt;
    ParticleData alpha{};
    //Rapidity check

    if (std::abs(d1.getRapidity()) > rap_nucl || std::abs(d2.getRapidity()) > rap_nucl) {
        return ParticleData{};
    }
    //alpha rapidity check
    const double px_total = d1.px + d2.px;
    const double py_total = d1.py + d2.py;
    const double pz_total = d1.pz + d2.pz;
    const double p0_total = d1.p0 + d2.p0;

    const double rapidity_alpha = calculateParticleRapidity(p0_total, pz_total);
    if (std::abs(rapidity_alpha) > rap_coal) {
        return ParticleData{};
    }
    //lorentz boost calculate
    const double beta_x   = px_total / p0_total;
    const double beta_y   = py_total / p0_total;
    const double beta_z   = pz_total / p0_total;
    ParticleData boost_d1 = d1.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boost_d2 = d2.lorentzBoost(beta_x, beta_y, beta_z);
    //find t_boost_max in d1 d2
    double t_boost_max = std::max(boost_d1.freeze_out_time, boost_d2.freeze_out_time);
    //update position
    boost_d1.updatePosition(t_boost_max);
    boost_d2.updatePosition(t_boost_max);
    //Calculate the position difference and momentum difference between two particles
    auto [diff_dx, diff_dy, diff_dz, diff_dpx, diff_dpy, diff_dpz] =
            TwoBodyJacobi(boost_d1, boost_d2);

    double diff_dr = sqrt(diff_dx * diff_dx + diff_dy * diff_dy + diff_dz * diff_dz);
    double diff_dp = sqrt(diff_dpx * diff_dpx + diff_dpy * diff_dpy + diff_dpz * diff_dpz);

    if (diff_dr > cut_dr * sig || diff_dp > cut_dp * 0.19733 / sig) {
        return ParticleData{};
    }
    alpha.probability = 1.0 / 4 * 8 *
                        exp(-diff_dr * diff_dr / sig / sig - diff_dp * diff_dp * sig * sig / hbar2);
    alpha.getTwobodyData(d1, d2);
    const double pt = sqrt(px_total * px_total + py_total * py_total);
    updateMomentumArray(pt, alpha.probability, d_pt, alphaConfig.ptBins, rapidity_alpha, pt_array,
                        rapidityRange);
    return alpha;
}


void processAlphaOneBatch2(const std::vector<ParticleData> &deutrons,
                           const reactionConfig &alphaConfig, ptArray &pt_array,
                           const RapidityMap &rapidityRange, double &batch_alpha, int mixEvents,
                           std::vector<ParticleData> &alpha) {
    alpha.clear();
    batch_alpha = 0.0;
    int ptBins  = alphaConfig.ptBins;
    double d_pt = alphaConfig.dpt;
    std::map<std::string, std::vector<double>> batch_pt_array;
    for (auto &[label, _]: rapidityRange) {
        batch_pt_array[label] = std::vector<double>(ptBins, 0.0);
    }

    std::vector<std::pair<ParticleData, double>> potential_alpha;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < deutrons.size(); i++) {
        for (size_t j = i + 1; j < deutrons.size(); j++) {
            ParticleData alpha_particle =
                    dDToAlpha(deutrons[i], deutrons[j], alphaConfig, batch_pt_array, rapidityRange);
            if (alpha_particle.probability > 0) {
                potential_alpha.emplace_back(alpha_particle, alpha_particle.probability);
                batch_alpha += alpha_particle.probability;
                cumulated_probabilities.push_back(batch_alpha);
            }
        }
    }

    weightedSampling(alpha, potential_alpha, cumulated_probabilities, batch_alpha);

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
void calculateAlphaAllBatch2(const std::string &deuteronFile, const std::string &alphaFile,
                             std::string &ptFile, const reactionConfig &alphaConfig,
                             ptArray &pt_array, const RapidityMap &rapidityRange) {
    double total_alpha = 0.0;
    int total_batches  = 0;
    int ptBins         = alphaConfig.ptBins;
    std::map<std::string, double> clusterCountByRapidity;
    for (auto &[label, _]: rapidityRange) {
        pt_array[label]               = std::vector<double>(ptBins, 0.0);
        clusterCountByRapidity[label] = 0.0;
    }

    BatchMap deuteronBatches;
    readBatchDeutrons(deuteronFile, deuteronBatches);

    std::ofstream alphaFileOut(alphaFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: deuteronBatches) {
        const auto &deutrons = batchData.particles;
        int eventInBatch     = batchData.eventCount;
        int mixEvents        = eventInBatch * eventInBatch * eventInBatch * eventInBatch;

        double batch_number_alpha = 0.0;
        std::vector<ParticleData> batch_alpha;
        processAlphaOneBatch2(deutrons, alphaConfig, pt_array, rapidityRange, batch_number_alpha,
                              mixEvents, batch_alpha);
        total_alpha += batch_number_alpha / mixEvents;
        total_batches++;
        if (alphaFileOut.is_open()) {
            outputCluster(batch_alpha, eventInBatch, alphaFileOut);
        }
        std::cout << "Average number of alpha per batch (Batch " << batchNumber
                  << "): " << batch_number_alpha / mixEvents << std::endl;
        for (const auto &alpha: batch_alpha) {
            double rapidity = alpha.getRapidity();
            for (const auto &[label, range]: rapidityRange) {
                if (rapidity >= range.min && rapidity < range.max) {
                    clusterCountByRapidity[label] += 1.0 / static_cast<double>(mixEvents);
                    break;
                }
            }
        }
    }

    alphaFileOut.close();
    const double average_alpha = total_batches > 0 ? total_alpha / total_batches : 0.0;
    std::cout << "Average number of alpha: " << average_alpha << std::endl;

    outputPt(pt_array, clusterCountByRapidity, alphaConfig, ptFile, total_batches);
}

ParticleData pPNNToAlpha(const ParticleData &p1, const ParticleData &p2, const ParticleData &n1,
                         const ParticleData &n2, const config_in &config_input, ptArray &pt_array,
                         const RapidityMap &rapidityRange) {
    double sig1                  = config_input.alpha.rms;
    double sig2                  = config_input.alpha.rms;
    double sig3                  = config_input.alpha.rms;
    constexpr const double hbar2 = 0.19733 * 0.19733;

    double rap_nucl = config_input.alpha.rap_cut_nucl;
    double rap_coal = config_input.alpha.rap_cut_coal;

    ParticleData alpha_particle{};
    //Rapidity check
    double rapidity_proton1  = calculateParticleRapidity(p1.p0, p1.pz);
    double rapidity_proton2  = calculateParticleRapidity(p2.p0, p2.pz);
    double rapidity_neutron1 = calculateParticleRapidity(n1.p0, n1.pz);
    double rapidity_neutron2 = calculateParticleRapidity(n2.p0, n2.pz);
    if (std::abs(rapidity_proton1) > rap_nucl || std::abs(rapidity_proton2) > rap_nucl ||
        std::abs(rapidity_neutron1) > rap_nucl || std::abs(rapidity_neutron2) > rap_nucl) {
        return ParticleData{};
    }
    //alpha rapidity check
    double px_total       = p1.px + p2.px + n1.px + n2.px;
    double py_total       = p1.py + p2.py + n1.py + n2.py;
    double pz_total       = p1.pz + p2.pz + n1.pz + n2.pz;
    double p0_total       = p1.p0 + p2.p0 + n1.p0 + n2.p0;
    double rapidity_alpha = calculateParticleRapidity(p0_total, pz_total);
    if (std::abs(rapidity_alpha) > rap_coal) {
        return ParticleData{};
    }
    //lorentz boost calculate
    const double beta_x   = px_total / p0_total;
    const double beta_y   = py_total / p0_total;
    const double beta_z   = pz_total / p0_total;
    ParticleData boost_p1 = p1.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boost_p2 = p2.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boost_n1 = n1.lorentzBoost(beta_x, beta_y, beta_z);
    ParticleData boost_n2 = n2.lorentzBoost(beta_x, beta_y, beta_z);
    //find t_max in p1 p2 n1 n2
    double t_max = std::max(std::max(boost_p1.freeze_out_time, boost_p2.freeze_out_time),
                            std::max(boost_n1.freeze_out_time, boost_n2.freeze_out_time));
    //update position
    boost_p1.updatePosition(t_max);
    boost_p2.updatePosition(t_max);
    boost_n1.updatePosition(t_max);
    boost_n2.updatePosition(t_max);
    //Calculate the position difference and momentum difference between two particles two by two
    auto [diff_dr_p1p2, diff_dp_p1p2, diff_dr_p1p2n1, diff_dp_p1p2n1, diff_dr_p1p2n1n2,
          diff_dp_p1p2n1n2] = fourBodyJacobi(boost_p1, boost_p2, boost_n1, boost_n2);

    if (diff_dr_p1p2 > config_input.cut_dr * sig1 ||
        diff_dp_p1p2 > config_input.cut_dp * 0.19733 / sig1 ||
        diff_dr_p1p2n1 > config_input.cut_dr * sig2 ||
        diff_dp_p1p2n1 > config_input.cut_dp * 0.19733 / sig2 ||
        diff_dr_p1p2n1n2 > config_input.cut_dr * sig3 ||
        diff_dp_p1p2n1n2 > config_input.cut_dp * 0.19733 / sig3) {
        return ParticleData{};
    }

    double diff_dr = diff_dr_p1p2 * diff_dr_p1p2 / sig1 / sig1 +
                     diff_dp_p1p2 * diff_dp_p1p2 * sig1 * sig1 / hbar2 +
                     diff_dr_p1p2n1 * diff_dr_p1p2n1 / sig2 / sig2 +
                     diff_dp_p1p2n1 * diff_dp_p1p2n1 * sig2 * sig2 / hbar2 +
                     diff_dr_p1p2n1n2 * diff_dr_p1p2n1n2 / sig3 / sig3 +
                     diff_dp_p1p2n1n2 * diff_dp_p1p2n1n2 * sig3 * sig3 / hbar2;
    if (diff_dr < 40) {
        alpha_particle.probability = 1.0 / 16 * 8 * 8 * 8 * exp(-diff_dr);
    } else {
        alpha_particle.probability = 0;
    }
    alpha_particle.getFourbodyData(boost_p1, boost_p2, boost_n1, boost_n2);
    ParticleData alpha_particle_boost = alpha_particle.lorentzBoost(-beta_x, -beta_y, -beta_z);
    const double pt                   = sqrt(px_total * px_total + py_total * py_total);
    updateMomentumArray(pt, alpha_particle_boost.probability, config_input.alpha.dpt,
                        config_input.alpha.ptBins, rapidity_alpha, pt_array, rapidityRange);
    return alpha_particle_boost;
}


void processAlphaOneBatch4(const std::vector<ParticleData> &protons,
                           const std::vector<ParticleData> &neutrons, const config_in &config_input,
                           ptArray &pt_array, const RapidityMap &rapidityRang, double &batch_alpha,
                           int eventsInBatch, std::vector<ParticleData> &alpha) {
    alpha.clear();
    batch_alpha = 0.0;
    int ptBins  = 10;
    double d_pt = config_input.alpha.dpt;

    std::vector<ParticleData> protons_fraction, neutrons_fraction;
    double factor = samplingAndScalingFactor(protons, neutrons, protons_fraction, neutrons_fraction,
                                             config_input);

    const int mixEvents = eventsInBatch * eventsInBatch * eventsInBatch * eventsInBatch;
    std::vector<std::pair<ParticleData, double>> potential_alpha;
    std::vector<double> cumulated_probabilities;

    for (size_t i = 0; i < protons_fraction.size(); i++) {
        for (size_t j = i + 1; j < protons_fraction.size(); j++) {
            for (size_t k = 0; k < neutrons_fraction.size(); k++) {
                for (size_t l = k + 1; l < neutrons_fraction.size(); l++) {
                    ParticleData alpha_particle = pPNNToAlpha(
                            protons_fraction[i], protons_fraction[j], neutrons_fraction[k],
                            neutrons_fraction[l], config_input, pt_array, rapidityRang);
                    if (alpha_particle.probability > 0) {
                        //                     potential_alpha.emplace_back(alpha_particle, alpha_particle.probability);
                        batch_alpha += alpha_particle.probability * factor;
                        //                        cumulated_probabilities.push_back(batch_alpha);
                    }
                }
            }
        }
    }
    //    weightedSampling(alpha, potential_alpha, cumulated_probabilities, batch_alpha);

    for (auto &[label, pts]: pt_array) {
        for (size_t k = 0; k < ptBins; ++k) {
            double pt = d_pt / 2 + static_cast<double>(k) * d_pt;
            if (pt > 0) {
                pts[k] = pts[k] / (2 * M_PI * pt * d_pt * mixEvents) * factor;
            }
        }
    }
}

void calculateAlphaAllBatch4(const std::string &protonFile, const std::string &neutronFile,
                             const std::string &alphaFile, std::string &ptFile,
                             const config_in &configInput, ptArray &pt_array,
                             const RapidityMap &rapidityRang) {
    BatchMap protonBatches, neutronBatches;
    double total_alpha = 0.0;
    int total_batches  = 0;
    int batchSize      = configInput.mix_events;
    int ptBins         = 10;
    std::map<std::string, double> clusterCountByRapidity;
    for (auto &[label, _]: rapidityRang) {
        pt_array[label]               = std::vector<double>(ptBins, 0.0);
        clusterCountByRapidity[label] = 0.0;
    }
    readBatchNuclei(protonFile, batchSize, protonBatches);
    readBatchNuclei(neutronFile, batchSize, neutronBatches);

    std::ofstream alphaFileOut(alphaFile, std::ios::out);

    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons      = batchData.particles;
        const auto &neutrons     = neutronBatches[batchNumber].particles;
        int eventsInBatch        = batchData.eventCount;
        double total_batch_alpha = 0.0;
        //        std::vector<ParticleData> batch_average_alpha;

        for (int repetition = 0; repetition < 10; ++repetition) {

            double batch_number_alpha = 0.0;
            std::vector<ParticleData> batch_alpha;
            processAlphaOneBatch4(protons, neutrons, configInput, pt_array, rapidityRang,
                                  batch_number_alpha, eventsInBatch, batch_alpha);
            total_batch_alpha += batch_number_alpha / eventsInBatch / eventsInBatch /
                                 eventsInBatch / eventsInBatch;
            std::cout << "Batch " << batchNumber << " repetition " << repetition << " alpha: "
                      << batch_number_alpha / eventsInBatch / eventsInBatch / eventsInBatch /
                                 eventsInBatch
                      << std::endl;
        }
        total_alpha += total_batch_alpha / 10;
        total_batches++;
        //        if (alphaFileOut.is_open()) {
        //            outputCluster(batch_average_alpha, alphaFileOut);
        //        }
        //        std::cout << "average number of alpha per batch:" << total_batch_alpha / 10 << std::endl;
        for (auto &[label, pts]: pt_array) {
            for (size_t i = 0; i < ptBins; ++i) {
                pts[i] /= 10;
            }
        }
    }
    alphaFileOut.close();
    const double average_alpha = total_alpha > 0 ? total_alpha / total_batches : 0.0;
    std::cout << "average number of alpha:" << average_alpha << std::endl;
    outputPt(pt_array, clusterCountByRapidity, configInput.alpha, ptFile, total_batches);
}
