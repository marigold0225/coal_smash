//
// Created by mafu on 12/22/2023.
//
#include "../include/general_function.h"
#include <algorithm>
#include <random>
double samplingAndScalingFactor(const std::vector<ParticleData> &protons,
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

    auto kp     = static_cast<double>(protons.size());
    auto kn     = static_cast<double>(neutrons.size());
    auto kp_mix = static_cast<double>(protons_fraction.size());
    auto kn_mix = static_cast<double>(neutrons_fraction.size());

    double factor = (kp * (kp - 1) / (kp_mix * (kp_mix - 1)) * 2.0) *
                    (kn * (kn - 1) / (kn_mix * (kn_mix - 1)) * 2.0);
    //    double factor = (kp * kn) / (kp_mix * kn_mix);

    return factor;
}
double calculateParticleRapidity(const double p0, const double pz) {
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
void updateMomentumArray(const ParticleData &particle, const reactionConfig &config,
                         ptArray &pt_array, const RapidityMap &rapidityRange,
                         std::map<std::string, double> &clusterCountByRapidity) {
    double d_pt = config.dpt;
    int ptBins  = config.ptBins;
    double pt   = std::sqrt(particle.px * particle.px + particle.py * particle.py);
    double rap  = particle.getRapidity();
    for (auto &[label, range]: rapidityRange) {
        if (rap >= range.min && rap < range.max) {
            int npt = static_cast<int>(pt / d_pt);
            if (npt < ptBins) {
                pt_array[label][npt] += particle.probability;
            }
            clusterCountByRapidity[label] += particle.probability;
        }
    }
}
void weightedSampling(std::vector<ParticleData> &sample_particles,
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
std::tuple<double, double, double, double, double, double> TwoBodyJacobi(const ParticleData &p1,
                                                                         const ParticleData &p2) {
    double diff_x  = (p1.x - p2.x) / sqrt(2.0);
    double diff_y  = (p1.y - p2.y) / sqrt(2.0);
    double diff_z  = (p1.z - p2.z) / sqrt(2.0);
    double diff_px = (p1.px - p2.px) / sqrt(2.0);
    double diff_py = (p1.py - p2.py) / sqrt(2.0);
    double diff_pz = (p1.pz - p2.pz) / sqrt(2.0);
    return {diff_x, diff_y, diff_z, diff_px, diff_py, diff_pz};
}
std::tuple<double, double, double, double, double, double> fourBodyJacobi(const ParticleData &p1,
                                                                          const ParticleData &p2,
                                                                          const ParticleData &n1,
                                                                          const ParticleData &n2) {
    //p1+p2
    double dx_p1p2      = (p1.x - p2.x) / sqrt(2.0);
    double dy_p1p2      = (p1.y - p2.y) / sqrt(2.0);
    double dz_p1p2      = (p1.z - p2.z) / sqrt(2.0);
    double dpx_p1p2     = (p2.mass * p1.px - p1.mass * p2.px) / (p1.mass + p2.mass) * sqrt(2.0);
    double dpy_p1p2     = (p2.mass * p1.py - p1.mass * p2.py) / (p1.mass + p2.mass) * sqrt(2.0);
    double dpz_p1p2     = (p2.mass * p1.pz - p1.mass * p2.pz) / (p1.mass + p2.mass) * sqrt(2.0);
    double diff_dr_p1p2 = sqrt(dx_p1p2 * dx_p1p2 + dy_p1p2 * dy_p1p2 + dz_p1p2 * dz_p1p2);
    double diff_dp_p1p2 = sqrt(dpx_p1p2 * dpx_p1p2 + dpy_p1p2 * dpy_p1p2 + dpz_p1p2 * dpz_p1p2);
    //p1+p2+n1
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
    double diff_dr_p1p2n1 =
            sqrt(dx_p1p2n1 * dx_p1p2n1 + dy_p1p2n1 * dy_p1p2n1 + dz_p1p2n1 * dz_p1p2n1);
    double diff_dp_p1p2n1 =
            sqrt(dpx_p1p2n1 * dpx_p1p2n1 + dpy_p1p2n1 * dpy_p1p2n1 + dpz_p1p2n1 * dpz_p1p2n1);
    //p1+p2+n1+n2
    double dx_p1p2n1n2 = (p1.mass * p1.x + p2.mass * p2.x + n1.mass * n1.x -
                          (p1.mass + p2.mass + n1.mass) * n2.x) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dy_p1p2n1n2 = (p1.mass * p1.y + p2.mass * p2.y + n1.mass * n1.y -
                          (p1.mass + p2.mass + n1.mass) * n2.y) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dz_p1p2n1n2 = (p1.mass * p1.z + p2.mass * p2.z + n1.mass * n1.z -
                          (p1.mass + p2.mass + n1.mass) * n2.z) /
                         (p1.mass + p2.mass + n1.mass) * sqrt(3.0 / 4.0);
    double dpx_p1p2n1n2 =
            (n2.mass * (p1.px + p2.px + n1.px) - (p1.mass + p2.mass + n1.mass) * n2.px) /
            (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double dpy_p1p2n1n2 =
            (n2.mass * (p1.py + p2.py + n1.py) - (p1.mass + p2.mass + n1.mass) * n2.py) /
            (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double dpz_p1p2n1n2 =
            (n2.mass * (p1.pz + p2.pz + n1.pz) - (p1.mass + p2.mass + n1.mass) * n2.pz) /
            (p1.mass + p2.mass + n1.mass + n2.mass) * sqrt(4.0 / 3.0);
    double diff_dr_p1p2n1n2 =
            sqrt(dx_p1p2n1n2 * dx_p1p2n1n2 + dy_p1p2n1n2 * dy_p1p2n1n2 + dz_p1p2n1n2 * dz_p1p2n1n2);
    double diff_dp_p1p2n1n2 = sqrt(dpx_p1p2n1n2 * dpx_p1p2n1n2 + dpy_p1p2n1n2 * dpy_p1p2n1n2 +
                                   dpz_p1p2n1n2 * dpz_p1p2n1n2);
    return std::make_tuple(diff_dr_p1p2, diff_dp_p1p2, diff_dr_p1p2n1, diff_dp_p1p2n1,
                           diff_dr_p1p2n1n2, diff_dp_p1p2n1n2);
}
void calculateProtonPt(const std::string &protonFileName, const std::string &ptFileName,
                       ptArray &protonPt, const RapidityMap &rapidityRanges) {
    double d_pt = 0.2;
    int ptBins  = 10;

    std::map<std::string, int> protonCountsByRapidity;
    for (auto &[label, _]: rapidityRanges) {
        protonPt[label]               = std::vector<double>(ptBins, 0.0);
        protonCountsByRapidity[label] = 0;
    }
    BatchMap protonBatches;
    int total_events = 0;
    readBatchNuclei(protonFileName, 1, protonBatches);
    for (const auto &[batchNumber, batchData]: protonBatches) {
        const auto &protons = batchData.particles;
        int eventsInBatch   = batchData.eventCount;
        for (const auto &proton: protons) {
            double pt       = std::sqrt(proton.px * proton.px + proton.py * proton.py);
            double rapidity = proton.getRapidity();

            for (const auto &[label, range]: rapidityRanges) {
                if (rapidity >= range.min && rapidity < range.max) {
                    int npt = static_cast<int>(pt / d_pt);
                    if (npt < ptBins) {
                        protonPt[label][npt] += 1;
                    }
                    protonCountsByRapidity[label]++;
                    break;
                }
            }
        }
        total_events += eventsInBatch;
    }
    std::cout << "Total events: " << total_events << std::endl;
    std::ofstream output(ptFileName, std::ios::out);
    if (!output.is_open()) {
        throw std::runtime_error("Could not open file " + ptFileName);
    }
    for (auto &[label, pts]: protonPt) {
        output << "Rapidity range: " << label << ", Proton yield: "
               << static_cast<double>(protonCountsByRapidity[label]) / total_events << std::endl;
        for (int i = 0; i < ptBins; ++i) {
            double pt = d_pt / 2 + static_cast<double>(i) * d_pt;
            pts[i] /= (2 * M_PI * pt * d_pt * total_events);
            output << pt << "\t" << pts[i] << std::endl;
        }
        output << std::endl;
    }
    output.close();
}
RapidityMap defineRapidityRange() {
    RapidityMap rapidityRanges = {{"-0.1<y<0.0", {-0.1, 0.0}},   {"-0.2<y<-0.1", {-0.2, -0.1}},
                                  {"-0.3<y<-0.2", {-0.3, -0.2}}, {"-0.4<y<-0.3", {-0.4, -0.3}},
                                  {"-0.5<y<-0.4", {-0.5, -0.4}}, {"-0.6<y<-0.5", {-0.6, -0.5}},
                                  {"-0.7<y<-0.6", {-0.7, -0.6}}, {"-0.8<y<-0.7", {-0.8, -0.7}},
                                  {"-0.9<y<-0.8", {-0.9, -0.8}}, {"-1.0<y<-0.9", {-1.0, -0.9}}};
    return rapidityRanges;
}
void calculateClusterPt(const std::string &clusterFileName, const std::string &ptFileName,
                        const RapidityMap &rapidityRanges, const reactionConfig &clusterConfig) {
    double d_pt = clusterConfig.dpt;
    int ptBins  = clusterConfig.ptBins;
    ptArray clusterPt;
    std::map<std::string, double> clusterCountsByRapidity;
    for (auto &[label, _]: rapidityRanges) {
        clusterPt[label]               = std::vector<double>(ptBins, 0.0);
        clusterCountsByRapidity[label] = 0.0;
    }
    BatchMap clusterBatches;
    int total_Batch = 0;
    readBatchDeutrons(clusterFileName, clusterBatches);
    for (const auto &[batchNumber, batchData]: clusterBatches) {
        const auto &clusters = batchData.particles;

        int eventsInBatch = batchData.eventCount;
        int mixEvents     = eventsInBatch * eventsInBatch;
        //        int mixEvents     = eventsInBatch * eventsInBatch * eventsInBatch * eventsInBatch *
        //                        eventsInBatch * eventsInBatch * eventsInBatch * eventsInBatch;
        for (const auto &cluster: clusters) {
            double pt       = std::sqrt(cluster.px * cluster.px + cluster.py * cluster.py);
            double rapidity = cluster.getRapidity();

            for (const auto &[label, range]: rapidityRanges) {
                if (rapidity >= range.min && rapidity < range.max) {
                    int npt = static_cast<int>(pt / d_pt);
                    if (npt < ptBins) {
                        clusterPt[label][npt] += 1 / static_cast<double>(mixEvents);
                    }
                    clusterCountsByRapidity[label] += 1 / static_cast<double>(mixEvents);
                    break;
                }
            }
        }
        std::cout << "number of cluster per batch (Batch " << batchNumber
                  << "): " << clusters.size() << std::endl;
        total_Batch++;
    }
    std::cout << "Total Batch: " << total_Batch << std::endl;
    std::ofstream output(ptFileName, std::ios::out);
    if (!output.is_open()) {
        throw std::runtime_error("Could not open file " + ptFileName);
    }
    for (auto &[label, pts]: clusterPt) {
        output << "Rapidity range: " << label
               << ", Cluster yield: " << clusterCountsByRapidity[label] / total_Batch << std::endl;
        for (int i = 0; i < ptBins; ++i) {
            double pt = d_pt / 2 + i * d_pt;
            pts[i] /= (2 * M_PI * pt * d_pt * total_Batch);
            output << pt << "\t" << pts[i] << std::endl;
        }
        output << std::endl;
    }
    output.close();
}
void weightedSample(std::deque<ParticleData> &sample_particles,
                    std::deque<std::pair<ParticleData, double>> &potential_particles,
                    std::deque<double> &cumulative_probabilities, double number_of_particles,
                    std::mt19937 &gen, std::uniform_real_distribution<> &dis) {
    if (potential_particles.empty() || cumulative_probabilities.empty()) {
        return;
    }

    int expected_particles = static_cast<int>(number_of_particles);
    double fraction        = number_of_particles - expected_particles;
    if (dis(gen) < fraction) {
        expected_particles++;
    }

    //Generate a low-discrepancy sequence using van der Corput method
    //https://en.wikipedia.org/wiki/Low-discrepancy_sequence
    std::vector<double> random_numbers(expected_particles);
    for (int i = 1; i <= expected_particles; i++) {
        int reversed_bits = 0;
        int n             = i;
        while (n > 0) {
            reversed_bits = (reversed_bits << 1) | (n & 1);
            n >>= 1;
        }
        random_numbers[i - 1] =
                (i + reversed_bits / static_cast<double>(expected_particles)) / expected_particles;
    }

    //Resample the particles using the low-discrepancy sequence
    double upper_bound = cumulative_probabilities.back();
    auto it            = cumulative_probabilities.begin();
    for (double random_number: random_numbers) {
        random_number *= upper_bound;
        while (it != cumulative_probabilities.end() && *it < random_number) {
            it++;
        }
        if (it != cumulative_probabilities.end()) {
            size_t index = std::distance(cumulative_probabilities.begin(), it);
            sample_particles.push_back(potential_particles[index].first);
            cumulative_probabilities.erase(it);
            potential_particles.erase(potential_particles.begin() +
                                      static_cast<std::ptrdiff_t>(index));
        }
    }
}
