//
// Created by mafu on 12/13/2023.
//

#include "../include/coal.h"
#include <iomanip>
#include <ranges>

void calculate_freeze_position(ParticleData &p) {
    const double vx = p.px / p.p0;
    const double vy = p.py / p.p0;
    const double vz = p.pz / p.p0;

    const double dt = p.t - p.freeze_out_time;
    const double dx = vx * dt;
    const double dy = vy * dt;
    const double dz = vz * dt;

    p.x -= dx;
    p.y -= dy;
    p.z -= dz;
}

ParticleData lorentz_boost(const double beta_x, const double beta_y, const double beta_z,
                           const ParticleData &p) {
    ParticleData boost_p = p;
    if (const double beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z; beta2 > 1.0e-5) {
        const double gamma = 1.0 / sqrt(1.0 - beta2);

        const double xlam00 = gamma;
        const double xlam01 = -gamma * beta_x;
        const double xlam02 = -gamma * beta_y;
        const double xlam03 = -gamma * beta_z;
        const double xlam11 = 1.0 + (gamma - 1.0) * beta_x * beta_x / beta2;
        const double xlam22 = 1.0 + (gamma - 1.0) * beta_y * beta_y / beta2;
        const double xlam33 = 1.0 + (gamma - 1.0) * beta_z * beta_z / beta2;
        const double xlam12 = (gamma - 1.0) * beta_x * beta_y / beta2;
        const double xlam13 = (gamma - 1.0) * beta_x * beta_z / beta2;
        const double xlam23 = (gamma - 1.0) * beta_y * beta_z / beta2;

        const double new_t =
            p.freeze_out_time * xlam00 + p.x * xlam01 + p.y * xlam02 + p.z * xlam03;
        const double new_x =
            p.freeze_out_time * xlam01 + p.x * xlam11 + p.y * xlam12 + p.z * xlam13;
        const double new_y =
            p.freeze_out_time * xlam02 + p.x * xlam12 + p.y * xlam22 + p.z * xlam23;
        const double new_z =
            p.freeze_out_time * xlam03 + p.x * xlam13 + p.y * xlam23 + p.z * xlam33;

        const double new_p0 = p.p0 * xlam00 + p.px * xlam01 + p.py * xlam02 + p.pz * xlam03;
        const double new_px = p.p0 * xlam01 + p.px * xlam11 + p.py * xlam12 + p.pz * xlam13;
        const double new_py = p.p0 * xlam02 + p.px * xlam12 + p.py * xlam22 + p.pz * xlam23;
        const double new_pz = p.p0 * xlam03 + p.px * xlam13 + p.py * xlam23 + p.pz * xlam33;

        boost_p.freeze_out_time = new_t;
        boost_p.x = new_x;
        boost_p.y = new_y;
        boost_p.z = new_z;
        boost_p.p0 = new_p0;
        boost_p.px = new_px;
        boost_p.py = new_py;
        boost_p.pz = new_pz;
    }
    return boost_p;
}

double calculate_pro(const double diff_r, const double diff_p, const double rms) {
    const double diff_r2 = diff_r * diff_r;
    const double diff_p2 = diff_p * diff_p;
    const double rms_2 = rms * rms;
    const double probility =
        3. / 4. * 8 * std::exp(-diff_r2 / rms_2 - diff_p2 * rms_2 / 0.19733 / 0.19733);
    return probility;
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

double performCalculations(const ParticleData &p1, const ParticleData &p2,
                           const config_in &config_input, std::vector<double> &d_mix_spv) {
    if (const double proton_rap = calculate_rap(p1.p0, p1.pz);
        std::abs(proton_rap) > config_input.rap_cut_nucl) {
        return 0;
    }
    if (const double neutron_rap = calculate_rap(p2.p0, p2.pz);
        std::abs(neutron_rap) > config_input.rap_cut_nucl) {
        return 0;
    }

    const double px_total = p1.px + p2.px;
    const double py_total = p1.py + p2.py;
    const double pz_total = p1.pz + p2.pz;
    const double p0_total = p1.p0 + p2.p0;

    if (const double deutron_rap = calculate_rap(p0_total, pz_total);
        std::abs(deutron_rap) > config_input.rap_cut_coal) {
        return 0;
    }

    const double beta_x = px_total / p0_total;
    const double beta_y = py_total / p0_total;
    const double beta_z = pz_total / p0_total;

    ParticleData boosted_proton = lorentz_boost(beta_x, beta_y, beta_z, p1);
    ParticleData boosted_neutron = lorentz_boost(beta_x, beta_y, beta_z, p2);

    const double tmax = std::max(boosted_proton.freeze_out_time, boosted_neutron.freeze_out_time);

    boosted_neutron.x +=
        (tmax - boosted_neutron.freeze_out_time) * boosted_neutron.px / boosted_neutron.p0;
    boosted_neutron.y +=
        (tmax - boosted_neutron.freeze_out_time) * boosted_neutron.py / boosted_neutron.p0;
    boosted_neutron.z +=
        (tmax - boosted_neutron.freeze_out_time) * boosted_neutron.pz / boosted_neutron.p0;
    boosted_proton.x +=
        (tmax - boosted_proton.freeze_out_time) * boosted_proton.px / boosted_proton.p0;
    boosted_proton.y +=
        (tmax - boosted_proton.freeze_out_time) * boosted_proton.py / boosted_proton.p0;
    boosted_proton.z +=
        (tmax - boosted_proton.freeze_out_time) * boosted_proton.pz / boosted_proton.p0;

    const double diff_x = (boosted_proton.x - boosted_neutron.x) / sqrt(2.0);
    const double diff_y = (boosted_proton.y - boosted_neutron.y) / sqrt(2.0);
    const double diff_z = (boosted_proton.z - boosted_neutron.z) / sqrt(2.0);
    const double diff_px = (boosted_proton.px - boosted_neutron.px) / sqrt(2.0);
    const double diff_py = (boosted_proton.py - boosted_neutron.py) / sqrt(2.0);
    const double diff_pz = (boosted_proton.pz - boosted_neutron.pz) / sqrt(2.0);

    const double diff_r = sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    const double diff_p = sqrt(diff_px * diff_px + diff_py * diff_py + diff_pz * diff_pz);
    if (diff_p > config_input.cut_dp * 0.19733 / config_input.rms ||
        diff_r > config_input.cut_dr * config_input.rms) {
        return 0;
    }

    const double probility = calculate_pro(diff_r, diff_p, config_input.rms);
    const double pt = sqrt(px_total * px_total + py_total * py_total);
    if (const int npt = static_cast<int>(pt / config_input.d_mix_dpt); npt <= d_mix_spv.size()) {
        d_mix_spv[npt] += probility;
    }
    return probility;
}

// void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
//                                 const std::string &protonFileName,
//                                 const std::string &neutronFileName) {
//     std::ofstream protonFile(protonFileName, std::ios::out);
//     std::ofstream neutronFile(neutronFileName, std::ios::out);
//
//     for (auto &[eventID, particlesByType]: all_Events | std::views::values) {
//         protonFile << "t x y z px py pz p0\n";
//         neutronFile << "t x y z px py pz p0\n";
//
//         for (auto &[pdgCode, particles]: particlesByType) {
//             if (pdgCode == 2212 || pdgCode == 2112) {
//                 for (auto &particle: particles) {
//                     calculate_freeze_position(particle);
//                     std::ofstream &outputFile = (pdgCode == 2212) ? protonFile : neutronFile;
//
//                     outputFile << std::fixed << std::setprecision(7) << particle.freeze_out_time
//                                << " " << particle.x << " " << particle.y << " " << particle.z << " "
//                                << particle.px << " " << particle.py << " " << particle.pz << " "
//                                << particle.p0 << "\n";
//                 }
//             }
//         }
//     }
//
//     protonFile.close();
//     neutronFile.close();
// }
void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
                                const std::string &protonFileName,
                                const std::string &neutronFileName) {
    std::ofstream protonFile(protonFileName, std::ios::out);
    std::ofstream neutronFile(neutronFileName, std::ios::out);
    for (auto &[eventID, particlesByType] : all_Events | std::views::values) {
        protonFile << "t x y z px py pz p0\n";
        neutronFile << "t x y z px py pz p0\n";
        for (auto &[pdgCode, particles] : particlesByType) {
            if (pdgCode == 2212 || pdgCode == 2112) {
                for (auto &particle : particles) {
                    calculate_freeze_position(particle);
                    std::ofstream &outputFile =
                        (pdgCode == 2212) ? protonFile : neutronFile;
                    outputFile << std::fixed << std::setprecision(7)
                               << particle.freeze_out_time << " " << particle.x
                               << " " << particle.y << " " << particle.z << " "
                               << particle.px << " " << particle.py << " "
                               << particle.pz << " " << particle.p0 << "\n";
                }
            }
        }
    }
    protonFile.close();
    neutronFile.close();
}

void processBatches(const std::vector<ParticleData> &protons,
                    const std::vector<ParticleData> &neutrons,
                    const std::map<int, int> &protonCounts, const std::map<int, int> &neutronCounts,
                    const config_in &config_input, std::vector<double> &d_mix_spv) {
    int currentBatchSize = 0;
    double totalDeuterons = 0.0;
    int totalBatches = 0;
    int total_events = 0;
    int protonIndex = 0, neutronIndex = 0;
    int batchProtonCount = 0, batchNeutronCount = 0;

    std::fill(d_mix_spv.begin(), d_mix_spv.end(), 0.0);

    for (const auto &[eventID, pCount] : protonCounts) {
        const int nCount = neutronCounts.at(eventID);

        batchProtonCount += pCount;
        batchNeutronCount += nCount;
        currentBatchSize++;

        if (currentBatchSize >= config_input.mix_events ||
            std::next(protonCounts.find(eventID)) == protonCounts.end()) {
            double batchDeuterons = 0.0;
            int nparak = 0;

            for (int i = protonIndex; i < protonIndex + batchProtonCount; ++i) {
                const ParticleData &proton = protons[i];
                if (std::abs(calculate_rap(proton.p0, proton.pz)) > config_input.rap_cut_nucl) {
                    continue;
                }
                for (int j = neutronIndex; j < neutronIndex + batchNeutronCount; ++j) {
                    const ParticleData &neutron = neutrons[j];
                    if (std::abs(calculate_rap(neutron.p0, neutron.pz)) >
                        config_input.rap_cut_nucl) {
                        continue;
                    }

                    const double probility =
                        performCalculations(proton, neutron, config_input, d_mix_spv);
                    batchDeuterons += probility;
                    nparak++;
                }
            }

            protonIndex += batchProtonCount;
            neutronIndex += batchNeutronCount;

            std::cout << "Processing batch: " << totalBatches << ","
                      << "each batch events:" << config_input.mix_events << std::endl;
            std::cout << "each batch number of deuterons: " << batchDeuterons << std::endl;
            totalDeuterons += batchDeuterons;
            totalBatches++;
            total_events += nparak;
            batchProtonCount = 0;
            batchNeutronCount = 0;
            currentBatchSize = 0;
        }
    }
    std::cout << "total events :" << total_events << std::endl;
    const double averageDeuterons = totalBatches > 0 ? totalDeuterons / 100 : 0.0;
    std::cout << "Average number of deuterons per batch: " << averageDeuterons << std::endl;
}

void calcute_one_batch(const std::vector<ParticleData> &protons,
                       const std::vector<ParticleData> &neutrons, const config_in &config_input,
                       std::vector<double> &d_mix_spv, const std::vector<double> &d_mix_ptv,
                       double &batch_deutrons, const int eventsInBatch) {
    batch_deutrons = 0.0;
    const int mixEvents = eventsInBatch * eventsInBatch;

    for (const auto &proton : protons) {
        for (const auto &neutron : neutrons) {
            const double probility = performCalculations(proton, neutron, config_input, d_mix_spv);
            batch_deutrons += probility / mixEvents;
        }
    }
    for (int k = 0; k < 100; ++k) {
        if (d_mix_ptv[k] > 0) {
            d_mix_spv[k] /= (2 * M_PI * d_mix_ptv[k] * config_input.d_mix_dpt * mixEvents);
        }
    }
}

void calculate_deuteron(const std::string &protonFile,
                        const std::string &neutronFile,
                        const config_in &configInput,
                        const int batchSize, std::vector<double> &dMixSpv,
                        const std::vector<double> &dMixPtv) {
    BatchMap protonBatches, neutronBatches;
    double total_deutrons = 0.0;
    int total_batches = 0;

    read_batch_size(protonFile, batchSize, protonBatches);
    read_batch_size(neutronFile, batchSize, neutronBatches);

    for (const auto &[batchNumber, batchData] : protonBatches) {
        const auto &protons = batchData.particles;
        const auto &neutrons = neutronBatches[batchNumber].particles;
        double batch_deutrons = 0.0;
        const int eventsInBatch = batchData.eventCount;

        calcute_one_batch(protons,
                          neutrons,
                          configInput,
                          dMixSpv,
                          dMixPtv,
                          batch_deutrons,
                          eventsInBatch);
        total_deutrons += batch_deutrons;
        total_batches++;
        std::cout << "Average number of deuteron per couple: " << batch_deutrons << std::endl;
    }

    const double average_deutrons = total_batches > 0 ? total_deutrons / total_batches : 0.0;
    std::cout << "Average number of deuteron per batch: " << average_deutrons << std::endl;

    std::ofstream output("tem/d_mix_spv.dat");
    if (output.is_open()) {
        for (int k = 0; k < dMixSpv.size(); ++k) {
            dMixSpv[k] /= total_batches;
            output << std::setw(15) << dMixPtv[k] << std::setw(15) << dMixSpv[k] << std::endl;
        }
    }
    output.close();
}

