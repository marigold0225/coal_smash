#pragma once

#include "read_particle_data.h"

double calculate_rap(double p0, double pz);

double calculate_pro_deutron(double diff_r, double diff_p, double rms);

std::tuple<double, double, double, double, double, double>
calculate_differences(const ParticleData &p1, const ParticleData &p2);

void update_d_mix_spv(double pt, double probability, const config_in &config_input,
                      std::vector<double> &d_mix_spv);

ParticleData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                    const config_in &config_input,
                                    std::vector<double> &d_mix_spv);

void calculate_one_batch(const std::vector<ParticleData> &protons, const std::vector<ParticleData> &neutrons,
                         const config_in &config_input, std::vector<double> &d_mix_spv,
                         const std::vector<double> &d_mix_ptv, double &batch_deutrons,
                         int eventsInBatch, std::vector<ParticleData> &deutrons);

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const config_in &configInput, int batchSize,
                        std::vector<double> &dMixSpv, const std::vector<double> &dMixPtv,
                        std::vector<ParticleData> &deutrons);
