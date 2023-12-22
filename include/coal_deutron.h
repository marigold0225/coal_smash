#pragma once

#include "read_particle_data.h"

double calculate_rap(double p0, double pz);

double calDeutronPro(double diff_r, double diff_p, double rms);

std::tuple<double, double, double, double, double, double>
calculateTwoBodyDifferences(const ParticleData &p1, const ParticleData &p2);

void update_momentum_array(double pt, double probability, double d_pt, double rap,
                           std::map<std::string, std::vector<double>> &pt_array,
                           const std::map<std::string, RapidityRange> &rapidityRange);

void calculateProtonPt(const std::string &protonFileName, const std::string &ptFileName,
                       std::map<std::string, std::vector<double>> &protonPtsByRapidity,
                       const std::map<std::string, RapidityRange> &rapidityRanges);

ParticleData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                    const config_in &config_input,
                                    std::map<std::string, std::vector<double>> &pt_array,
                                    const std::map<std::string, RapidityRange> &rapidityRange);

void calculate_deutron_batch(const std::vector<ParticleData> &protons,
                             const std::vector<ParticleData> &neutrons,
                             const config_in &config_input,
                             std::map<std::string, std::vector<double>> &pt_array,
                             const std::map<std::string, RapidityRange> &rapidityRange,
                             double &batch_deutrons, int eventsInBatch,
                             std::vector<ParticleData> &deutrons);

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const std::string &deutronFile, const std::string &ptFile,
                        const config_in &configInput,
                        std::map<std::string, std::vector<double>> &pt_array,
                        std::map<std::string, RapidityRange> &rapidityRanges);

void weighted_sampling(std::vector<ParticleData> &sample_particles,
                       std::vector<std::pair<ParticleData, double>> &potential_particles,
                       std::vector<double> &cumulative_probabilities, double number_of_particles);