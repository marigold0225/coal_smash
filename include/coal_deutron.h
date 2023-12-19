#pragma once

#include "read_particle_data.h"

double calculate_rap(double p0, double pz);

double calculate_pro_deutron(double diff_r, double diff_p, double rms);

std::tuple<double, double, double, double, double, double>
calculateTwoBodyDifferences(const ParticleData &p1, const ParticleData &p2);

void update_momentum_array(double pt, double probability, const config_in &config_input,
                           std::vector<double> &pt_array);

ParticleData calculate_info_deutron(const ParticleData &p1, const ParticleData &p2,
                                    const config_in &config_input, std::vector<double> &pt_array);

void calculate_deutron_batch(const std::vector<ParticleData> &protons,
                             const std::vector<ParticleData> &neutrons,
                             const config_in &config_input, std::vector<double> &pt_array,
                             const std::vector<double> &d_pt, double &batch_deutrons,
                             int eventsInBatch, std::vector<ParticleData> &deutrons);

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const std::string &deutronFile, const std::string &ptFile,
                        const config_in &configInput, std::vector<double> &pt_array);

void weighted_sampling(std::vector<ParticleData> &sample_particles,
                       std::vector<std::pair<ParticleData, double>> &potential_particles,
                       std::vector<double> &cumulative_probabilities, double number_of_particles);