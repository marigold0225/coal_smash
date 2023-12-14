#pragma once

#include "read_config.h"
#include "read_particle_data.h"

void calculate_freeze_position(ParticleData &p);

double calculate_rap(double p0, double pz);

ParticleData lorentz_boost(double beta_x, double beta_y, double beta_z,
                           const ParticleData &p);

double calculate_pro(double diff_r, double diff_p, double rms);

double performCalculations(const ParticleData &p1, const ParticleData &p2,
                           const config_in &config_input,
                           std::vector<double> &d_mix_spv);

void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
                                const std::string &protonFileName,
                                const std::string &neutronFileName);

void calculate_one_batch(const std::vector<ParticleData> &protons,
                         const std::vector<ParticleData> &neutrons, const config_in &config_input,
                         std::vector<double> &d_mix_spv, const std::vector<double> &d_mix_ptv,
                         double &batch_deutrons, int eventsInBatch);

void calculate_deuteron(const std::string &protonFile, const std::string &neutronFile,
                        const config_in &configInput, int batchSize,
                        std::vector<double> &dMixSpv, const std::vector<double> &dMixPtv);
