#pragma once

#include "general_function.h"
#include "read_particle_data.h"


double calculateDeutronProbability(double diff_r, double diff_p, double rms);


ParticleData pNToDeutron(const ParticleData &p1, const ParticleData &p2,
                         const config_in &config_input,
                         std::map<std::string, std::vector<double>> &pt_array,
                         const std::map<std::string, RapidityRange> &rapidityRange);

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const config_in &config_input,
                     std::map<std::string, std::vector<double>> &pt_array,
                     const std::map<std::string, RapidityRange> &rapidityRange,
                     double &batch_deutrons, int eventsInBatch,
                     std::vector<ParticleData> &deutrons);

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput,
                      std::map<std::string, std::vector<double>> &pt_array,
                      std::map<std::string, RapidityRange> &rapidityRanges);
