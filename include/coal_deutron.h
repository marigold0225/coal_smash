#pragma once

#include "general_function.h"
#include "read_particle_data.h"

ParticleData pNToDeutron(const ParticleData &p1, const ParticleData &p2,
                         const reactionConfig &deutronConfig, ptArray &pt_array,
                         const RapidityMap &rapidityRange,
                         std::map<std::string, double> &clusterCountByRapidity);

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const reactionConfig &deutronConfig,
                     ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_deutrons,
                     int eventsInBatch, std::deque<ParticleData> &deutrons,
                     std::map<std::string, double> &clusterCountByRapidity, std::mt19937 &gen,
                     std::uniform_real_distribution<> &dis);

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput, const RapidityMap &rapidityRanges,
                      std::mt19937 &gen, std::uniform_real_distribution<> &dis);

void deutronAllBatchRandom(const std::string &protonFile, const std::string &neutronFile,
                           const std::string &deutronFile, const std::string &ptFile,
                           const config_in &configInput, const RapidityMap &rapidityRanges,
                           std::mt19937 &gen, std::uniform_real_distribution<> &dis);
