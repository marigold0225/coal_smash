#pragma once

#include "general_function.h"
#include "read_particle_data.h"

ParticleData pNToDeutron(const ParticleData &p1, const ParticleData &p2,
                         const reactionConfig &deutronConfig, ptArray &pt_array,
                         const RapidityMap &rapidityRange);

void DeutronOneBatch(const std::vector<ParticleData> &protons,
                     const std::vector<ParticleData> &neutrons, const reactionConfig &deutronConfig,
                     ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_deutrons,
                     int eventsInBatch, std::vector<ParticleData> &deutrons);

void DeuteronAllBatch(const std::string &protonFile, const std::string &neutronFile,
                      const std::string &deutronFile, const std::string &ptFile,
                      const config_in &configInput, ptArray &pt_array,
                      const RapidityMap &rapidityRanges);
