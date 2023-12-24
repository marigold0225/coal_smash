//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "general_function.h"
#include "read_particle_data.h"

ParticleData dDToAlpha(const ParticleData &d1, const ParticleData &d2,
                       const reactionConfig &alphaConfig, ptArray &pt_array,
                       const RapidityMap &rapidityRange);

void processAlphaOneBatch2(const std::vector<ParticleData> &deutrons,
                           const reactionConfig &alphaConfig, ptArray &pt_array,
                           const RapidityMap &rapidityRange, double &batch_alpha, int mixEvents,
                           std::vector<ParticleData> &alpha);

void calculateAlphaAllBatch2(const std::string &deuteronFile, const std::string &alphaFile,
                             std::string &ptFile, const reactionConfig &alphaConfig,
                             const RapidityMap &rapidityRange);

ParticleData pPNNToAlpha(const ParticleData &p1, const ParticleData &p2, const ParticleData &n1,
                         const ParticleData &n2, const config_in &config_input, ptArray &pt_array,
                         const RapidityMap &rapidityRange);

void processAlphaOneBatch4(const std::vector<ParticleData> &protons,
                           const std::vector<ParticleData> &neutrons, const config_in &config_input,
                           ptArray &pt_array, const RapidityMap &rapidityRang, double &batch_alpha,
                           int eventsInBatch, std::vector<ParticleData> &alpha);

void calculateAlphaAllBatch4(const std::string &protonFile, const std::string &neutronFile,
                             const std::string &alphaFile, std::string &ptFile,
                             const config_in &configInput, const RapidityMap &rapidityRang);