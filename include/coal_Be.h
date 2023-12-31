//
// Created by mafu on 12/23/2023.
//
#pragma once
#include "general_function.h"
#include "output_file.h"
#include "read_particle_data.h"

ParticleData he4He4ToBe(const ParticleData &he4_1, const ParticleData &he4_2,
                        const reactionConfig &BeConfig, ptArray &pt_array,
                        const RapidityMap &rapidityRange,
                        std::map<std::string, double> &clusterCountOneBatch);

void processBeOneBatch(const std::vector<ParticleData> &alpha, const reactionConfig &BeConfig,
                       ptArray &pt_array, const RapidityMap &rapidityRange, double &batch_be,
                       int mixEvents, std::deque<ParticleData> &be,
                       std::map<std::string, double> &clusterCountByRapidity, std::mt19937 &gen,
                       std::uniform_real_distribution<> &dis);

void calculateBeAllBatch(const std::string &alphaFile, const std::string &beFile,
                         const std::string &ptFile, const reactionConfig &BeConfig,
                         const RapidityMap &rapidityRange, std::mt19937 &gen,
                         std::uniform_real_distribution<> &dis);
