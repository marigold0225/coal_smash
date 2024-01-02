//
// Created by mafu on 12/23/2023.
//

#pragma once

#include "coal_Be.h"
#include "coal_alpha.h"
#include "coal_deutron.h"
#include "read_particle_data.h"

void handleCentralityCalculations(const std::string &input_particle_filename,
                                  const std::string &dataOutputDir, const config_in &config,
                                  std::mt19937 &gen, std::uniform_real_distribution<> &dis);

void handleNoCentralityCalculations(const std::string &input_particle_filename,
                                    const std::string &dataOutputDir, const config_in &config,
                                    std::mt19937 &gen, std::uniform_real_distribution<> &dis);
