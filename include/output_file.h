//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "particle.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

void output_deutrons(const std::vector<ParticleData> &deutrons, const std::string &filename);


//void output_nuclei(const std::vector<ParticleData> &nuclei, const std::string &filename);

void output_d_mix_spv(const std::vector<double> &d_mix_spv,
                      const std::vector<double> &d_mix_ptv,
                      const std::string &filename,
                      int total_batch);
