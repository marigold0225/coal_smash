//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "particle.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

void output_cluster(const std::vector<ParticleData> &deutrons, std::ofstream &output);


//void output_nuclei(const std::vector<ParticleData> &nuclei, const std::string &filename);

void output_spv(std::map<std::string, std::vector<double>> &pt_array,
                std::map<std::string, double> clusterCountByRapidity, double d_pt, int ptBins,
                const std::string &filename, int total_batch);
