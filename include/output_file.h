//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "particle.h"
#include "read_config.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

void outputCluster(const std::vector<ParticleData> &clusters, int events, std::ofstream &output);


void writeNucleiData(ParticleData &nuclei, std::ofstream &outputFile);

void outputPt(ptArray &pt_array, std::map<std::string, double> clusterCountByRapidity,
              const reactionConfig &clusterConfig, const std::string &filename, int total_batch);
