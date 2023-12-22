//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "particle.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

void outputCluster(const std::vector<ParticleData> &clusters, std::ofstream &output);


void writeNucleiData(ParticleData &nuclei, std::ofstream &outputFile);

void outputPt(std::map<std::string, std::vector<double>> &pt_array,
              std::map<std::string, double> clusterCountByRapidity, double d_pt, int ptBins,
              const std::string &filename, int total_batch);
