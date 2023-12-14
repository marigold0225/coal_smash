//
// Created by mafu on 12/13/2023.
//

#pragma once

#include "particle.h"
#include <fstream>
#include <vector>

typedef struct BatchData {
    std::vector<ParticleData> particles;
    int eventCount;
} BatchData;
using BatchMap = std::map<int, BatchData>;

void ReadLine(const std::string &line, EventData &currentEvent);

void readFile(const std::string &filename, std::map<int, EventData> &all_Events);

std::pair<std::map<int, int>, std::map<int, int>>
countParticlesInEvents(const std::map<int, EventData> &allEvents);

void readParticleData(const std::string &filename,
                      std::vector<ParticleData> &particles);

void loadParticleData(std::vector<ParticleData> &batchProtons,
                      std::vector<ParticleData> &batchNeutrons);

void read_batch_size(const std::string &filename, int batchSize,
                     BatchMap &batches);

