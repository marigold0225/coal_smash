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

void readFile_smash(const std::string &filename, std::map<int, EventData> &all_Events);

void read_batch_size(const std::string &filename, int batchSize,
                     BatchMap &batches);
