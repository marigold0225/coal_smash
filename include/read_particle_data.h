//
// Created by mafu on 12/13/2023.
//

#pragma once

#include "output_file.h"
#include "read_config.h"
#include <fstream>
#include <vector>

typedef struct BatchData {
    std::vector<ParticleData> particles;
    int eventCount;
} BatchData;
using BatchMap = std::map<int, BatchData>;

void ReadLine(const std::string &line, EventData &currentEvent);

void calculate_freeze_position(ParticleData &p);

void readFile_smash(const std::string &filename, std::map<int, EventData> &all_Events);

void read_batch_size(const std::string &filename, int batchSize,
                     BatchMap &batches);

void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
                                const std::string &protonFileName,
                                const std::string &neutronFileName);