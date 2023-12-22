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

void checkAndCreateDataOutputDir(const std::string &outputDir);

bool fileExistsInCurrentDir(const std::string &filename);

void ReadLine(const std::string &line, EventData &currentEvent);

// handle datafile with centrality information
std::string constructFilename(const std::string &outputDir, const std::string &fileType,
                              const std::string &label);

std::map<int, int> calculateMultiplicity(const std::map<int, EventData> &all_Events);

void calculateCentralityBounds(const std::map<int, int> &multiplicities,
                               std::map<std::string, int> &centralityBounds);

void classifyAndCountEvents(const std::map<int, int> &multiplicities,
                            const std::map<std::string, int> &centralityBounds,
                            std::map<int, std::string> &eventCentrality,
                            std::map<std::string, int> &centralityEventCounts);

double percentile(const std::vector<int> &data, double percent);

void writeParticlesByCentrality(std::map<int, EventData> &allEvents,
                                const std::map<int, std::string> &eventCentrality,
                                const std::string &outputDir);

void processParticleData(const std::string &input_particle_filename, const std::string &outputDir);

void readFileSmash(const std::string &filename, std::map<int, EventData> &all_Events);

bool checkFileExists(const std::string &filename, std::vector<std::string> &labels,
                     const std::string &fileType);

void readBatchNuclei(const std::string &filename, int batchSize, BatchMap &batch_nuclei);

void readBatchDeutrons(const std::string &filename, std::vector<BatchData> &batche_deutron);

void writeParticlesNoCentrality(std::map<int, EventData> &all_Events,
                                const std::string &protonFileName,
                                const std::string &neutronFileName);
