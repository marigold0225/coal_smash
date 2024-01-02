//
// Created by mafu on 12/22/2023.
//
#pragma once
#include "particle.h"
#include "read_particle_data.h"
#include <deque>
#include <functional>
#include <random>
#include <string>

RapidityMap defineRapidityRange();

double samplingAndScalingFactor(const std::vector<ParticleData> &protons,
                                const std::vector<ParticleData> &neutrons,
                                std::vector<ParticleData> &protons_fraction,
                                std::vector<ParticleData> &neutrons_fraction,
                                const config_in &config_input);

double calculateParticleRapidity(double p0, double pz);

void updateMomentumArray(const ParticleData &particle, const reactionConfig &config,
                         ptArray &pt_array, const RapidityMap &rapidityRange,
                         std::map<std::string, double> &clusterCountByRapidity);

void weightedSampling(std::vector<ParticleData> &sample_particles,
                      std::vector<std::pair<ParticleData, double>> &potential_particles,
                      std::vector<double> &cumulative_probabilities, double number_of_particles);

void weightedSample(std::deque<ParticleData> &sample_particles,
                    std::deque<std::pair<ParticleData, double>> &potential_particles,
                    std::deque<double> &cumulative_probabilities, double number_of_particles,
                    std::mt19937 &gen, std::uniform_real_distribution<> &dis);

std::tuple<double, double, double, double, double, double> TwoBodyJacobi(const ParticleData &p1,
                                                                         const ParticleData &p2);

std::tuple<double, double, double, double, double, double> fourBodyJacobi(const ParticleData &p1,
                                                                          const ParticleData &p2,
                                                                          const ParticleData &n1,
                                                                          const ParticleData &n2);

void calculateProtonPt(const std::string &protonFileName, const std::string &ptFileName,
                       ptArray &protonPt, const RapidityMap &rapidityRanges);

void calculateClusterPt(const std::string &clusterFileName, const std::string &ptFileName,
                        const RapidityMap &rapidityRanges, const reactionConfig &clusterConfig);