//
// Created by mafu on 12/15/2023.
//
#pragma once
#include "coal_deutron.h"


double calculate_fraction_batch_factor(const std::vector<ParticleData> &protons,
                                       const std::vector<ParticleData> &neutrons,
                                       std::vector<ParticleData> &protons_fraction,
                                       std::vector<ParticleData> &neutrons_fraction,
                                       const config_in &config_input);

ParticleData calculate_info_alpha_twobody(const ParticleData &d1, const ParticleData &d2,
                                          const config_in &config_input,
                                          std::map<std::string, std::vector<double>> &pt_array,
                                          std::map<std::string, RapidityRange> &rapidityRange);

ParticleData calculate_info_alpha_fourbody(const ParticleData &p1, const ParticleData &p2,
                                           const ParticleData &n1, const ParticleData &n2,
                                           const config_in &config_input,
                                           std::map<std::string, std::vector<double>> &pt_array,
                                           std::map<std::string, RapidityRange> &rapidityRangv);

std::tuple<double, double, double, double, double, double>
calculate_diff_fourbody(const ParticleData &p1, const ParticleData &p2, const ParticleData &n1,
                        const ParticleData &n2);

void calculate_alpha_twobody_batch(const std::vector<ParticleData> &deutrons,
                                   const config_in &config_input,
                                   std::map<std::string, std::vector<double>> &pt_array,
                                   std::map<std::string, RapidityRange> &rapidityRange,
                                   double &batch_alpha, int eventsInBatch,
                                   std::vector<ParticleData> &alpha);

void calculate_alpha_twobody(const std::string &deuteronFile, const std::string &alphaFile,
                             std::string &momentumFile, const config_in &configInput,
                             std::map<std::string, std::vector<double>> &pt_array,
                             std::map<std::string, RapidityRange> &rapidityRange);

void calculate_alpha_fourbody_batch(const std::vector<ParticleData> &protons,
                                    const std::vector<ParticleData> &neutrons,
                                    const config_in &config_input,
                                    std::map<std::string, std::vector<double>> &pt_array,
                                    std::map<std::string, RapidityRange> &rapidityRang,
                                    double &batch_alpha, int eventsInBatch,
                                    std::vector<ParticleData> &alpha);


void calculate_alpha_fourbody(const std::string &protonFile, const std::string &neutronFile,
                              const std::string &alphaFile, std::string &ptFile,
                              const config_in &configInput,
                              std::map<std::string, std::vector<double>> &pt_array,
                              std::map<std::string, RapidityRange> &rapidityRang);