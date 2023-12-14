//
// Created by mafu on 12/13/2023.
//
#pragma once
#include "config_parser.h"
typedef struct config_input {
    std::string calculation_mode;
    int mix_events;
    double rap_cut_nucl;
    double rap_cut_coal;
    double cut_dr;
    double cut_dp;
    double d_mix_dpt;
    double rms;
    // reaction rates
    bool rac_deuteron;
    bool rac_tritium;
    bool rac_helium;
} config_in;
void initialize_config_from_parser(config_in &config,
                                   const config_parser &read_config_parser);
