//
// Created by mafu on 12/13/2023.
//
#pragma once
#include "config_parser.h"
typedef struct config_input {
    std::string calculation_mode;
    int mix_events;
    //select particles
    bool set_on;
    int proton_limit;
    int neutron_limit;
    //constant
    double cut_dr;
    double cut_dp;
    //deutron
    double d_mix_dpt;
    double d_rap_cut_nucl;
    double d_rap_cut_coal;
    double rms;
    //alpha
    double alpha_mix_dpt;
    double alpha_rap_cut_nucl;
    double alpha_rap_cut_coal;
    double alpha_rms;
    // reaction rates
    bool rac_deuteron;
    bool rac_tritium;
    bool rac_helium4;
} config_in;
void initialize_config_from_parser(config_in &config,
                                   const config_parser &read_config_parser);
