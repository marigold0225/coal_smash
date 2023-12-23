//
// Created by mafu on 12/13/2023.
//
#pragma once
#include "config_parser.h"
typedef struct reactionConfig {
    double dpt;
    int ptBins;
    double rap_cut_nucl;
    double rap_cut_coal;
    double rms;
} reactionConfig;

typedef struct coalSwitch {
    bool rac_deuteron;
    bool rac_tritium;
    bool rac_helium4;
    bool rac_Be;
} coalSwitch;

typedef struct config_input {
    std::string calculation_mode;
    int mix_events;
    bool needCentrality;
    //select particles
    bool set_on;
    int proton_limit;
    int neutron_limit;
    //constant
    double cut_dr;
    double cut_dp;
    //deutron
    reactionConfig deuteron;
    reactionConfig alpha;
    reactionConfig Be;
    coalSwitch coalescenceSwitch;

} config_in;
void initializeConfigFromParser(config_in &config, const config_parser &read_config_parser);

config_in readAndParseConfig(const std::string &filename);
