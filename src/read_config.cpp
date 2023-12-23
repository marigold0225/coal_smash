//
// Created by mafu on 12/3/2023.
//
#include "../include/read_config.h"

void initializeConfigFromParser(config_in &config, const config_parser &read_config_parser) {
    config.calculation_mode = read_config_parser.get("Mode");// "smash" or "urqmd"
    config.mix_events       = read_config_parser.get_int("MixEvents");
    config.needCentrality   = read_config_parser.get_bool("needCentrality");
    //select particles
    config.set_on        = read_config_parser.get_bool("switch");
    config.proton_limit  = read_config_parser.get_int("proton_limit");
    config.neutron_limit = read_config_parser.get_int("neutron_limit");
    //constants
    config.cut_dr = read_config_parser.get_double("cut_dr");
    config.cut_dp = read_config_parser.get_double("cut_dp");
    //deutron
    config.deuteron.dpt          = read_config_parser.get_double("d_mix_dpt");
    config.deuteron.ptBins       = read_config_parser.get_int("d_ptBins");
    config.deuteron.rap_cut_nucl = read_config_parser.get_double("d_Rap_cut_nucl");
    config.deuteron.rap_cut_coal = read_config_parser.get_double("d_Rap_cut_coal");
    config.deuteron.rms          = read_config_parser.get_double("d_Rms");
    //alpha
    config.alpha.dpt          = read_config_parser.get_double("alpha_mix_dpt");
    config.alpha.ptBins       = read_config_parser.get_int("alpha_ptBins");
    config.alpha.rap_cut_nucl = read_config_parser.get_double("alpha_Rap_cut_nucl");
    config.alpha.rap_cut_coal = read_config_parser.get_double("alpha_Rap_cut_coal");
    config.alpha.rms          = read_config_parser.get_double("alpha_Rms");
    //Be
    config.Be.dpt          = read_config_parser.get_double("Be_mix_dpt");
    config.Be.ptBins       = read_config_parser.get_int("Be_ptBins");
    config.Be.rap_cut_nucl = read_config_parser.get_double("Be_Rap_cut_nucl");
    config.Be.rap_cut_coal = read_config_parser.get_double("Be_Rap_cut_coal");
    config.Be.rms          = read_config_parser.get_double("Be_Rms");
    //switches
    config.coalescenceSwitch.rac_deuteron = read_config_parser.get_bool("Deuteron");
    config.coalescenceSwitch.rac_tritium  = read_config_parser.get_bool("Tritium");
    config.coalescenceSwitch.rac_helium4  = read_config_parser.get_bool("Helium4");
    config.coalescenceSwitch.rac_Be       = read_config_parser.get_bool("Be");
}
config_in readAndParseConfig(const std::string &filename) {
    const config_parser read_config_parser(filename);
    config_in config;
    initializeConfigFromParser(config, read_config_parser);
    return config;
}
