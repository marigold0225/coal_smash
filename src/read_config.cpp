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
    config.d_mix_dpt      = read_config_parser.get_double("d_mix_dpt");
    config.d_rap_cut_nucl = read_config_parser.get_double("d_Rap_cut_nucl");
    config.d_rap_cut_coal = read_config_parser.get_double("d_Rap_cut_coal");
    config.rms            = read_config_parser.get_double("d_Rms");
    //alpha
    config.alpha_mix_dpt      = read_config_parser.get_double("alpha_mix_dpt");
    config.alpha_rap_cut_nucl = read_config_parser.get_double("alpha_Rap_cut_nucl");
    config.alpha_rap_cut_coal = read_config_parser.get_double("alpha_Rap_cut_coal");
    config.alpha_rms          = read_config_parser.get_double("alpha_Rms");
    //switches
    config.rac_deuteron = read_config_parser.get_bool("Deuteron");
    config.rac_tritium  = read_config_parser.get_bool("Tritium");
    config.rac_helium4  = read_config_parser.get_bool("Helium4");
}
