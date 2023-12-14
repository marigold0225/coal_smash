//
// Created by mafu on 12/3/2023.
//
#include "../include/read_config.h"

void initialize_config_from_parser(config_in &config,
                                   const config_parser &read_config_parser) {
    config.calculation_mode = read_config_parser.get("Mode");// "smash" or "urqmd"
    config.mix_events       = read_config_parser.get_int("MixEvents");
    config.cut_dr           = read_config_parser.get_double("cut_dr");
    config.cut_dp           = read_config_parser.get_double("cut_dp");
    config.rap_cut_nucl     = read_config_parser.get_double("Rap_cut_nucl");
    config.rap_cut_coal     = read_config_parser.get_double("Rap_cut_coal");
    config.rms              = read_config_parser.get_double("Rms");
    config.d_mix_dpt        = read_config_parser.get_double("d_mix_dpt");
    config.rac_deuteron     = read_config_parser.get_bool("Deuteron");
    config.rac_tritium      = read_config_parser.get_bool("Tritium");
    config.rac_helium       = read_config_parser.get_bool("Helium");
}
