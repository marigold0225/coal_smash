//
// Created by mafu on 12/3/2023.
//
#include "../include/read_config.h"

void initialize_config_from_parser(config_in &config1,
                                   const config_parser &config) {
    config1.mix_events   = config.get_int("MixEvents");
    config1.cut_dr       = config.get_double("cut_dr");
    config1.cut_dp       = config.get_double("cut_dp");
    config1.rap_cut_nucl = config.get_double("Rap_cut_nucl");
    config1.rap_cut_coal = config.get_double("Rap_cut_coal");
    config1.rms          = config.get_double("Rms");
    config1.d_mix_dpt    = config.get_double("d_mix_dpt");
    config1.rac_deuteron = config.get_bool("Deuteron");
    config1.rac_tritium  = config.get_bool("Tritium");
    config1.rac_helium   = config.get_bool("Helium");
}
