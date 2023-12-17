//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_alpha.h"
#include "../include/coal_deutron.h"
#include <filesystem>

int main() {
    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::string config_file_name = "input/config.ini";
    const std::string particle_file    = "data/particle_lists.oscar";
    const std::string protonFileName   = "data/proton.dat";
    const std::string neutronFileName  = "data/neutron.dat";
    const std::string deuteronFileName = "data/deuteron.dat";
    const std::string alphaFileName    = "data/alpha_two.dat";
    const std::string alphaFileName4   = "data/alpha_four.dat";

    const config_parser config(config_file_name);
    config_in config_input;
    initialize_config_from_parser(config_input, config);
    const int batchSize = config_input.mix_events;

    const bool proton_exists  = std::filesystem::exists(protonFileName);
    const bool neutron_exists = std::filesystem::exists(neutronFileName);

    if (!proton_exists || !neutron_exists) {
        std::cout << "'proton.dat' or 'neutron.dat' not found, using "
                     "'particle_lists.oscar' for calculations."
                  << std::endl;
        std::map<int, EventData> allEvents;
        readFile_smash(particle_file, allEvents);
        std::cout << "number of events: " << allEvents.size() << std::endl;
        extractParticlesFromEvents(allEvents, protonFileName, neutronFileName);
    } else {
        std::cout << "Using existing 'proton.dat' and 'neutron.dat' for calculations."
                  << std::endl;
    }

    int steps = 100;
    //Deutron calculations
    if (config_input.rac_deuteron) {
        std::vector<double> deutron_mix_spv(steps, 0.0);
        std::vector<double> deutron_mix_ptv(steps);
        std::vector<ParticleData> deutrons;
        for (int i = 0; i < steps; ++i) {
            deutron_mix_ptv[i] = config_input.d_mix_dpt / 2 + static_cast<double>(i) * config_input.d_mix_dpt;
        }
        calculate_deuteron(protonFileName, neutronFileName, deuteronFileName,
                           config_input, batchSize,
                           deutron_mix_spv, deutron_mix_ptv, deutrons);
    }
    //alpha calculations
    if (config_input.rac_helium4) {
        std::vector<double> alpha_mix_spv(steps, 0.0);
        std::vector<double> alpha_mix_ptv(steps);
        std::vector<ParticleData> alpha;
        for (int i = 0; i < steps; ++i) {
            alpha_mix_ptv[i] = config_input.d_mix_dpt / 2 + i * config_input.d_mix_dpt;
        }
        //        calculate_alpha_fourbody(protonFileName, neutronFileName, alphaFileName4,
        //                                 config_input, batchSize,
        //                                 alpha_mix_spv, alpha_mix_ptv, alpha);
        calculate_alpha_twobody(deuteronFileName, alphaFileName,
                                config_input, batchSize,
                                alpha_mix_spv, alpha_mix_ptv, alpha);
    }
    return 0;
}
