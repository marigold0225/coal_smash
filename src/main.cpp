//
// Created by mafu on 12/13/2023.
//
#include "../include/coal.h"
#include <filesystem>

int main() {
    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::string config_file_name = "input/config.ini";
    const std::string particle_file = "data/particle_lists.oscar";
    const std::string protonFileName = "tem/proton.dat";
    const std::string neutronFileName = "tem/neutron.dat";

    const config_parser config(config_file_name);
    config_in config_input;
    initialize_config_from_parser(config_input, config);
    const int batchSize = config_input.mix_events;
    std::vector<double> d_mix_spv(100, 0.0);
    std::vector<double> d_mix_ptv(100);
    for (int i = 0; i < 100; ++i) {
        d_mix_ptv[i] = config_input.d_mix_dpt / 2 + i * config_input.d_mix_dpt;
    }

    const bool proton_exists = std::filesystem::exists(protonFileName);
    const bool neutron_exists = std::filesystem::exists(neutronFileName);

    if (proton_exists && neutron_exists) {
        std::cout << "Using existing 'proton.dat' and 'neutron.dat' for "
                     "calculations."
                  << std::endl;
        calculate_deuteron(protonFileName, neutronFileName, config_input,
                           batchSize, d_mix_spv, d_mix_ptv);

    } else {
        std::cout << "'proton.dat' or 'neutron.dat' not found, using "
                     "'particle_lists.oscar' for calculations."
                  << std::endl;

        std::map<int, EventData> allEvents;
        readFile(particle_file, allEvents);

        std::cout << "Number of events:" << allEvents.size() << std::endl;

        extractParticlesFromEvents(allEvents, protonFileName, neutronFileName);
        calculate_deuteron(protonFileName, neutronFileName, config_input,
                           batchSize, d_mix_spv, d_mix_ptv);
    }

    return 0;
}

