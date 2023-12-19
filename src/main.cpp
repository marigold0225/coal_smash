//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_alpha.h"
#include <filesystem>
#include <stdexcept>

int main() {
    //    try {
    //        const std::string input_config_filename   = "config.ini";
    //        const std::string input_particle_filename = "particle_lists.oscar";
    //        const std::string dataOutputDir           = "data";
    //        if (!fileExistsInCurrentDir(input_config_filename) ||
    //            !fileExistsInCurrentDir(input_particle_filename)) {
    //            throw std::runtime_error("Required files not found in the current directory.");
    //        }
    //        checkAndCreateDataOutputDir(dataOutputDir);
    //    } catch (const std::exception &e) {
    //        std::cout << "Error" << e.what() << std::endl;
    //        return 1;
    //    }

    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::string input_config_filename   = "input/config.ini";
    const std::string input_particle_filename = "data/1000/particle_lists.oscar";
    const std::string dataOutputDir           = "data/1000";
    if (!fileExistsInCurrentDir(input_config_filename) ||
        !fileExistsInCurrentDir(input_particle_filename)) {
        throw std::runtime_error("Required files not found in the current directory.");
    }

    const config_parser read_config(input_config_filename);
    config_in config;
    initialize_config_from_parser(config, read_config);

    std::vector<std::string> centralityLabels = {"0-10", "10-20", "20-40", "40-80"};
    if (!checkFileExists(dataOutputDir, centralityLabels, "proton") ||
        !checkFileExists(dataOutputDir, centralityLabels, "neutron")) {
        std::cout << "Particle data being generated..." << std::endl;
        processParticleData(input_particle_filename, dataOutputDir);
    } else {
        std::cout << "Calculations using existing data" << std::endl;
    }
    //Deutron calculations
    if (config.rac_deuteron) {
        for (const auto &label: centralityLabels) {
            std::string protonFileName   = constructFilename(dataOutputDir, "proton", label);
            std::string neutronFileName  = constructFilename(dataOutputDir, "neutron", label);
            std::string deuteronFileName = constructFilename(dataOutputDir, "deuteron", label);
            std::string ptFileName       = constructFilename(dataOutputDir, "d_mix_spv", label);

            std::vector<double> deutron_pt(100, 0.0);
            calculate_deuteron(protonFileName, neutronFileName, deuteronFileName, ptFileName,
                               config, deutron_pt);
        }
    }
    //    //alpha calculations
    if (config.rac_helium4) {
        for (const auto &label: centralityLabels) {
            std::string protonFileName   = constructFilename(dataOutputDir, "proton", label);
            std::string neutronFileName  = constructFilename(dataOutputDir, "neutron", label);
            std::string deuteronFileName = constructFilename(dataOutputDir, "deuteron", label);
            std::string alphaFileName    = constructFilename(dataOutputDir, "alpha", label);
            std::string ptFileName       = constructFilename(dataOutputDir, "alpha_pt", label);

            std::vector<double> alpha_pt(100, 0.0);
            calculate_alpha_fourbody(protonFileName, neutronFileName, alphaFileName, config,
                                     alpha_pt);
            //            calculate_alpha_twobody(deuteronFileName, alphaFileName, ptFileName, config,
            //                                    alpha_pt);
        }
    }

    return 0;
}
