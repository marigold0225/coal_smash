//
// Created by mafu on 12/13/2023.
//
#include "../include/run.h"
#include <filesystem>
#include <stdexcept>

int main() {

    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::string input_config_filename   = "input/config.ini";
    const std::string input_particle_filename = "data/1000/particle_lists.oscar";
    const std::string dataOutputDir           = "data/1000";

    checkAndCreateDataOutputDir(dataOutputDir);
    if (!fileExistsInCurrentDir(input_config_filename) ||
        !fileExistsInCurrentDir(input_particle_filename)) {
        throw std::runtime_error("Required files not found in the current directory.");
    }

    config_in config = readAndParseConfig(input_config_filename);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    if (config.needCentrality) {
        handleCentralityCalculations(input_particle_filename, dataOutputDir, config, gen, dis);

    } else {
        handleNoCentralityCalculations(input_particle_filename, dataOutputDir, config, gen, dis);
    }

    return 0;
}
