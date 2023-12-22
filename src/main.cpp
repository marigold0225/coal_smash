//
// Created by mafu on 12/13/2023.
//
#include "../include/coal_alpha.h"
#include <filesystem>
#include <stdexcept>

int main() {

    std::cout << "Current path is " << std::filesystem::current_path() << std::endl;
    const std::string input_config_filename   = "input/config.ini";
    const std::string input_particle_filename = "data/6/particle_lists.oscar";
    const std::string dataOutputDir           = "data/6";

    checkAndCreateDataOutputDir(dataOutputDir);
    if (!fileExistsInCurrentDir(input_config_filename) ||
        !fileExistsInCurrentDir(input_particle_filename)) {
        throw std::runtime_error("Required files not found in the current directory.");
    }

    const config_parser read_config(input_config_filename);
    config_in config;
    initialize_config_from_parser(config, read_config);
    std::map<std::string, RapidityRange> rapidityRanges = {
            {"-0.1<y<0.0", {-0.1, 0.0}},   {"-0.2<y<-0.1", {-0.2, -0.1}},
            {"-0.3<y<-0.2", {-0.3, -0.2}}, {"-0.4<y<-0.3", {-0.4, -0.3}},
            {"-0.5<y<-0.4", {-0.5, -0.4}}, {"-0.6<y<-0.5", {-0.6, -0.5}},
            {"-0.7<y<-0.6", {-0.7, -0.6}}, {"-0.8<y<-0.7", {-0.8, -0.7}},
            {"-0.9<y<-0.8", {-0.9, -0.8}}, {"-1.0<y<-0.9", {-1.0, -0.9}}};

    if (config.needCentrality) {
        std::cout << "Centrality calculations being performed..." << std::endl;
        std::vector<std::string> centralityLabels = {"0-10", "10-20", "20-40", "40-80"};
        if (!checkFileExists(dataOutputDir, centralityLabels, "proton") ||
            !checkFileExists(dataOutputDir, centralityLabels, "neutron")) {
            std::cout << "Particle data being generated..." << std::endl;
            processParticleData(input_particle_filename, dataOutputDir);
        } else {
            std::cout << "Calculations using existing data" << std::endl;
        }
        //calculate proton transverse momentum
        if (!checkFileExists(dataOutputDir, centralityLabels, "p_pt")) {
            std::cout << "Calculating proton transverse momentum..." << std::endl;
            for (const auto &label: centralityLabels) {
                std::string protonFileName = constructFilename(dataOutputDir, "proton", label);
                std::string ptFileName     = constructFilename(dataOutputDir, "p_pt", label);

                std::map<std::string, std::vector<double>> proton_pt;
                calculateProtonPt(protonFileName, ptFileName, proton_pt, rapidityRanges);
            }
        } else {
            std::cout << "proton transverse momentum already exists.." << std::endl;
            //            for (const auto &label: centralityLabels) {
            //                std::string protonFileName = constructFilename(dataOutputDir, "proton", label);
            //                std::string ptFileName     = constructFilename(dataOutputDir, "p_pt", label);
            //
            //                std::map<std::string, std::vector<double>> proton_pt;
            //                calculateProtonPt(protonFileName, ptFileName, proton_pt, rapidityRanges);
            //            }
        }

        //Deutron calculations
        if (config.rac_deuteron) {
            for (const auto &label: centralityLabels) {
                std::string protonFileName   = constructFilename(dataOutputDir, "proton", label);
                std::string neutronFileName  = constructFilename(dataOutputDir, "neutron", label);
                std::string deuteronFileName = constructFilename(dataOutputDir, "deuteron", label);
                std::string ptFileName       = constructFilename(dataOutputDir, "d_pt", label);
                std::map<std::string, std::vector<double>> deutron_pt;
                std::cout << "Calculating deutron for centrality " << label << std::endl;
                calculate_deuteron(protonFileName, neutronFileName, deuteronFileName, ptFileName,
                                   config, deutron_pt, rapidityRanges);
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

                std::map<std::string, std::vector<double>> alpha_pt;
                std::cout << "Calculating alpha for centrality " << label << std::endl;
                //                calculate_alpha_fourbody(protonFileName, neutronFileName, alphaFileName, ptFileName,
                //                                         config, alpha_pt, rapidityRanges);
                calculate_alpha_twobody(deuteronFileName, alphaFileName, ptFileName, config,
                                        alpha_pt, rapidityRanges);
            }
        }
    } else {
        std::cout << "Centrality calculations not being performed..." << std::endl;

        auto protonFileName  = dataOutputDir + "/proton_no_centrality.dat";
        auto neutronFileName = dataOutputDir + "/neutron_no_centrality.dat";
        if (!fileExistsInCurrentDir(protonFileName) || !fileExistsInCurrentDir(neutronFileName)) {
            std::cout << "Particle data being generated..." << std::endl;
            std::map<int, EventData> all_Events;
            readFile_smash(input_particle_filename, all_Events);
            writeParticlesNoCentrality(all_Events, protonFileName, neutronFileName);
        } else {
            std::cout << "Calculations using existing data" << std::endl;
        }

        if (config.rac_deuteron) {
            std::string deuteronFileName = dataOutputDir + "/deuteron_no_centrality.dat";
            std::string ptFileName       = dataOutputDir + "/d_mix_spv_no_centrality.dat";
            std::map<std::string, std::vector<double>> deutron_pt;
            std::cout << "Calculating deutron for no centrality..." << std::endl;
            calculate_deuteron(protonFileName, neutronFileName, deuteronFileName, ptFileName,
                               config, deutron_pt, rapidityRanges);
        }
        if (config.rac_helium4) {
            std::string alphaFileName = dataOutputDir + "/alpha_no_centrality.dat";
            std::string ptFileName    = dataOutputDir + "/alpha_pt_no_centrality.dat";
            std::map<std::string, std::vector<double>> alpha_pt;
            std::cout << "Calculating alpha for no centrality..." << std::endl;
            calculate_alpha_fourbody(protonFileName, neutronFileName, alphaFileName, ptFileName,
                                     config, alpha_pt, rapidityRanges);
        }
    }

    return 0;
}
