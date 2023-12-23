//
// Created by mafu on 12/23/2023.
//
#include "../include/run.h"


void handleCentralityCalculations(const std::string &input_particle_filename,
                                  const std::string &dataOutputDir, const config_in &config) {

    std::cout << "Centrality calculations being performed..." << std::endl;

    RapidityMap rapidityRanges                = defineRapidityRange();
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
            ptArray proton_pt;
            calculateProtonPt(protonFileName, ptFileName, proton_pt, rapidityRanges);
        }
    } else {
        std::cout << "proton transverse momentum already exists.." << std::endl;
    }

    //Deutron calculations
    if (config.coalescenceSwitch.rac_deuteron) {
        for (const auto &label: centralityLabels) {
            std::string protonFileName   = constructFilename(dataOutputDir, "proton", label);
            std::string neutronFileName  = constructFilename(dataOutputDir, "neutron", label);
            std::string deuteronFileName = constructFilename(dataOutputDir, "deuteron", label);
            std::string ptFileName       = constructFilename(dataOutputDir, "d_pt", label);
            ptArray deutron_pt;
            std::cout << "Calculating deutron for centrality " << label << std::endl;
            DeuteronAllBatch(protonFileName, neutronFileName, deuteronFileName, ptFileName, config,
                             deutron_pt, rapidityRanges);
        }
    }
    //    //alpha calculations
    if (config.coalescenceSwitch.rac_helium4) {
        for (const auto &label: centralityLabels) {
            std::string protonFileName   = constructFilename(dataOutputDir, "proton", label);
            std::string neutronFileName  = constructFilename(dataOutputDir, "neutron", label);
            std::string deuteronFileName = constructFilename(dataOutputDir, "deuteron", label);
            std::string alphaFileName    = constructFilename(dataOutputDir, "alpha", label);
            std::string ptFileName       = constructFilename(dataOutputDir, "alpha_pt", label);

            ptArray alpha_pt;
            std::cout << "Calculating alpha for centrality " << label << std::endl;
            //                calculateAlphaAllBatch4(protonFileName, neutronFileName, alphaFileName, ptFileName,
            //                                         config, alpha_pt, rapidityRanges);
            calculateAlphaAllBatch2(deuteronFileName, alphaFileName, ptFileName, config.alpha,
                                    alpha_pt, rapidityRanges);
        }
    }
    if (config.coalescenceSwitch.rac_Be) {
        for (const auto &label: centralityLabels) {
            std::string alphaFileName = constructFilename(dataOutputDir, "alpha", label);
            std::string beFileName    = constructFilename(dataOutputDir, "Be", label);
            std::string ptFileName    = constructFilename(dataOutputDir, "Be_pt", label);
            ptArray be_pt;
            std::cout << "Calculating Be for centrality " << label << std::endl;
            calculateBeAllBatch(alphaFileName, beFileName, ptFileName, config.Be, be_pt,
                                rapidityRanges);
        }
    }
    std::cout << "Centrality calculations completed." << std::endl;
}
void handleNoCentralityCalculations(const std::string &input_particle_filename,
                                    const std::string &dataOutputDir, const config_in &config) {

    std::cout << "Centrality calculations not being performed..." << std::endl;
    auto rapidityRanges = defineRapidityRange();

    auto protonFileName  = dataOutputDir + "/proton.dat";
    auto neutronFileName = dataOutputDir + "/neutron.dat";
    if (!fileExistsInCurrentDir(protonFileName) || !fileExistsInCurrentDir(neutronFileName)) {
        std::cout << "Particle data being generated..." << std::endl;
        std::map<int, EventData> all_Events;
        readFileSmash(input_particle_filename, all_Events);
        writeParticlesNoCentrality(all_Events, protonFileName, neutronFileName);
    } else {
        std::cout << "Calculations using existing data" << std::endl;
    }

    //calculate proton transverse momentum
    std::string protonPtFileName = dataOutputDir + "/p_pt.dat";
    if (!fileExistsInCurrentDir(protonPtFileName)) {
        std::cout << "Calculating proton transverse momentum..." << std::endl;
        ptArray proton_pt;
        calculateProtonPt(protonFileName, protonPtFileName, proton_pt, rapidityRanges);
    } else {
        std::cout << "proton transverse momentum already exists.." << std::endl;
    }

    if (config.coalescenceSwitch.rac_deuteron) {
        std::string deuteronFileName = dataOutputDir + "/deutron.dat";
        std::string ptFileName       = dataOutputDir + "/d_pt.dat";
        ptArray deutron_pt;
        std::cout << "Calculating deutron for no centrality..." << std::endl;
        DeuteronAllBatch(protonFileName, neutronFileName, deuteronFileName, ptFileName, config,
                         deutron_pt, rapidityRanges);
    }
    if (config.coalescenceSwitch.rac_helium4) {
        std::string deuteronFileName = dataOutputDir + "/deuteron.dat";
        std::string alphaFileName    = dataOutputDir + "/alpha.dat";
        std::string ptFileName       = dataOutputDir + "/alpha_pt.dat";
        ptArray alpha_pt;
        std::cout << "Calculating alpha for no centrality..." << std::endl;
        calculateAlphaAllBatch2(deuteronFileName, alphaFileName, ptFileName, config.alpha, alpha_pt,
                                rapidityRanges);
    }
    if (config.coalescenceSwitch.rac_Be) {
        std::string alphaFileName = dataOutputDir + "/alpha.dat";
        std::string beFileName    = dataOutputDir + "/Be.dat";
        std::string ptFileName    = dataOutputDir + "/Be_pt.dat";
        ptArray be_pt;
        std::cout << "Calculating Be for no centrality..." << std::endl;
        calculateBeAllBatch(alphaFileName, beFileName, ptFileName, config.Be, be_pt,
                            rapidityRanges);
    }
    std::cout << "No centrality calculations completed." << std::endl;
}
