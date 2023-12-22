//
// Created by mafu on 12/13/2023.
//

#include "../include/read_particle_data.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <ranges>
#include <sstream>

bool checkFileExists(const std::string &filename, std::vector<std::string> &labels,
                     const std::string &fileType) {
    for (const auto &label: labels) {
        std::string fileName = filename;
        fileName.append("/")
                .append(label)
                .append("/")
                .append(fileType)
                .append("_")
                .append(label)
                .append(".dat");
        if (!std::filesystem::exists(fileName)) {
            return false;
        }
    }
    return true;
}

void ReadLine(const std::string &line, EventData &currentEvent) {
    std::istringstream stream(line);
    ParticleData particle{};
    double n_coll = NAN, form_time = NAN, xsec_fac = NAN;
    int id = 0, proc_id_origin = 0, proc_type_origin = 0;
    int pdg_mother1   = 0;
    int pdg_mother2   = 0;
    int baryon_number = 0;

    stream >> particle.t >> particle.x >> particle.y >> particle.z >> particle.mass >>
            particle.p0 >> particle.px >> particle.py >> particle.pz >> particle.pdg >> id >>
            particle.charge >> n_coll >> form_time >> xsec_fac >> proc_id_origin >>
            proc_type_origin >> particle.freeze_out_time >> pdg_mother1 >> pdg_mother2 >>
            baryon_number;

    if (!stream.fail()) {
        currentEvent.particlesByType[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}

void readFile_smash(const std::string &filename, std::map<int, EventData> &all_Events) {
    std::ifstream file(filename);
    std::string line;
    EventData currentEvent;
    bool eventStarted = false;
    bool inEvent      = false;
    bool eventValid   = false;

    while (std::getline(file, line)) {
        if (line.starts_with("# event")) {
            if (line.find("out") != std::string::npos) {
                if (eventStarted && eventValid) {
                    all_Events[currentEvent.eventID] = currentEvent;
                    currentEvent                     = EventData{};
                }
                std::istringstream stream(line);
                std::string temp;
                stream >> temp >> temp >> currentEvent.eventID;
                eventStarted = true;
                inEvent      = true;
                eventValid   = false;
            } else if (line.find("end") != std::string::npos && inEvent) {
                if (line.find("scattering_projectile_target yes") != std::string::npos) {
                    eventValid = true;
                }
                inEvent = false;
            }
        } else if (inEvent) {
            ReadLine(line, currentEvent);
        }
    }

    if (eventStarted && eventValid) {
        all_Events[currentEvent.eventID] = currentEvent;
    }
}

void read_batch_nuclei(const std::string &filename, int batchSize, BatchMap &batches) {
    std::ifstream file(filename);
    std::string line;
    ParticleData particle{};
    int currentBatch = 0;
    int eventCount   = 0;

    while (std::getline(file, line)) {
        if (line.find("t x y z px py pz p0 mass t_out") != std::string::npos) {
            eventCount++;
            if (eventCount > batchSize) {
                batches[currentBatch].eventCount = eventCount - 1;
                currentBatch++;
                eventCount = 1;
            }
            continue;
        }

        if (std::istringstream iss(line); iss >> particle.t >> particle.x >> particle.y >>
                                          particle.z >> particle.px >> particle.py >> particle.pz >>
                                          particle.p0 >> particle.mass >>
                                          particle.freeze_out_time) {

            if (filename.find("proton") != std::string::npos) {
                particle.pdg    = 2212;
                particle.charge = 1;
            } else if (filename.find("neutron") != std::string::npos) {
                particle.pdg    = 2112;
                particle.charge = 0;
            }
            batches[currentBatch].particles.push_back(particle);
        }
    }
    if (eventCount > 0) {
        batches[currentBatch].eventCount = eventCount;
    }
}

void writeParticlesNoCentrality(std::map<int, EventData> &all_Events,
                                const std::string &protonFileName,
                                const std::string &neutronFileName) {
    std::ofstream protonFile(protonFileName, std::ios::out);
    std::ofstream neutronFile(neutronFileName, std::ios::out);
    for (auto &[eventID, particlesByType]: all_Events | std::views::values) {
        protonFile << "t x y z px py pz p0 mass t_out\n";
        neutronFile << "t x y z px py pz p0 mass t_out\n";
        for (auto &[pdgCode, particles]: particlesByType) {
            if (pdgCode == 2212 || pdgCode == 2112) {
                for (auto &particle: particles) {
                    particle.get_freeze_out_position();
                    std::ofstream &outputFile = (pdgCode == 2212) ? protonFile : neutronFile;
                    outputFile << std::fixed << particle.t << " " << particle.x << " " << particle.y
                               << " " << particle.z << " " << particle.px << " " << particle.py
                               << " " << particle.pz << " " << particle.p0 << " " << particle.mass
                               << " " << particle.freeze_out_time << "\n";
                }
            }
        }
    }
    protonFile.close();
    neutronFile.close();
}
void read_batch_deutrons(const std::string &filename, std::vector<BatchData> &batches) {
    std::ifstream file(filename);
    std::string line;
    ParticleData particle;
    int currentBatchNumber = 0;
    BatchData currentBatch;

    while (std::getline(file, line)) {
        if (line.find("t x y z px py pz p0 mass probability") != std::string::npos) {
            if (!currentBatch.particles.empty()) {
                currentBatch.eventCount = currentBatchNumber;
                batches.push_back(currentBatch);
                currentBatch.particles.clear();
                currentBatchNumber++;
            }
            continue;
        }

        if (std::istringstream iss(line); iss >> particle.freeze_out_time >> particle.x >>
                                          particle.y >> particle.z >> particle.px >> particle.py >>
                                          particle.pz >> particle.p0 >> particle.mass >>
                                          particle.probability) {
            currentBatch.particles.push_back(particle);
        }
    }
    if (!currentBatch.particles.empty()) {
        currentBatch.eventCount = currentBatchNumber;
        batches.push_back(currentBatch);
    }
}

std::string constructFilename(const std::string &outputDir, const std::string &fileType,
                              const std::string &label) {
    return outputDir + "/" + label + "/" + fileType + "_" + label + ".dat";
}

std::map<int, int> calculateMultiplicity(const std::map<int, EventData> &all_Events) {
    std::map<int, int> multiplicities;
    for (const auto &[eventID, eventData]: all_Events) {
        int multiplicity        = eventData.countChargeParticles();
        multiplicities[eventID] = multiplicity;
    }
    return multiplicities;
}

double percentile(const std::vector<int> &data, double percent) {
    if (data.empty()) {
        return 0.0;
    }
    if (percent <= 0) return data.front();
    if (percent >= 100) return data.back();

    std::vector<int> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    double idx      = static_cast<double>((data.size() - 1)) * percent / 100.0;
    int idx_lo      = static_cast<int>(idx);
    int idx_hi      = idx_lo + 1;
    double fraction = idx - idx_lo;

    if (idx_hi >= sortedData.size()) {
        return sortedData[idx_lo];
    }
    return sortedData[idx_lo] * (1.0 - fraction) + sortedData[idx_hi] * fraction;
}

void calculateCentralityBounds(const std::map<int, int> &multiplicities,
                               std::map<std::string, int> &centralityBounds) {
    std::vector<int> multiplicitiesVector;
    for (const auto &[eventID, multiplicity]: multiplicities) {
        multiplicitiesVector.push_back(multiplicity);
    }
    centralityBounds["40-80"] = static_cast<int>(percentile(multiplicitiesVector, 20));
    centralityBounds["20-40"] = static_cast<int>(percentile(multiplicitiesVector, 60));
    centralityBounds["10-20"] = static_cast<int>(percentile(multiplicitiesVector, 80));
    centralityBounds["0-10"]  = static_cast<int>(percentile(multiplicitiesVector, 90));
}

void classifyAndCountEvents(const std::map<int, int> &multiplicities,
                            const std::map<std::string, int> &centralityBounds,
                            std::map<int, std::string> &eventCentrality,
                            std::map<std::string, int> &centralityEventCounts) {
    for (const auto &[eventID, multiplicity]: multiplicities) {
        bool classified = false;
        for (const auto &[label, bound]: centralityBounds) {
            if (multiplicity >= bound) {
                eventCentrality[eventID] = label;
                centralityEventCounts[label]++;
                classified = true;
                break;
            }
        }
        if (!classified) {
            eventCentrality[eventID] = "unknown";
            centralityEventCounts["unknown"]++;
        }
    }
}

void writeParticlesByCentrality(std::map<int, EventData> &allEvents,
                                const std::map<int, std::string> &eventCentrality,
                                const std::string &outputDir) {
    // Define centrality labels
    std::vector<std::string> centralityLabels = {"0-10", "10-20", "20-40", "40-80"};
    std::map<std::string, std::ofstream> protonFiles, neutronFiles;
    for (const auto &label: centralityLabels) {
        std::string protonFileName = outputDir;
        protonFileName.append("/").append(label).append("/proton_").append(label).append(".dat");
        std::string neutronFileName = outputDir;
        neutronFileName.append("/").append(label).append("/neutron_").append(label).append(".dat");
        protonFiles[label].open(protonFileName, std::ios::out);
        neutronFiles[label].open(neutronFileName, std::ios::out);
    }
    // Process and write particle data
    for (auto &[eventID, eventData]: allEvents) {
        const std::string &centralityLabel = eventCentrality.at(eventID);
        std::ofstream &protonFile          = protonFiles[centralityLabel];
        std::ofstream &neutronFile         = neutronFiles[centralityLabel];

        std::string header = "t x y z px py pz p0 mass t_out\n";
        for (auto &[pdgCode, particles]: eventData.particlesByType) {
            if (pdgCode == 2212 || pdgCode == 2112) {
                std::ofstream &outputFile = (pdgCode == 2212) ? protonFile : neutronFile;
                outputFile << header;
                for (auto &particle: particles) {
                    particle.get_freeze_out_position();
                    outputFile << std::fixed << particle.t << " " << particle.x << " " << particle.y
                               << " " << particle.z << " " << particle.px << " " << particle.py
                               << " " << particle.pz << " " << particle.p0 << " " << particle.mass
                               << " " << particle.freeze_out_time << "\n";
                }
            }
        }
    }
    for (auto &[label, file]: protonFiles) file.close();
    for (auto &[label, file]: neutronFiles) file.close();
}

void processParticleData(const std::string &particle_file, const std::string &outputDir) {
    std::map<int, EventData> allEvents;
    readFile_smash(particle_file, allEvents);

    std::cout << "number of events: " << allEvents.size() << std::endl;

    std::map<int, int> multiplicity = calculateMultiplicity(allEvents);

    std::map<std::string, int> centralityBounds;
    calculateCentralityBounds(multiplicity, centralityBounds);

    std::map<int, std::string> eventCentrality;
    std::map<std::string, int> centralityEventCounts;
    classifyAndCountEvents(multiplicity, centralityBounds, eventCentrality, centralityEventCounts);
    for (const auto &[label, count]: centralityEventCounts) {
        std::cout << "Number of events in centrality range " << label << ": " << count << std::endl;
        std::string centralityDir = outputDir;
        centralityDir.append("/").append(label);
        checkAndCreateDataOutputDir(centralityDir);
    }
    writeParticlesByCentrality(allEvents, eventCentrality, outputDir);
}
void checkAndCreateDataOutputDir(const std::string &outputDir) {
    if (!std::filesystem::exists(outputDir)) {
        std::filesystem::create_directory(outputDir);
    }
}
bool fileExistsInCurrentDir(const std::string &filename) {
    return std::filesystem::exists(std::filesystem::current_path() / filename);
}
