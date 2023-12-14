//
// Created by mafu on 12/13/2023.
//

#include "../include/read_particle_data.h"
#include <cmath>
#include <iostream>
#include <sstream>

void ReadLine(const std::string &line, EventData &currentEvent) {
    std::istringstream stream(line);
    ParticleData particle{};
    double ncoll = NAN, form_time = NAN, xsecfac = NAN;
    int id = 0, proc_id_origin = 0, proc_type_origin = 0;
    int pdg_mother1 = 0;
    int pdg_mother2 = 0;
    int baryon_number = 0;

    stream >> particle.t >> particle.x >> particle.y >> particle.z >> particle.mass >>
           particle.p0 >> particle.px >> particle.py >> particle.pz >> particle.pdg >> id >>
           particle.charge >> ncoll >> form_time >> xsecfac >> proc_id_origin >> proc_type_origin >>
           particle.freeze_out_time >> pdg_mother1 >> pdg_mother2 >> baryon_number;

    if (!stream.fail()) {
        currentEvent.particlesByType[particle.pdg].push_back(particle);
    } else {
        std::cerr << "Failed to parse line: " << line << std::endl;
    }
}

void readFile(const std::string &filename, std::map<int, EventData> &all_Events) {
    std::ifstream file(filename);
    std::string line;
    EventData currentEvent;
    bool eventStarted = false;
    bool inEvent = false;

    while (std::getline(file, line)) {
        if (line.starts_with("# event")) {
            if (line.find("out") != std::string::npos) {
                if (eventStarted) {
                    all_Events[currentEvent.eventID] = currentEvent;
                    currentEvent = EventData();
                }
                std::istringstream stream(line);
                std::string temp;
                stream >> temp >> temp >> currentEvent.eventID;
                eventStarted = true;
                inEvent = true;
            } else if (line.find("end") != std::string::npos && inEvent) {
                inEvent = false;
            }
        } else if (inEvent) {
            ReadLine(line, currentEvent);
        }
    }

    if (eventStarted) {
        all_Events[currentEvent.eventID] = currentEvent;
    }
}

// void readFile(const std::string &filename) {
//     std::ifstream file(filename);
//     std::string line;
//     EventData currentEvent;
//     bool eventStarted = false;
//     bool inEvent      = false;
//
//     while (std::getline(file, line)) {
//         if (line.starts_with("# event")) {
//             if (line.find("out") != std::string::npos) {
//                 if (eventStarted) {
//                     allEvents[currentEvent.eventID] = currentEvent;
//                     currentEvent                    = EventData();
//                 }
//                 std::istringstream stream(line);
//                 std::string temp;
//                 stream >> temp >> temp >> currentEvent.eventID;
//                 eventStarted = true;
//                 inEvent      = true;
//             } else if (line.find("end") != std::string::npos && inEvent) {
//                 inEvent = false;
//             }
//         } else if (inEvent) {
//             ReadLine(line, currentEvent);
//         }
//     }
//
//     if (eventStarted) {
//         allEvents[currentEvent.eventID] = currentEvent;
//     }
// }
std::pair<std::map<int, int>, std::map<int, int>>
countParticlesInEvents(const std::map<int, EventData> &allEvents) {
    std::map<int, int> protonCountPerEvent;
    std::map<int, int> neutronCountPerEvent;

    for (const auto &[eventID, eventData] : allEvents) {
        for (const auto &[pdgCode, particles] : eventData.particlesByType) {
            const auto size = static_cast<int>(particles.size());
            if (size < 0) {
                throw std::runtime_error("Negative particle count");
            }
            if (pdgCode == 2212) {
                protonCountPerEvent[eventID] += size;
            } else if (pdgCode == 2112) {
                neutronCountPerEvent[eventID] += size;
            }
        }
    }

    return {protonCountPerEvent, neutronCountPerEvent};
}

void readParticleData(const std::string &filename, std::vector<ParticleData> &particles) {
    std::ifstream file(filename);
    std::string line;
    ParticleData particle{};

    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.find("t x y z px py pz p0") != std::string::npos) {
            continue;
        }

        if (std::istringstream iss(line); iss >> particle.freeze_out_time >> particle.x >>
                                              particle.y >> particle.z >> particle.px >> particle.py
                                              >>
                                              particle.pz >> particle.p0) {
            particle.t = 50;
            particle.mass = 0.938;
            particle.pdg = filename == "tem/proton.dat" ? 2212 : 2112;
            particle.charge = filename == "tem/proton.dat" ? 1 : 0;
            particles.push_back(particle);
        }
    }

    file.close();
}

void loadParticleData(std::vector<ParticleData> &batchProtons,
                      std::vector<ParticleData> &batchNeutrons) {
    readParticleData("tem/proton.dat", batchProtons);
    readParticleData("tem/neutron.dat", batchNeutrons);
}

void read_batch_size(const std::string &filename, int batchSize, BatchMap &batches) {
    std::ifstream file(filename);
    std::string line;
    ParticleData particle{};
    int currentBatch = 0;
    int eventCount = 0;

    while (std::getline(file, line)) {
        if (line.find("t x y z px py pz p0") != std::string::npos) {
            eventCount++;
            if (eventCount > batchSize) {
                batches[currentBatch].eventCount = eventCount - 1;
                currentBatch++;
                eventCount = 1;
            }
            continue;
        }

        if (std::istringstream iss(line); iss >> particle.freeze_out_time >> particle.x >>
                                              particle.y >> particle.z >> particle.px >> particle.py
                                              >>
                                              particle.pz >> particle.p0) {

            particle.t = 50;
            particle.mass = 0.938;
            particle.pdg = filename.find("tem/proton.dat") != std::string::npos ? 2212 : 2112;
            particle.charge = filename.find("tem/proton.dat") != std::string::npos ? 1 : 0;
            batches[currentBatch].particles.push_back(particle);
        }
    }
    if (eventCount > 0) {
        batches[currentBatch].eventCount = eventCount;
    }
}
