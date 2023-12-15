//
// Created by mafu on 12/13/2023.
//

#include "../include/read_particle_data.h"
#include <cmath>
#include <iostream>
#include <ranges>
#include <sstream>

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
            particle.charge >> n_coll >> form_time >> xsec_fac >> proc_id_origin >> proc_type_origin >>
            particle.freeze_out_time >> pdg_mother1 >> pdg_mother2 >> baryon_number;

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

    while (std::getline(file, line)) {
        if (line.starts_with("# event")) {
            if (line.find("out") != std::string::npos) {
                if (eventStarted) {
                    all_Events[currentEvent.eventID] = currentEvent;
                    currentEvent                     = EventData();
                }
                std::istringstream stream(line);
                std::string temp;
                stream >> temp >> temp >> currentEvent.eventID;
                eventStarted = true;
                inEvent      = true;
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

void read_batch_size(const std::string &filename, int batchSize, BatchMap &batches) {
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

        if (std::istringstream iss(line); iss >> particle.t >> particle.x >>
                                          particle.y >> particle.z >> particle.px >> particle.py >>
                                          particle.pz >> particle.p0 >> particle.mass >> particle.freeze_out_time) {

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
void calculate_freeze_position(ParticleData &p) {
    const double vx = p.px / p.p0;
    const double vy = p.py / p.p0;
    const double vz = p.pz / p.p0;
    const double dt = p.t - p.freeze_out_time;
    const double dx = vx * dt;
    const double dy = vy * dt;
    const double dz = vz * dt;
    p.x -= dx;
    p.y -= dy;
    p.z -= dz;
}
void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
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
                    calculate_freeze_position(particle);
                    std::ofstream &outputFile =
                            (pdgCode == 2212) ? protonFile : neutronFile;
                    outputFile << std::fixed << std::setprecision(7)
                               << particle.t << " " << particle.x
                               << " " << particle.y << " " << particle.z << " "
                               << particle.px << " " << particle.py << " "
                               << particle.pz << " " << particle.p0 << " "
                               << particle.mass << " " << particle.freeze_out_time << "\n";
                }
            }
        }
    }
    protonFile.close();
    neutronFile.close();
}
//void extractParticlesFromEvents(std::map<int, EventData> &all_Events,
//                                const std::string &protonFileName,
//                                const std::string &neutronFileName) {
//    std::vector<ParticleData> protons;
//    std::vector<ParticleData> neutrons;
//    for (auto &[eventID, particlesByType]: all_Events | std::views::values) {
//        for (auto &[pdgCode, particles]: particlesByType) {
//            if (pdgCode == 2212 || pdgCode == 2112) {
//                for (auto &particle: particles) {
//                    calculate_freeze_position(particle);
//                    if (pdgCode == 2212) {
//                        protons.push_back(particle);
//                    } else {
//                        neutrons.push_back(particle);
//                    }
//                }
//            }
//        }
//    }
//    output_nuclei(protons, protonFileName);
//    output_nuclei(neutrons, neutronFileName);
//}