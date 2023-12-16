//
// Created by mafu on 12/13/2023.
//

#pragma once

#include <cmath>
#include <map>
#include <vector>

typedef struct ParticleData {
    double t, x, y, z;
    double p0, px, py, pz;
    double mass;
    int pdg, charge;
    double freeze_out_time;
    [[maybe_unused]] double probability;

    bool operator!=(const ParticleData &other) const;

    void update_position(double t_max);

    [[nodiscard]] ParticleData lorentz_boost(double beta_x, double beta_y, double beta_z) const;

    void calculate_deutron_data(const ParticleData &proton, const ParticleData &neutron);

    void calculate_fourbody_data(const ParticleData &p1, const ParticleData &p2,
                                 const ParticleData &n1, const ParticleData &n2);


} ParticleData;

typedef struct EventData {
    int eventID;
    std::map<int, std::vector<ParticleData>> particlesByType;
} EventData;
