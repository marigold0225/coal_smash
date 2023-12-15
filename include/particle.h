//
// Created by mafu on 12/13/2023.
//

#pragma once

#include <cmath>
#include <map>
#include <vector>

typedef struct ParticleData {
    double t, x, y, z;
    double mass, p0, px, py, pz;
    int pdg, charge;
    double freeze_out_time;

    bool operator!=(const ParticleData &other) const {
        constexpr double tolerance = 1e-6;
        return std::abs(t - other.t) > tolerance ||
               std::abs(x - other.x) > tolerance ||
               std::abs(y - other.y) > tolerance ||
               std::abs(z - other.z) > tolerance ||
               std::abs(p0 - other.p0) > tolerance ||
               std::abs(px - other.px) > tolerance ||
               std::abs(py - other.py) > tolerance ||
               std::abs(pz - other.pz) > tolerance;
    }
    void update_position(double t_max) {
        x += (t_max - freeze_out_time) * px / p0;
        y += (t_max - freeze_out_time) * py / p0;
        z += (t_max - freeze_out_time) * pz / p0;
    }
} ParticleData;

typedef struct DeutronData {
    double t, x, y, z;
    double p0, px, py, pz;
    double probability;
    void calculate_deutron_data(const ParticleData &proton, const ParticleData &neutron) {
        t  = std::max(proton.freeze_out_time, neutron.freeze_out_time);
        x  = (proton.x + neutron.x) / 2.0;
        y  = (proton.y + neutron.y) / 2.0;
        z  = (proton.z + neutron.z) / 2.0;
        p0 = proton.p0 + neutron.p0;
        px = proton.px + neutron.px;
        py = proton.py + neutron.py;
        pz = proton.pz + neutron.pz;
    }
} DeutronData;

typedef struct EventData {
    int eventID;
    std::map<int, std::vector<ParticleData>> particlesByType;
} EventData;
