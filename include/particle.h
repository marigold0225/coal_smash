//
// Created by mafu on 12/13/2023.
//

#pragma once

#include <cmath>
#include <map>
#include <vector>

struct ParticleData {
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
};


struct EventData {
    int eventID;
    std::map<int, std::vector<ParticleData>> particlesByType;
};
