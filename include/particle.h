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

    void updatePosition(double t_max);

    [[nodiscard]] ParticleData lorentzBoost(double beta_x, double beta_y, double beta_z) const;

    [[nodiscard]] double getRapidity() const;

    [[nodiscard]] double getArtifactRapidity() const;

    //Calculating Transverse Momentum
    //    [[nodiscard]] double get_pt() const;

    //get freeze out position
    void getFreezeOutPosition();

    void getTwobodyData(const ParticleData &p1, const ParticleData &p2);

    void getFourbodyData(const ParticleData &p1, const ParticleData &p2, const ParticleData &n1,
                         const ParticleData &n2);


} ParticleData;

typedef struct RapidityRange {
    double min, max;
} RapidityRange;

typedef struct EventData {
    int eventID;
    std::map<int, std::vector<ParticleData>> particlesByType;
    [[nodiscard]] int countChargeParticles() const;
} EventData;
