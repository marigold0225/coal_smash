//
// Created by mafu on 12/15/2023.
//
#include "../include/particle.h"

void ParticleData::updatePosition(double t_max) {
    x += (t_max - freeze_out_time) * px / p0;
    y += (t_max - freeze_out_time) * py / p0;
    z += (t_max - freeze_out_time) * pz / p0;
}

ParticleData ParticleData::lorentzBoost(const double beta_x, const double beta_y,
                                        const double beta_z) const {
    ParticleData boost_p = *this;
    if (const double beta2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z; beta2 > 1.0e-5) {
        const double gamma   = 1.0 / sqrt(1.0 - beta2);
        const double x_lam00 = gamma;
        const double x_lam01 = -gamma * beta_x;
        const double x_lam02 = -gamma * beta_y;
        const double x_lam03 = -gamma * beta_z;
        const double x_lam11 = 1.0 + (gamma - 1.0) * beta_x * beta_x / beta2;
        const double x_lam22 = 1.0 + (gamma - 1.0) * beta_y * beta_y / beta2;
        const double x_lam33 = 1.0 + (gamma - 1.0) * beta_z * beta_z / beta2;
        const double x_lam12 = (gamma - 1.0) * beta_x * beta_y / beta2;
        const double x_lam13 = (gamma - 1.0) * beta_x * beta_z / beta2;
        const double x_lam23 = (gamma - 1.0) * beta_y * beta_z / beta2;
        const double new_t   = freeze_out_time * x_lam00 + x * x_lam01 + y * x_lam02 + z * x_lam03;
        const double new_x   = freeze_out_time * x_lam01 + x * x_lam11 + y * x_lam12 + z * x_lam13;
        const double new_y   = freeze_out_time * x_lam02 + x * x_lam12 + y * x_lam22 + z * x_lam23;
        const double new_z   = freeze_out_time * x_lam03 + x * x_lam13 + y * x_lam23 + z * x_lam33;
        const double new_p0  = p0 * x_lam00 + px * x_lam01 + py * x_lam02 + pz * x_lam03;
        const double new_px  = p0 * x_lam01 + px * x_lam11 + py * x_lam12 + pz * x_lam13;
        const double new_py  = p0 * x_lam02 + px * x_lam12 + py * x_lam22 + pz * x_lam23;
        const double new_pz  = p0 * x_lam03 + px * x_lam13 + py * x_lam23 + pz * x_lam33;
        boost_p.freeze_out_time = new_t;
        boost_p.x               = new_x;
        boost_p.y               = new_y;
        boost_p.z               = new_z;
        boost_p.p0              = new_p0;
        boost_p.px              = new_px;
        boost_p.py              = new_py;
        boost_p.pz              = new_pz;
    }
    return boost_p;
}

bool ParticleData::operator!=(const ParticleData &other) const {
    constexpr double tolerance = 1e-6;
    return std::abs(t - other.t) > tolerance || std::abs(x - other.x) > tolerance ||
           std::abs(y - other.y) > tolerance || std::abs(z - other.z) > tolerance ||
           std::abs(p0 - other.p0) > tolerance || std::abs(px - other.px) > tolerance ||
           std::abs(py - other.py) > tolerance || std::abs(pz - other.pz) > tolerance;
}
void ParticleData::getTwobodyData(const ParticleData &p1, const ParticleData &p2) {
    freeze_out_time = std::max(p1.freeze_out_time, p2.freeze_out_time);
    x               = (p1.x + p2.x) / 2.0;
    y               = (p1.y + p2.y) / 2.0;
    z               = (p1.z + p2.z) / 2.0;
    p0              = p1.p0 + p2.p0;
    px              = p1.px + p2.px;
    py              = p1.py + p2.py;
    pz              = p1.pz + p2.pz;
    mass            = std::sqrt(p0 * p0 - px * px - py * py - pz * pz);
    charge          = p1.charge + p2.charge;
    pdg             = 1000010020;
}
void ParticleData::getFourbodyData(const ParticleData &p1, const ParticleData &p2,
                                   const ParticleData &n1, const ParticleData &n2) {
    freeze_out_time = std::max(std::max(p1.freeze_out_time, p2.freeze_out_time),
                               std::max(n1.freeze_out_time, n2.freeze_out_time));
    x               = (p1.x + p2.x + n1.x + n2.x) / 4.0;
    y               = (p1.y + p2.y + n1.y + n2.y) / 4.0;
    z               = (p1.z + p2.z + n1.z + n2.z) / 4.0;
    p0              = p1.p0 + p2.p0 + n1.p0 + n2.p0;
    px              = p1.px + p2.px + n1.px + n2.px;
    py              = p1.py + p2.py + n1.py + n2.py;
    pz              = p1.pz + p2.pz + n1.pz + n2.pz;
    mass            = std::sqrt(p0 * p0 - px * px - py * py - pz * pz);
    charge          = p1.charge + p2.charge + n1.charge + n2.charge;
    pdg             = 1000020030;
}
double ParticleData::getRapidity() const {
    if (p0 <= pz) {
        return 0.0;
    }
    return 0.5 * std::log((p0 + pz) / (p0 - pz));
}

double ParticleData::getArtifactRapidity() const {
    double p_total = std::sqrt(px * px + py * py + pz * pz);
    if (p_total <= pz) {
        return 0.0;
    }
    return 0.5 * std::log((p_total + pz) / (p_total - pz));
}

//double ParticleData::get_pt() const { return std::sqrt(px * px + py * py); }

void ParticleData::getFreezeOutPosition() {

    if (freeze_out_time == 0.0) {
        return;
    } else {
        double vx = px / p0;
        double vy = py / p0;
        double vz = pz / p0;
        double dt = t - freeze_out_time;
        double dx = vx * dt;
        double dy = vy * dt;
        double dz = vz * dt;
        x -= dx;
        y -= dy;
        z -= dz;
    }
}

int EventData::countChargeParticles() const {
    int chargeParticleCount = 0;
    for (const auto &[type, particles]: particlesByType) {
        for (const auto &particle: particles) {
            if (std::abs(particle.getArtifactRapidity()) < 0.5 && particle.charge != 0 &&
                std::abs(particle.pdg) > 100) {
                chargeParticleCount++;
            }
        }
    }
    return chargeParticleCount;
}
