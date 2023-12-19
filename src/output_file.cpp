//
// Created by mafu on 12/15/2023.
//
#include "../include/output_file.h"

//void output_nuclei(const std::vector<ParticleData> &nuclei, const std::string &filename) {
//    std::ofstream output_file(filename, std::ios::out);
//    output_file << "t x y z px py pz p0 mass t_out\n";
//    for (const auto &nucleus: nuclei) {
//        output_file << nucleus.t << " " << nucleus.x << " " << nucleus.y << " " << nucleus.z << " "
//                    << nucleus.px << " " << nucleus.py << " " << nucleus.pz << " " << nucleus.p0 << " "
//                    << nucleus.mass << " " << nucleus.freeze_out_time << "\n";
//    }
//    output_file.close();
//}

void output_cluster(const std::vector<ParticleData> &clusters, std::ofstream &output) {
    output << "t x y z px py pz p0 mass probability\n";
    for (const auto &cluster: clusters) {
        output << std::fixed << std::setprecision(15) << cluster.freeze_out_time << " " << cluster.x
               << " " << cluster.y << " " << cluster.z << " " << cluster.px << " " << cluster.py
               << " " << cluster.pz << " " << cluster.p0 << " " << cluster.mass << " "
               << cluster.probability << "\n";
    }
}


void output_spv(std::vector<double> &d_mix_spv, const std::vector<double> &d_mix_ptv,
                const std::string &filename, int total_batch) {
    std::ofstream output_file(filename, std::ios::out);
    if (output_file.is_open()) {
        for (int i = 0; i < d_mix_spv.size(); ++i) {
            d_mix_spv[i] /= total_batch;
            output_file << std::setw(15) << d_mix_ptv[i] << std::setw(15) << d_mix_spv[i]
                        << std::endl;
        }
    }
    output_file.close();
}