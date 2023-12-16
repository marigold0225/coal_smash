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

void output_deutrons(const std::vector<ParticleData> &deutrons, std::ofstream &output) {
    output << "t x y z px py pz p0 mass probability\n";
    for (const auto &deutron: deutrons) {
        output << std::fixed << std::setprecision(15)
               << deutron.freeze_out_time << " " << deutron.x << " "
               << deutron.y << " " << deutron.z << " "
               << deutron.px << " " << deutron.py << " "
               << deutron.pz << " " << deutron.p0 << " "
               << deutron.mass << " " << deutron.probability << "\n";
    }
}


void output_d_mix_spv(const std::vector<double> &d_mix_spv,
                      const std::vector<double> &d_mix_ptv,
                      const std::string &filename,
                      int total_batch) {
    std::ofstream output_file(filename, std::ios::out);
    if (output_file.is_open()) {
        for (int i = 0; i < d_mix_spv.size(); i++) {
            double spv = d_mix_spv[i] / total_batch;
            output_file << std::setw(15) << d_mix_ptv[i] << std::setw(15) << spv << std::endl;
        }
    }
    output_file.close();
}