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
        output << cluster.freeze_out_time << " " << cluster.x << " " << cluster.y << " "
               << cluster.z << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
               << cluster.p0 << " " << cluster.mass << " " << cluster.probability << "\n";
    }
}


void output_spv(std::map<std::string, std::vector<double>> &pt_array,
                std::map<std::string, double> clusterCountByRapidity, double d_pt, int ptBins,
                const std::string &filename, int total_batch) {
    std::ofstream output_file(filename, std::ios::out);
    if (!output_file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }
    for (auto &[label, pts]: pt_array) {
        output_file << "Rapidity range: " << label
                    << ", cluster yield: " << clusterCountByRapidity[label] / total_batch
                    << std::endl;
        for (size_t i = 0; i < ptBins; ++i) {
            double pt = d_pt / 2 + static_cast<double>(i) * d_pt;
            pts[i] /= total_batch;
            output_file << pt << "\t" << pts[i] << std::endl;
        }
        output_file << std::endl;
    }
    output_file.close();
}