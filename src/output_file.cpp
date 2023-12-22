//
// Created by mafu on 12/15/2023.
//
#include "../include/output_file.h"

void writeNucleiData(ParticleData &nuclei, std::ofstream &outputFile) {
    nuclei.getFreezeOutPosition();
    outputFile << nuclei.t << " " << nuclei.x << " " << nuclei.y << " " << nuclei.z << " "
               << nuclei.px << " " << nuclei.py << " " << nuclei.pz << " " << nuclei.p0 << " "
               << nuclei.mass << " " << nuclei.freeze_out_time << "\n";
}
void outputCluster(const std::vector<ParticleData> &clusters, std::ofstream &output) {
    output << "t x y z px py pz p0 mass probability\n";
    for (const auto &cluster: clusters) {
        output << cluster.freeze_out_time << " " << cluster.x << " " << cluster.y << " "
               << cluster.z << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
               << cluster.p0 << " " << cluster.mass << " " << cluster.probability << "\n";
    }
}


void outputPt(std::map<std::string, std::vector<double>> &pt_array,
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
