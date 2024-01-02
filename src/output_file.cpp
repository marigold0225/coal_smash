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
//void outputCluster(const std::vector<ParticleData> &clusters, int events, std::ofstream &output) {
//    output << "Number of events: " << events << "\n";
//    output << "t x y z px py pz p0 mass probability\n";
//    for (const auto &cluster: clusters) {
//        output << cluster.freeze_out_time << " " << cluster.x << " " << cluster.y << " "
//               << cluster.z << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
//               << cluster.p0 << " " << cluster.mass << " " << cluster.probability << "\n";
//    }
//}
void outputCluster(const std::vector<ParticleData> &clusters, int events, std::ofstream &output) {
    output << events << " " << clusters.size() << " " << 0 << " " << 0 << " " << 0 << " " << 0
           << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";
    for (const auto &cluster: clusters) {
        output << 4 << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
               << cluster.mass << " " << cluster.x << " " << cluster.y << " " << cluster.z << " "
               << cluster.freeze_out_time << " " << cluster.probability << "\n";
    }
}

void outputPt(ptArray &pt_array, std::map<std::string, double> clusterCountByRapidity,
              const reactionConfig &clusterConfig, const std::string &filename, int total_batch) {
    std::ofstream output_file(filename, std::ios::out);
    if (!output_file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }
    for (auto &[label, pts]: pt_array) {
        output_file << "Rapidity range: " << label
                    << ", cluster yield: " << clusterCountByRapidity[label] / total_batch
                    << std::endl;
        for (int i = 0; i < clusterConfig.ptBins; ++i) {
            double pt = clusterConfig.dpt / 2 + i * clusterConfig.dpt;
            pts[i] /= total_batch;
            output_file << pt << "\t" << pts[i] << std::endl;
        }
        output_file << std::endl;
    }
    output_file.close();
}
void outputClusterOther(const std::pair<std::deque<ParticleData>, int> &cluster_event,
                        std::ofstream &output) {
    //    output << "Number of events: " << cluster_event.second << "\n";
    //    output << "t x y z px py pz p0 mass probability\n";
    output << cluster_event.second << " " << cluster_event.first.size() << " " << 0 << " " << 0
           << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << "\n";
    for (const auto &cluster: cluster_event.first) {
        //        output << cluster.freeze_out_time << " " << cluster.x << " " << cluster.y << " "
        //               << cluster.z << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
        //               << cluster.p0 << " " << cluster.mass << " " << cluster.probability << "\n";
        output << 4 << " " << cluster.px << " " << cluster.py << " " << cluster.pz << " "
               << cluster.mass << " " << cluster.x << " " << cluster.y << " " << cluster.z << " "
               << cluster.freeze_out_time << " " << cluster.probability << "\n";
    }
}
