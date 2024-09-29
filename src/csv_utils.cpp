
#include "csv_utils.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <filesystem>
#include <vector>


std::vector<std::vector<double>> loadCSV(const std::string &filename, int dim)
{
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::vector<double> row(dim);
        for (int i = 0; i < dim; ++i)
        {
            ss >> row[i];
            if (ss.peek() == ',')
                ss.ignore();
        }
        data.push_back(row);
    }
    return data;
}

// In the corresponding .cpp file, implement the function:
void saveClusterInfo(const std::vector<ClusterData>& clusters, const std::vector<Point>& centroids, const std::string& filename) {
    std::ofstream file(filename);
    file << "type,cluster_id,x,y,z\n";

    for (size_t i = 0; i < clusters.size(); ++i) {
        // Save centroid
        file << "centroid," << clusters[i].clusterIndex;
        for (double coord : centroids[i].coordinates) {
            file << "," << coord;
        }
        file << "\n";

        // Save cluster points
        for (const auto& point : clusters[i].clustersLocations) {
            file << "point," << clusters[i].clusterIndex;
            for (double coord : point) {
                file << "," << coord;
            }
            file << "\n";
        }

        // Save nearest neighbors
        for (const auto& neighbor : clusters[i].nearestNeighborsLocations) {
            file << "neighbor," << clusters[i].clusterIndex;
            for (double coord : neighbor) {
                file << "," << coord;
            }
            file << "\n";
        }
    }
    file.close();
}

void writeResultsToCSV(const std::vector<ClusterData>& clusters, const std::vector<double>& theta, double mspe, int k, int m, int seed) {
    // save the mean and variance of the conditional simulation as a csv file
    // create log folder
    std::string log_folder = "log";
    if (!std::filesystem::exists(log_folder)) {
        std::filesystem::create_directory(log_folder);
    }
    std::ofstream csvFile("log/conditional_simulation_k_" + std::to_string(k) + "_m_" + std::to_string(m) + "_theta_" + std::to_string(theta[0]) + "_" + std::to_string(theta[1]) + "_" + std::to_string(theta[2]) + "_seed_" + std::to_string(seed) + ".csv");
    // header
    csvFile << "smean,svariance,x,y,z,mspe,true_value\n";
    for (const auto& cluster : clusters) {
        for (int i = 0; i < cluster.numPoints; ++i) {
            csvFile << cluster.mean[i] << "," << cluster.variance[i] << "," << cluster.clustersLocations[i][0] << "," << cluster.clustersLocations[i][1] << "," << cluster.clustersLocations[i][2] << "," << mspe << "," << cluster.observations[i]<< std::endl;
        }
    }
    csvFile.close();
}