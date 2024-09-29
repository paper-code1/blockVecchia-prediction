#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H

#include <vector>
#include <array>
#include "Point.h"

// Define T as double
using T = double;

class ClusterData {
public:
    std::vector<std::vector<T>> clustersLocations;
    std::vector<T> observations;
    int numPoints;
    std::vector<T> centroids;
    std::vector<std::vector<T>> nearestNeighborsLocations;
    std::vector<T> nearestNeighborsObservations;
    int dimension;  // 2 for 2D, 3 for 3D
    std::vector<std::vector<T>> covarianceMat_cluster;  // 2D array for covariance matrix for cluster
    std::vector<std::vector<T>> covarianceMat_nearestNeighbors; 
    std::vector<std::vector<T>> covarianceMat_cross;
    // predicted values, observations and uncertainties, mspe
    std::vector<T> predictedValues;
    std::vector<std::vector<T>> predictedUncertainties;
    T mspe;
    int clusterIndex;
    // calcuate the mean and variance of the conditional simulation
    std::vector<double> mean;
    std::vector<double> variance;


    // constructor and destructor
    ClusterData(int dim): numPoints(0), dimension(dim), mspe(0), mean(0), variance(0) {};
    ~ClusterData() {};

    void generateCovarianceMatrix(const std::vector<T>& theta, int distance_metric, const double scale_factor=1.0);
    void krigingPredict();
    void conditionalSimulate(int m_replicates);
    // Add a new method to save cluster and neighbor locations to CSV
    void saveToCSV(const std::string& filename) const;

};

// Add a new function declaration
void saveAllClustersToCSV(const std::vector<ClusterData>& clusters, const std::string& filename);


// Function declarations
double euclideanDistance(const std::vector<T>& p1, const std::vector<T>& p2);

std::vector<int> findNearestNeighborsIndices(const std::vector<T>& centroid, 
                                             const std::vector<std::vector<T>>& trainLocs, 
                                             int m);

std::vector<ClusterData> constructClusterData(
    const std::vector<Point>& points, 
    const std::vector<T>& observations, 
    const std::vector<Point>& centroids_points, 
    const std::vector<std::vector<T>>& trainLocs, 
    const std::vector<T>& trainObservations, 
    int k, int m, int num_threads, int dimension);

#endif // CLUSTER_UTILS_H