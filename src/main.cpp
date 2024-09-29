
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include "csv_utils.h"
#include "kmeans.h"
#include "cluster_utils.h"

// Argument parsing library (you can use argparse or getopt)
#include "cxxopts.hpp"


int main(int argc, char *argv[])
{
    // Parse command-line arguments
    cxxopts::Options options("Prediction using the block Vecchia approximation", "Perform prediction using the block Vecchia approximation");

    options.add_options()
    ("train_locs", "Locations path to the training CSV file", cxxopts::value<std::string>())
    ("test_locs", "Locations path to the testing CSV file", cxxopts::value<std::string>())
    ("train_data", "Data path to the training CSV file", cxxopts::value<std::string>())
    ("test_data", "Data path to the testing CSV file", cxxopts::value<std::string>())
    ("seed", "Seed for the random number generator", cxxopts::value<int>()->default_value("42"))
    ("n", "Number of locations", cxxopts::value<int>())
    ("k", "Number of clusters", cxxopts::value<int>()->default_value("100"))
    ("m", "Number of nearest neighbors", cxxopts::value<int>()->default_value("60"))
    ("num_threads", "Number of threads", cxxopts::value<int>()->default_value("20"))
    ("dim", "Dimension of the location data (2D or 3D)", cxxopts::value<int>()->default_value("3"))
    ("kmeans_iter", "Number of iterations for kmeans", cxxopts::value<int>()->default_value("50"))
    ("theta", "Theta for the covariance matrix, e.g., 1.0,0.5,0.1, no space between the numbers", cxxopts::value<std::vector<double>>()->default_value("1.0, 0.5, 0.1"))
    ("distance_metric", "Distance metric (1 for earth, 2 for euclidean)", cxxopts::value<int>()->default_value("2"))
    ("scale_factor", "Scale factor for the covariance matrix, For real 2D dataset (e.g., wind data, 2523.64 or soil data, 9348.317", cxxopts::value<double>()->default_value("1.0"))
    ("conditional_sim", "number of conditional simulations, e.g., 1000", cxxopts::value<int>()->default_value("1000"))
    ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }
    // Load parameters
    std::string trainLocsFile = result["train_locs"].as<std::string>();
    std::string testLocsFile = result["test_locs"].as<std::string>();
    std::string trainDataFile = result["train_data"].as<std::string>();
    std::string testDataFile = result["test_data"].as<std::string>();
    int seed = result["seed"].as<int>();
    int n = result["n"].as<int>();
    int k = result["k"].as<int>();
    int m = result["m"].as<int>();
    int dim = result["dim"].as<int>();
    int omp_numthreads = result["num_threads"].as<int>();
    int kmeans_iter = result["kmeans_iter"].as<int>();
    std::vector<double> theta = result["theta"].as<std::vector<double>>();
    int distance_metric = result["distance_metric"].as<int>();
    double scale_factor = result["scale_factor"].as<double>();
    int conditional_sim = result["conditional_sim"].as<int>();

    // print the parameters
    std::cout << "trainLocsFile: " << trainLocsFile << std::endl;
    std::cout << "testLocsFile: " << testLocsFile << std::endl;
    std::cout << "trainDataFile: " << trainDataFile << std::endl;
    std::cout << "testDataFile: " << testDataFile << std::endl;
    std::cout << "seed: " << seed << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "m: " << m << std::endl;

    omp_set_num_threads(omp_numthreads);

    // Load the datasets
    std::vector<double> trainData = loadOneDimensionalData(trainDataFile);
    std::vector<double> testData = loadOneDimensionalData(testDataFile);
    std::vector<std::vector<double>> trainLocs = loadCSV(trainLocsFile, dim);
    std::vector<std::vector<double>> testLocs = loadCSV(testLocsFile, dim);

    // Perform K-means clustering and find nearest neighbors
    // kmeans
    std::vector<Point> points;
    std::vector<Point> centroids;
    // transform the locations into points
    points = convertToPoints(testLocs, n);
    // init the centroids
    centroids = random_initializer(points, k, seed);
    // kmeans_iter, kmeans iterations
    kmean_par(points, centroids, kmeans_iter, k, omp_numthreads);

    // Construct data for each cluster
    std::vector<ClusterData> clusters = constructClusterData(points, testData, centroids, trainLocs, trainData, k, m, omp_numthreads, dim);
    saveAllClustersToCSV(clusters, "cluster_data.csv");

    // for each cluster, generate the covariance matrix
// #pragma omp parallel for num_threads(omp_numthreads) schedule(dynamic)
    for (int i = 0; i < int(clusters.size()); i++)
    {
        clusters[i].generateCovarianceMatrix(theta, distance_metric, scale_factor);
        clusters[i].krigingPredict();
    }

    // calculate the mspe in total
    double mspe = 0.0;
    for (int i = 0; i < int(clusters.size()); i++)
    {
        mspe += clusters[i].mspe * clusters[i].numPoints;
    }
    mspe /= n;
    if (mspe <= 0 || std::isnan(mspe)){
        std::cout << "MSPE calculation error" << std::endl;
        std::cout << "MSPE: " << mspe << std::endl;
    }else{
        std::cout << "MSPE: " << mspe << std::endl;   
    }

    // conditional simulation
#pragma omp parallel for num_threads(omp_numthreads) schedule(dynamic)
    for (int i = 0; i < int(clusters.size()); i++)
    {
        clusters[i].conditionalSimulate(conditional_sim);
    }

    writeResultsToCSV(clusters, theta, mspe, k, m, seed);

    return 0;
}
