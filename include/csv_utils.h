
#ifndef CSV_UTILS_H
#define CSV_UTILS_H

#include <vector>
#include <string>
#include "Point.h"
#include "cluster_utils.h"
std::vector<std::vector<double>> loadCSV(const std::string& filename, int dim);

// Add this function declaration
void saveClusterInfo(const std::vector<ClusterData>& clusters, const std::vector<Point>& centroids, const std::string& filename);

void writeResultsToCSV(const std::vector<ClusterData>& clusters, const std::vector<double>& theta, double mspe, int k, int m, int seed);
#endif
