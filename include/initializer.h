//
// Created by dragos on 31/10/22.
//
#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <random>
#include <cmath>
#include "Point.h"

// Function definition for loadOneDimensionalData (to be placed in an include file)

inline std::vector<double> loadOneDimensionalData(const std::string &filename)
{
    std::vector<double> data;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        double value;
        if (ss >> value)
        {
            data.push_back(value);
        }
    }
    return data;
}

inline std::vector<Point> convertToPoints(std::vector<std::vector<double>> &loc, int n)
{
    std::vector<Point> points;
    points.reserve(n); // Optimize memory allocation

    for (int i = 0; i < n; ++i)
    {
        Point p; // Use the default constructor
        for (int j = 0; j < DIM; j++)
        {
            p.coordinates[j] = loc[i][j];
        }
        p.cluster = -1;      // Assuming you want to initialize the cluster to -1 or some default value
        points.push_back(p); // Add the point to the vector
    }

    return points;
}

/*
 * This function checks if a point is already in a vector of points
 * @param points: the vector of points
 * @param p: the point to be checked
 * @return true if the point is in the vector, false otherwise
 */
inline bool contatins(std::vector<Point> &points, Point &p)
{
    for (Point el : points)
    {
        if (el == p)
        {
            return true;
        }
    }
    return false;
}

/*
 * This function randomly initializes the centroids
 * @param points: the vector of points
 * @param k: the number of clusters
 * @return a vector of centroids
 */
inline std::vector<Point> random_initializer(const std::vector<Point> &points, const int k, int seed)
{
    std::vector<Point> centroids;
    Point p;
    centroids.reserve(k);
    // select points with equal interval
    for (int i = 0; i < int(points.size()); i += int(points.size() / k))
    {
        p = points[i];
        p.cluster = i;
        centroids.push_back(p);
    }
    // std::random_device rand_dev;
    // std::mt19937 gen(rand_dev());
    // std::mt19937 gen(seed);
    // std::uniform_int_distribution<> distrib(0, points.size());
    // int i = 0;

    // for (auto & point: points) point.print();

    // while (int(centroids.size()) < k)
    // {
    //     Point new_centroid = points[distrib(gen)];
    //     if (!contatins(centroids, new_centroid))
    //     {
    //         new_centroid.cluster = i;
    //         i++;
    //         // new_centroid.print();
    //         centroids.push_back(new_centroid);
    //     }
    // }
    return centroids;
}

/*
 * This function chooses the next centroid based on the weighted probability
 * @param points: the vector of points
 * @param distances: the vector of distances
 * @return the next centroid
 */
inline Point next_centroid(const std::vector<Point> &points, const std::vector<double> &distances)
{
    std::random_device rand_dev;
    std::mt19937 gen(rand_dev());
    std::discrete_distribution<> distrib(distances.begin(), distances.end());
    int index = distrib(gen);
    return points[index];
}

/*
 * This function initializes centroids by choosing the farthest points from each other
 * @param points: the vector of points
 * @param k: the number of clusters
 * @return a vector of centroids
 */
inline std::vector<Point> kmeanpp_initializer(const std::vector<Point> &points, const int k, const int threads)
{
    std::vector<Point> centroids;
    centroids.reserve(k);
    int block = ceil(points.size() / threads);

    std::random_device rand_dev;
    std::mt19937 gen(rand_dev());
    std::uniform_int_distribution<> distrib(0, points.size() - 1);

    // the first centroid is chosen randomly
    centroids.push_back(points[distrib(gen)]);

    // the other centroids are chosen based on the weighted probability
    while (int(centroids.size()) < k)
    {
        std::vector<std::vector<double>> privateDistances(threads, std::vector<double>(points.size(), std::numeric_limits<double>::max()));

#pragma omp parallel for num_threads(threads) schedule(static, block)
        for (size_t i = 0; i < points.size(); i++)
        {
            for (const Point &centroid : centroids)
            {
                double distance = euclidean_dist(points[i], centroid);
                privateDistances[omp_get_thread_num()][i] = std::min(privateDistances[omp_get_thread_num()][i], distance);
            }
        }

        // Merge private distances from all threads
        std::vector<double> distances(points.size(), std::numeric_limits<double>::max());
        for (int t = 0; t < threads; t++)
        {
            for (size_t i = 0; i < points.size(); i++)
            {
                distances[i] = std::min(distances[i], privateDistances[t][i]);
            }
        }
        Point nextCentroid = next_centroid(points, distances);
        centroids.push_back(nextCentroid);
    }
    return centroids;
}
#endif