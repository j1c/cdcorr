#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>
#include <vector>
#include <tuple>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>

double square_Euclidean_distance(std::vector<double> &vector1, std::vector<double> &vector2)
{
    double distance_value = 0.0;
    for (uint i = 0; i < vector1.size(); ++i)
    {
        distance_value += pow(vector1[i] - vector2[i], 2.0);
    }
    return distance_value;
}

std::vector<std::vector<double>> compute_gaussian_kernel_estimate(
    std::vector<std::vector<double>> &condition_variable, double bandwidth)
{
    double CDC_PI = 3.14159265358979323846;
    uint num = (uint)condition_variable.size();
    uint d = (uint)condition_variable[0].size();
    double det = pow(bandwidth, (double)d);
    double density = 1.0 / (pow(2 * CDC_PI, d / 2.0) * sqrt(det));

    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    for (uint i = 0; i < num; i++)
    {
        kernel_density_estimate[i][i] = density;
        for (uint j = 0; j < i; j++)
        {
            kernel_density_estimate[j][i] = exp(
                -0.5 * square_Euclidean_distance(condition_variable[i], condition_variable[j]) / bandwidth);
            kernel_density_estimate[j][i] *= density;
            kernel_density_estimate[i][j] = kernel_density_estimate[j][i];
        }
    }

    return kernel_density_estimate;
}

/**
 * Compute the Euclidean distance matrix
 * @param matrix: interpret x as an matrix with size n*d
 * @param index: see energy distance for reference
 */
std::vector<std::vector<double>> Euclidean_distance(std::vector<std::vector<double>> &matrix, double index)
{
    uint n = (uint)matrix.size();
    uint d = (uint)matrix[0].size();
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));

    double diff_sum, diff;
    for (uint i = 1; i < n; i++)
    {
        distance_matrix[i][i] = 0.0;
        for (uint j = 0; j < i; j++)
        {
            diff_sum = 0.0;
            for (uint k = 0; k < d; k++)
            {
                diff = matrix[i][k] - matrix[j][k];
                diff_sum += diff * diff;
            }
            distance_matrix[i][j] = distance_matrix[j][i] = pow(sqrt(diff_sum), index);
        }
    }
    return distance_matrix;
}

void printMat(std::vector<std::vector<double>> &A)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void printMatInt(std::vector<std::vector<uint>> &A)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVec(std::vector<double> &A)
{
    for (int i = 0; i < A.size(); i++)
    {

        std::cout << A[i] << " ";
    }
    std::cout << std::endl;
}

/**
 * Compute the summation of 1D std::vector
 * @param vector1 : 1D vector
 * @return summation value
 */
double vector_sum(std::vector<double> &vector1)
{
    double sum_value = 0.0;
    for (double value : vector1)
    {
        sum_value += value;
    }
    return sum_value;
}

double vector_mean(std::vector<double> &vector1)
{
    double vector_sum_value = vector_sum(vector1);
    return vector_sum_value / vector1.size();
}

double vector_weight_sum(std::vector<double> &vector1, std::vector<double> &weight)
{
    double sum_value = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i)
    {
        sum_value += vector1[i] * weight[i];
    }
    return sum_value;
}

std::vector<std::vector<double>> weight_distance_anova(std::vector<std::vector<double>> &distance_matrix,
                                                       std::vector<double> &weight)
{
    double weight_sum = vector_sum(weight);
    uint num = (uint)distance_matrix.size();

    std::vector<double> marginal_weight_distance(num);
    for (uint i = 0; i < num; ++i)
    {
        marginal_weight_distance[i] = vector_weight_sum(distance_matrix[i], weight);
    }

    double weight_distance_sum = 0.0;
    for (uint i = 0; i < num; ++i)
    {
        weight_distance_sum = vector_weight_sum(marginal_weight_distance, weight);
    }
    weight_distance_sum /= weight_sum * weight_sum;

    for (uint i = 0; i < num; ++i)
    {
        marginal_weight_distance[i] /= weight_sum;
    }

    std::vector<std::vector<double>> weight_distance_anova_table(num, std::vector<double>(num));
    for (uint k = 0; k < num; k++)
    {
        for (uint j = k; j < num; j++)
        {
            weight_distance_anova_table[k][j] =
                distance_matrix[k][j] - marginal_weight_distance[k] - marginal_weight_distance[j] +
                weight_distance_sum;
            weight_distance_anova_table[j][k] = weight_distance_anova_table[k][j];
        }
    }

    return weight_distance_anova_table;
}

uint sample_multinomial_distribution(std::vector<double> &probability, std::mt19937_64 &random_number_generator)
{
    std::discrete_distribution<uint> multinomial_sampler(probability.begin(), probability.end());
    return multinomial_sampler(random_number_generator);
}

std::vector<std::vector<uint>> generate_random_sample_index(uint replication_number,
                                                            std::vector<std::vector<double>> probability_matrix,
                                                            std::mt19937_64 &random_number_generator)
{
    std::vector<std::vector<uint>> random_sample_index(replication_number,
                                                       std::vector<uint>(probability_matrix.size()));
    for (uint i = 0; i < probability_matrix.size(); ++i)
    {
        for (uint j = 0; j < replication_number; ++j)
        {
            random_sample_index[j][i] = sample_multinomial_distribution(probability_matrix[i], random_number_generator);
        }
    }

    return random_sample_index;
}

#endif // SRC_UTILITY_H
