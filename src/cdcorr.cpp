#include <cmath>
#include <vector>
#include <tuple>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>

using namespace std;

void impMat(std::vector<std::vector<double>> &A)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void printVec(std::vector<double> &v)
{
    for (auto i : v)
        std::cout << i << ' ';
}

std::vector<std::vector<double>> makeOnes(double n)
{
    std::vector<std::vector<double>> mat;
    std::vector<double> temp;

    for (double i = 0; i < n; i++)
    {
        for (double j = 0; j < n; j++)
        {
            temp.push_back(i);
        }
        mat.push_back(temp);
        temp.clear();
    }
    return mat;
}

std::vector<std::vector<double>> makeCovariate(double nrow, double ncol)
{
    std::vector<std::vector<double>> mat;
    std::vector<double> temp;

    double counter = 0;
    for (double i = 0; i < nrow; i++)
    {
        for (double j = 0; j < ncol; j++)
        {
            temp.push_back(counter);
            counter++;
        }
        mat.push_back(temp);
        temp.clear();
    }
    return mat;
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

double vector_prod(std::vector<double> &vector1)
{
    double sum_value = 1.0;
    for (double value : vector1)
    {
        sum_value *= value;
    }
    return sum_value;
}

double vector_mean(std::vector<double> &vector1)
{
    double vector_sum_value = vector_sum(vector1);
    return vector_sum_value / vector1.size();
}

/**
 * V-type conditional distance covariance
 * @param distance_x : pairwise distance of variable x
 * @param distance_y : pairwise distance of variable y
 * @param kernel_density_estimation : KDE result of conditional variable
 * @return the value of conditional distance covariance test statistic
 * @refitem Xueqin Wang, Wenliang Pan, Wenhao Hu, Yuan Tian & Heping Zhang (2015) Conditional Distance Correlation,
 * Journal of the American Statistical Association, 110:512, 1726-1734, DOI: 10.1080/01621459.2014.993081
 */
std::vector<double> compute_condition_distance_covariance_crude(std::vector<std::vector<double>> &distance_x,
                                                                std::vector<std::vector<double>> &distance_y,
                                                                std::vector<std::vector<double>> &kernel_density_estimation)
{
    size_t num = distance_x.size();
    std::vector<double> condition_distance_covariance(num, 0.0);
    double cross_dx1, cross_dy1, cross_dx2, cross_dy2, cross_dx3, cross_dy3, d_ijkl;
    for (size_t u = 0; u < num; ++u)
    {
        for (size_t i = 0; i < num; ++i)
        {
            for (size_t j = 0; j < num; ++j)
            {
                for (size_t k = 0; k < num; ++k)
                {
                    for (size_t l = 0; l < num; ++l)
                    {
                        cross_dx1 = distance_x[i][j] + distance_x[k][l] - distance_x[i][k] - distance_x[j][l];
                        cross_dy1 = distance_y[i][j] + distance_y[k][l] - distance_y[i][k] - distance_y[j][l];
                        cross_dx2 = distance_x[i][j] + distance_x[k][l] - distance_x[i][l] - distance_x[j][k];
                        cross_dy2 = distance_y[i][j] + distance_y[k][l] - distance_y[i][l] - distance_y[j][k];
                        cross_dx3 = distance_x[i][l] + distance_x[k][j] - distance_x[i][k] - distance_x[j][l];
                        cross_dy3 = distance_y[i][l] + distance_y[k][j] - distance_y[i][k] - distance_y[j][l];
                        d_ijkl = cross_dx1 * cross_dy1 + cross_dx2 * cross_dy2 + cross_dx3 * cross_dy3;

                        condition_distance_covariance[u] += d_ijkl *
                                                            kernel_density_estimation[i][u] *
                                                            kernel_density_estimation[j][u] *
                                                            kernel_density_estimation[k][u] *
                                                            kernel_density_estimation[l][u];
                    }
                }
            }
        }
        condition_distance_covariance[u] /= (num * num * num * num);
    }
    return condition_distance_covariance;
}

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

int main()
{
    // vector<vector<double>> A;
    // vector<vector<double>> B;
    // vector<double> temp;

    // double n = 5;

    // for (double j = 0; j < n; j++)
    // {
    //     for (double i = 0; i < n; i++)
    //     {
    //         temp.push_back(i);
    //     }
    //     A.push_back(temp);
    //     B.push_back(temp);
    //     temp.clear();
    // }

    vector<vector<double>> mat = makeCovariate(5, 3);
    impMat(mat);

    vector<vector<double>> res = compute_gaussian_kernel_estimate(mat, 0.5);

    impMat(res);

    // printVec(cdcov);
    return 0;
}
