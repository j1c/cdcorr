#include <cmath>
#include <vector>
#include <tuple>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include "csv.h"
#include "utils.h"

using namespace std;

std::vector<double> compute_condition_distance_covariance(
    std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
    std::vector<std::vector<double>> &kernel_density_estimation)
{
    // initializing some variables
    uint num = (uint)distance_x.size();
    std::vector<std::vector<double>> anova_x(num, std::vector<double>(num)); // num x num matrix of zeros
    std::vector<std::vector<double>> anova_y(num, std::vector<double>(num)); // num x num matrix of zeros
    std::vector<double> condition_distance_covariance(num);                  // vector of zeros of size num

    double kernel_sum_square;
    for (uint i = 0; i < num; i++)
    {
        anova_x = weight_distance_anova(distance_x, kernel_density_estimation[i]); // fill the contents
        anova_y = weight_distance_anova(distance_y, kernel_density_estimation[i]); // fill the contents

        kernel_sum_square = vector_sum(kernel_density_estimation[i]);
        kernel_sum_square = kernel_sum_square * kernel_sum_square;

        for (uint k = 0; k < num; k++)
        {
            for (uint j = 0; j < num; j++)
            {
                condition_distance_covariance[i] += anova_x[k][j] * anova_y[k][j] *
                                                    kernel_density_estimation[i][k] *
                                                    kernel_density_estimation[i][j];
            }
        }
        condition_distance_covariance[i] /= kernel_sum_square;
    }

    return condition_distance_covariance;
}

std::vector<double> compute_condition_distance_correlation(
    std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
    std::vector<std::vector<double>> &kernel_density_estimation)
{

    uint num = (uint)distance_x.size();
    std::vector<std::vector<double>> anova_x(num, std::vector<double>(num));
    std::vector<std::vector<double>> anova_y(num, std::vector<double>(num));
    std::vector<double> condition_distance_covariance_xy(num);
    std::vector<double> condition_distance_covariance_xx(num);
    std::vector<double> condition_distance_covariance_yy(num);

    for (uint i = 0; i < num; i++)
    {
        anova_x = weight_distance_anova(distance_x, kernel_density_estimation[i]);
        anova_y = weight_distance_anova(distance_y, kernel_density_estimation[i]);

        for (uint k = 0; k < num; k++)
        {
            for (uint j = 0; j < num; j++)
            {
                condition_distance_covariance_xy[i] += anova_x[k][j] * anova_y[k][j] *
                                                       kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
                condition_distance_covariance_xx[i] += anova_x[k][j] * anova_x[k][j] *
                                                       kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
                condition_distance_covariance_yy[i] += anova_y[k][j] * anova_y[k][j] *
                                                       kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
            }
        }
    }

    double dcor_denominator;
    for (uint i = 0; i < num; i++)
    {
        dcor_denominator = condition_distance_covariance_xx[i] * condition_distance_covariance_yy[i];
        if (dcor_denominator > 0.0)
        {
            condition_distance_covariance_xy[i] /= sqrt(dcor_denominator);
        }
        else
        {
            condition_distance_covariance_xy[i] = 0.0;
        }
    }

    return condition_distance_covariance_xy;
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
double compute_condition_distance_covariance_crude(std::vector<std::vector<double>> &distance_x,
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
    double condition_distance_covariance_stat = vector_mean(condition_distance_covariance);
    return condition_distance_covariance_stat;
    // return condition_distance_covariance;
}

std::vector<std::vector<double>> read_csv(std::string path, int nrow, int ncol)
{
    io::CSVReader<3> in(path);
    in.read_header(io::ignore_extra_column, "a", "b", "c");
    std::vector<std::vector<double>> data(nrow, std::vector<double>(ncol));

    double a;
    double b;
    double c;
    int i = 0;
    while (in.read_row(a, b, c))
    {
        data[i][0] = a;
        data[i][1] = b;
        data[i][2] = c;
        i++;
    }

    return data;
}

int main()
{
    // std::vector<std::vector<double>> x = read_csv("../data/x_10.csv", 10, 3);
    // std::vector<std::vector<double>> y = read_csv("../data/y_10.csv", 10, 3);
    // std::vector<std::vector<double>> z = read_csv("../data/z_10.csv", 10, 3);

    std::vector<std::vector<double>> x = read_csv("../data/x.csv", 100, 3);
    std::vector<std::vector<double>> y = read_csv("../data/y.csv", 100, 3);
    std::vector<std::vector<double>> z = read_csv("../data/z.csv", 100, 3);

    std::vector<std::vector<double>> dx = Euclidean_distance(x, 1);
    std::vector<std::vector<double>> dy = Euclidean_distance(y, 1);
    std::vector<std::vector<double>> dz = compute_gaussian_kernel_estimate(z, 0.5);

    // uint num = (uint)dx.size();
    // std::vector<std::vector<double>> anova_x(num, std::vector<double>(num));
    // std::vector<std::vector<double>> anova_y(num, std::vector<double>(num));
    // std::vector<double> condition_distance_covariance(num);

    // double kernel_sum_square;

    // for (uint i = 0; i < num; i++)
    // {
    //     anova_x = weight_distance_anova(dx, dz[i]);
    //     anova_y = weight_distance_anova(dy, dz[i]);

    //     kernel_sum_square = vector_sum(dz[i]);
    //     kernel_sum_square = kernel_sum_square * kernel_sum_square;

    //     for (uint k = 0; k < num; k++)
    //     {
    //         for (uint j = 0; j < num; j++)
    //         {
    //             condition_distance_covariance[i] += anova_x[k][j] * anova_y[k][j] *
    //                                                 dz[i][k] *
    //                                                 dz[i][j];
    //         }
    //     }
    //     condition_distance_covariance[i] /= kernel_sum_square;
    //     condition_distance_covariance[i] *= 12 * pow(vector_sum(dz[i]) / num, 4.0);
    // }
    // printVec(condition_distance_covariance);

    // std::cout << "\n"
    //           << vector_mean(condition_distance_covariance);

    std::mt19937_64 random_number_generator;
    std::vector<std::vector<uint>> random_sample_index = generate_random_sample_index(5, dz,
                                                                                      random_number_generator);

    printMatInt(random_sample_index);

    return 0;
}
