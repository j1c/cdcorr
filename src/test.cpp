#include "csv.h"
#include <iostream>

using namespace std;

void printMat(std::vector<std::vector<double>> &A)
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
    io::CSVReader<3> in("./x.csv");
    in.read_header(io::ignore_extra_column, "a", "b", "c");
    std::vector<std::vector<double>> data(20, std::vector<double>(3));

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

    printMat(data);

    return 0;
}
