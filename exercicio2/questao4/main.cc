#include <cstdint>
#include <fstream>
#include <iostream>
#include <matplot/matplot.h>
#include <sstream>
#include <vector>

int main() {
  std::ifstream file("exercicio2/questao4/sunspots.txt");
  if (!file.is_open()) {
    std::cerr << "Failed to open the file." << std::endl;
    return 1;
  }

  std::vector<int> indices;
  std::vector<double> values;
  std::string line;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    int index;
    double value;
    if (iss >> index >> value) {
      indices.push_back(index);
      values.push_back(value);
    }
  }

  file.close();

  int r = 5;
  int N = 1000;

  std::vector<double> running_avg;

  for (auto i = 0; i < values.size(); i++) {
    auto sum = 0;

    for (auto j = -r; j <= r; j++) {
      if (i + j >= 0 && i + j < values.size()) {
        sum += values[i + j];
      }
    }
    running_avg.push_back(sum / (2 * r + 1));
  }

  auto fig1 = matplot::figure(true);
  matplot::plot(indices, values);
  matplot::hold(matplot::on);
  matplot::plot(indices, running_avg);
  matplot::show();

  auto fig2 = matplot::figure(true);
  matplot::plot(std::vector<int>(indices.begin(), indices.begin() + N),
                std::vector<double>(values.begin(), values.begin() + N));
  matplot::hold(matplot::on);
  matplot::plot(
      std::vector<int>(indices.begin(), indices.begin() + N),
      std::vector<double>(running_avg.begin(), running_avg.begin() + N));
  matplot::show();

  return 0;
}