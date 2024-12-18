#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <matplot/matplot.h>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void plot_results(Eigen::VectorXd &omega_squared,
                  Eigen::MatrixXd &displacements, int n_masses, int index) {
  // Exibir os resultados
  std::cout << "Frequências de vibração (ω²):\n" << omega_squared << "\n\n";

  // Exibir gráficos
  // Densidade de estados para ω²
  matplot::figure();
  matplot::hist(omega_squared);
  matplot::title("Densidade de estados para ω²");
  matplot::xlabel("ω²");
  matplot::ylabel("Densidade de estados");
  // matplot::show();
  auto filename = "densidade_de_estados_" + std::to_string(index) + "_" +
                  std::to_string(n_masses);
  matplot::save("output/" + filename, "svg");

  // Deslocamento relativo para uma frequência

  auto h = matplot::figure();
  h->size(1920, 1080);
  for (auto i = 0; i < 5; i++) {
    
    matplot::subplot(2, 5, i);
    matplot::plot(VectorXd(displacements.col(i)));
  }

  for (auto i = n_masses - 5; i < n_masses; i++) {
    matplot::subplot(2, 5, i - n_masses + 10);
    matplot::plot(VectorXd(displacements.col(i)));
  }
  // matplot::show();
  filename = "deslocamento_relativo_" + std::to_string(index) + "_" +
             std::to_string(n_masses);
  matplot::save("output/" + filename, "svg");
}

MatrixXd dynamic_matrix(int n, std::vector<double> &k, std::vector<double> &m) {

  Eigen::MatrixXd D(n, n);

  for (int i = 0; i < n; i++) {
    // D[i][i] = (k[i-1][i] + k[i][i+1])/m[i] for i == j
    // D[i][i+1] = -k[i][i+1]/m[i]
    // D[i][i-1] = -k[i-1][i]/m[i]
    if (i == 0) {
      D(i, i) = (k[i]) / m[i];
      D(i, i + 1) = -k[i] / m[i];
    } else if (i == n - 1) {
      D(i, i) = (k[i - 1]) / m[i];
      D(i, i - 1) = -k[i - 1] / m[i];
    } else {
      D(i, i) = (k[i - 1] + k[i]) / m[i];
      D(i, i + 1) = -k[i] / m[i];
      D(i, i - 1) = -k[i - 1] / m[i];
    }
  }

  return D;
}

int main() {
  // Calcule as frequências de vibração e os conjuntos de deslocamentos para uma
  // cadeia unidimensional de 100, 500 e depois de 1000 corpos de massas iguais
  // (m = 1 kg), todos unidos por molas iguais (k = 1 kg/s2), exceto para os k
  // das duas molas que ligam a massa central com a cadeia, que deve ter
  // k=5kg/s2. Apresente os resultados em gráficos da densidade de estados para
  // ω2 e os deslocamentos relativos para as cinco primeiras e cinco últimas
  // frequências (em ordem crescente). Investigue outras frequência que achar
  // interessante. Compare seus resultados para a cadeia unidimensional de 100,
  // 500 e depois de 1000 corpos de massas iguais.

  auto ns = {100, 500, 1000}; // number of masses

  for (auto n : ns) {
    std::vector<double> k(n - 1, 1.0);
    std::vector<double> m(n, 1.0);

    MatrixXd D_1 = dynamic_matrix(n, k, m);

    // Resolver autovalores e autovetores
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver_1(D_1);

    VectorXd omega_squared = solver_1.eigenvalues(); // ω²
    MatrixXd displacements = solver_1.eigenvectors();

    plot_results(omega_squared, displacements, n, 1);

    k[(n / 2) - 1] = 5.0;
    k[(n / 2)] = 5.0;
    k[(n / 2) - 2] = 5.0;

    MatrixXd D_2 = dynamic_matrix(n, k, m);

    // Resolver autovalores e autovetores
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver_2(D_2);

    omega_squared = solver_2.eigenvalues(); // ω²
    displacements = solver_2.eigenvectors();

    plot_results(omega_squared, displacements, n, 2);
  }
  return 0;
}
