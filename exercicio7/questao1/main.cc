#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

using namespace std;
using Matrix = vector<vector<double>>;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Função para multiplicar duas matrizes
Matrix multiplyMatrices(const Matrix &A, const Matrix &B) {
  int n = A.size();
  Matrix result(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return result;
}

// Função para transpor uma matriz
Matrix transpose(const Matrix &A) {
  int n = A.size();
  Matrix T(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      T[i][j] = A[j][i];
    }
  }
  return T;
}

// Função para realizar a decomposição QR usando o método de Gram-Schmidt
void qrDecomposition(const Matrix &A, Matrix &Q, Matrix &R) {
  int n = A.size();
  Q = Matrix(n, vector<double>(n, 0.0));
  R = Matrix(n, vector<double>(n, 0.0));

  for (int k = 0; k < n; ++k) {
    // Copiar a coluna k de A para Q
    for (int i = 0; i < n; ++i) {
      Q[i][k] = A[i][k];
    }

    // Ortogonalizar contra as colunas anteriores de Q
    for (int j = 0; j < k; ++j) {
      double dot = 0.0;
      for (int i = 0; i < n; ++i) {
        dot += Q[i][j] * A[i][k];
      }
      R[j][k] = dot;
      for (int i = 0; i < n; ++i) {
        Q[i][k] -= dot * Q[i][j];
      }
    }

    // Normalizar a coluna k de Q
    double norm = 0.0;
    for (int i = 0; i < n; ++i) {
      norm += Q[i][k] * Q[i][k];
    }
    norm = sqrt(norm);
    R[k][k] = norm;
    for (int i = 0; i < n; ++i) {
      Q[i][k] /= norm;
    }
  }
}

// Função para imprimir uma matriz
void printMatrix(const Matrix &matrix) {
  for (const auto &row : matrix) {
    for (double val : row) {
      cout << fixed << setprecision(6) << val << "\t";
    }
    cout << endl;
  }
}

void printMatrixLatexStyle(const Matrix &matrix) {
  cout << "$$\\begin{vmatrix}" << endl;
  for (const auto &row : matrix) {
    for (double val : row) {
      cout << fixed << setprecision(6) << val << " & ";
    }
    cout << "\\\\" << endl;
  }
  cout << "\\end{vmatrix}$$" << endl;
}

void printEnergyLevelsDiagramLatex(vector<double> eigenvalues) {
  vector<double> aux = eigenvalues;
  sort(aux.begin(), aux.end(), std::greater<double>());

  cout << "% ENERGY LEVELS" << endl;
  cout << "\\begin{tikzpicture}[xscale=2.8,yscale=0.8]" << endl;
  for (int n = 0; n < aux.size(); ++n) {

    cout << "\\draw[thick] (0," << n + 1 << ") --++ (1,0) node[right=1pt] {$E_{"
         << n + 1 << "}=" << aux[n] << " eV$};" << endl;
  }
  cout << "\\draw[->,thick] (-0.1,-0.2) --++ (0," << aux.size() + 1
       << ") node[left] {$E$};" << endl;
  cout << "\\end{tikzpicture}" << endl;
}

void printOrbitalPsiLatex(const Matrix &matrix) {
  cout << "$$" << endl;
  for (int i = 0; i < matrix.size(); ++i) {
    cout << "\\psi_{" << i + 1 << "}  = ";
    for (int j = 0; j < matrix[i].size(); ++j) {
      cout << fixed << setprecision(3) << "(" << matrix[j][i]
           << ") \\cdot \\phi_{" << j + 1 << "}";
      if (j < matrix[i].size() - 1) {
        cout << " + ";
      }
    }

    cout << "\\\\" << endl;
  }
  cout << "$$" << endl;
}

void sortEigenvaluesAndEigenvectors(vector<double> &eigenvalues,
                                    Matrix &eigenvectors) {
  Matrix aux = transpose(eigenvectors);
  vector<pair<double, vector<double>>> eigenvalues_eigenvectors;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    eigenvalues_eigenvectors.push_back({eigenvalues[i], aux[i]});
  }

  sort(eigenvalues_eigenvectors.begin(), eigenvalues_eigenvectors.end(),
       [](const pair<double, vector<double>> &a,
          const pair<double, vector<double>> &b) { return a.first < b.first; });

  for (int i = 0; i < eigenvalues.size(); ++i) {
    eigenvalues[i] = eigenvalues_eigenvectors[i].first;
    eigenvectors[i] = eigenvalues_eigenvectors[i].second;
  }
  eigenvectors = transpose(eigenvectors);
}

vector<double> eletronicPopulationTotal(const Matrix &eigenvectors, int n) {
  vector<double> population;
  for (int i = 0; i < eigenvectors.size(); ++i) {
    auto aux_n = n;
    double sum = 0.0;
    auto j = 0;
    while (aux_n > 0) {
      if (aux_n >= 2) {
        sum += eigenvectors[i][j] * eigenvectors[i][j] * 2;
        aux_n -= 2;
      } else {
        sum += eigenvectors[i][j] * eigenvectors[i][j];
        aux_n -= 1;
      }
      j++;
    }
    population.push_back(sum);
  }
  return population;
}

vector<double> bondOrderCyclic(const Matrix &eigenvectors, int n_electrons) {
  vector<double> bond_order;
  for (int i = 0; i < eigenvectors.size(); ++i) {
    auto aux_n = n_electrons;
    auto j = 0;
    auto sum = 0.0;
    while (aux_n > 0) {
      if (aux_n >= 2) {
        if (i == eigenvectors.size() - 1) {
          sum += eigenvectors[i][j] * eigenvectors[0][j] * 2;
        } else {
          sum += eigenvectors[i][j] * eigenvectors[i + 1][j] * 2;
        }
        aux_n -= 2;
      } else {
        if (i == eigenvectors.size() - 1) {
          sum += eigenvectors[i][j] * eigenvectors[0][j];
        } else {
          sum += eigenvectors[i][j] * eigenvectors[i + 1][j];
        }
        aux_n -= 1;
      }
      j++;
    }
    bond_order.push_back(sum);
  }
  return bond_order;
}

void plotVector(const vector<double> &v) {
  // Using gnuplot to plot the vector
  // Suavized curve with points
  FILE *pipe = popen("gnuplot -persist", "w");
  if (pipe) {
    fprintf(pipe, "plot '-' with points title 'Points', '-' smooth csplines title 'Levels'\n");
    for (int i = 0; i < v.size(); ++i) {
      fprintf(pipe, "%d %f\n", i, v[i]);
    }
    fprintf(pipe, "e\n");
    for (int i = 0; i < v.size(); ++i) {
      fprintf(pipe, "%d %f\n", i, v[i]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    pclose(pipe);
  }
}

int main() {
  auto a_0 = 0.0;
  auto b_0 = -2.5; // eV

  auto h_n = 1.5;
  auto h_c = 0.0;

  auto k_cc = 1.0;
  // auto k_cc2 = 1.1;
  // auto k_nc1 = 0.8;
  auto k_nc = 1.0;

  auto a_n = a_0 + h_n * b_0;
  auto a_c = a_0 + h_c * b_0;

  // auto b_nc1 = k_nc1 * b_0;
  auto b_nc = k_nc * b_0;
  auto b_cc = k_cc * b_0;
  // auto b_cc2 = k_cc2 * b_0;

  Matrix huckel_matrix = {// {a_n, b_nc1, 0, 0, 0, 0, b_nc1},
                          {a_n, b_nc, 0, 0, 0, 0, b_nc},
                          // {b_nc1, a_c, b_cc2, 0, 0, 0, 0},
                          {b_nc, a_c, b_cc, 0, 0, 0, 0},
                          // {0, b_cc2, a_c, b_nc1, 0, 0, 0},
                          {0, b_cc, a_c, b_nc, 0, 0, 0},
                          {0, 0, b_nc, a_n, b_nc, 0, 0},
                          {0, 0, 0, b_nc, a_c, b_cc, 0},
                          {0, 0, 0, 0, b_cc, a_c, b_cc},
                          {b_nc, 0, 0, 0, 0, b_cc, a_c}};

  cout << "Matriz de Huckel inicial:" << endl;
  printMatrix(huckel_matrix);
  // printMatrixLatexStyle(huckel_matrix);

  // Iteração do método QR
  Matrix Q, R;
  Matrix A = huckel_matrix;
  Matrix eigenvectors(huckel_matrix.size(),
                      vector<double>(huckel_matrix.size(), 0.0));
  for (auto i = 0; i < huckel_matrix.size(); ++i) {
    eigenvectors[i][i] =
        1.0; // Matriz identidade inicial para acumular os autovetores
  }

  const int max_iterations = 7000;
  const double tolerance = 1e-9;

  for (int iter = 0; iter < max_iterations; ++iter) {
    qrDecomposition(A, Q, R);
    A = multiplyMatrices(R, Q);
    eigenvectors = multiplyMatrices(eigenvectors, Q);

    // Convergência: verificar os elementos fora da diagonal
    double off_diagonal_norm = 0.0;
    for (auto i = 0; i < A.size(); ++i) {
      for (auto j = 0; j < A.size(); ++j) {
        if (i != j) {
          off_diagonal_norm += A[i][j] * A[i][j];
        }
      }
    }
    if (sqrt(off_diagonal_norm) < tolerance) {
      break;
    }
  }

  // Armazenar os autovalores em um vetor
  vector<double> eigenvalues;
  for (auto i = 0; i < A.size(); ++i) {
    eigenvalues.push_back(A[i][i]);
  }

  // Ordenar os autovalores e autovetores
  sortEigenvaluesAndEigenvectors(eigenvalues, eigenvectors);

  vector<double> population = eletronicPopulationTotal(eigenvectors, 8);
  vector<double> bond_order = bondOrderCyclic(eigenvectors, 8);

  cout << "Autovalores:" << endl;
  for (auto eigenvalue : eigenvalues) {
    cout << fixed << setprecision(6) << eigenvalue << "\n";
  }

  cout << "Autovetores:" << endl;
  printMatrix(eigenvectors);

  cout << "População eletrônica total:" << endl;
  for (auto pop : population) {
    cout << fixed << setprecision(6) << pop << "\n";
  }

  cout << "Ordem de ligação:" << endl;
  for (auto i = 0; i < bond_order.size(); ++i) {
        if (i == bond_order.size() - 1) {
        cout << fixed << setprecision(6) << "P_{" << i + 1 << "," << 1
             << "} = " << bond_order[i] << "\n";
      } else {
        cout << fixed << setprecision(6) << "P_{" << i + 1 << "," << i + 2
             << "} = " << bond_order[i] << "\n";
      }
  }

  // printEnergyLevelsDiagramLatex(eigenvalues);

  // cout << "Autovalores (LaTeX):" << endl;
  // cout << "$$\\begin{vmatrix}" << endl;
  // for (auto eigenvalue : eigenvalues) {
  //   cout << fixed << setprecision(6) << eigenvalue << " \\\\ ";
  // }
  // cout << "\\end{vmatrix}$$" << endl;

  // cout << "Autovetores (LaTeX):" << endl;
  // printMatrixLatexStyle(eigenvectors);

  // cout << "Orbitais (LaTeX):" << endl;
  // printOrbitalPsiLatex(eigenvectors);

  // cout << "População eletrônica total (LaTeX):" << endl;
  // cout << "$$";
  // for (auto i = 0; i < population.size(); ++i) {
  //   cout << fixed << setprecision(6) << "PET_{" << i + 1
  //        << "} = " << population[i] << " \\\\ ";
  // }
  // cout << "$$" << endl;

  // cout << "Ordem de ligação (LaTeX):" << endl;
  // cout << "$$";
  // for (auto i = 0; i < bond_order.size(); ++i) {
  //   if (i == bond_order.size() - 1) {
  //     cout << fixed << setprecision(6) << "P_{" << i + 1 << "," << 1
  //          << "} = " << bond_order[i] << "\\\\";
  //   } else {
  //     cout << fixed << setprecision(6) << "P_{" << i + 1 << "," << i + 2
  //          << "} = " << bond_order[i] << "\\\\";
  //   }
  // }
  // cout << "$$" << endl;

  // Matrix eigenvector_transpose = transpose(eigenvectors);

  // for (int i = 0; i < eigenvector_transpose.size(); ++i) {
  //   plotVector(eigenvector_transpose[i]);
  // }

  // MatrixXd A_eigen(7, 7);
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     A_eigen(i, j) = huckel_matrix[i][j];
  //   }
  // }

  // cout << "Matriz de Huckel (Eigen):" << endl;
  // cout << A_eigen << endl;

  // cout << "Solução usando Biblioteca Eigen:" << endl;
  // // Resolver autovalores e autovetores
  // Eigen::SelfAdjointEigenSolver<MatrixXd> solver_1(A_eigen);

  // VectorXd omega_squared = solver_1.eigenvalues(); // ω²
  // MatrixXd displacements = solver_1.eigenvectors();

  // cout << "Autovalores (ω²):" << endl;
  // cout << omega_squared << endl;

  // cout << "Autovetores (deslocamentos):" << endl;
  // cout << displacements << endl;

  return 0;
}
