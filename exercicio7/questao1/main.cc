#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
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

int main() {
    auto a_0 = 0.0;
    auto b_0 = -2.5; // eV

    auto h_n = 1.5;
    auto h_c = 0.0;

    auto k_cc = 1.0;
    auto k_cc2 = 1.1;
    auto k_nc1 = 0.8;
    auto k_nc = 1.0;

    auto a_n = a_0 + h_n * b_0;
    auto a_c = a_0 + h_c * b_0;

    auto b_nc1 = k_nc1 * b_0;
    auto b_nc = k_nc * b_0;
    auto b_cc = k_cc * b_0;
    auto b_cc2 = k_cc2 * b_0;

    Matrix huckel_matrix = {
        {a_n, b_nc1, 0, 0, 0, 0, b_nc1}, {b_nc1, a_c, b_cc2, 0, 0, 0, 0},
        {0, b_cc2, a_c, b_nc1, 0, 0, 0}, {0, 0, b_nc, a_n, b_nc, 0, 0},
        {0, 0, 0, b_nc, a_c, b_cc, 0},   {0, 0, 0, 0, b_cc, a_c, b_cc},
        {b_nc, 0, 0, 0, 0, b_cc, a_c}};

    cout << "Matriz de Huckel inicial:" << endl;
    printMatrix(huckel_matrix);

    // Iteração do método QR
    Matrix Q, R;
    Matrix A = huckel_matrix;
    Matrix eigenvectors(huckel_matrix.size(), vector<double>(huckel_matrix.size(), 0.0));
    for (auto i = 0; i < huckel_matrix.size(); ++i) {
        eigenvectors[i][i] = 1.0; // Matriz identidade inicial para acumular os autovetores
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

    cout << "Autovalores:" << endl;
    for (auto i = 0; i < A.size(); ++i) {
        cout << fixed << setprecision(6) << A[i][i] << endl;
    }

    cout << "Autovetores:" << endl;
    printMatrix(eigenvectors);

    MatrixXd A_eigen(7, 7);
    for (int i = 0; i < 7; i++) {
      for (int j = 0; j < 7; j++) {
        A_eigen(i, j) = huckel_matrix[i][j];
      }
    }
  
    cout << "Matriz de Huckel (Eigen):" << endl;
    cout << A_eigen << endl;
  
    // Resolver autovalores e autovetores
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver_1(A_eigen);
  
    VectorXd omega_squared = solver_1.eigenvalues(); // ω²
    MatrixXd displacements = solver_1.eigenvectors();
  
    cout << "Autovalores (ω²):" << endl;
    cout << omega_squared << endl;
  
    cout << "Autovetores (deslocamentos):" << endl;
    cout << displacements << endl;
  

    return 0;
}
