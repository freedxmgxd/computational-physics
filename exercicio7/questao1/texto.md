<!-- Usando o método QR, se possível*

Apropriadamente numere os átomos de sua molécula e nos conte o contexto de obtenção dessa molécula ou onde está mais presente. Usando o método de Huckel (apenas para os orbitias pz), construa o determinante secular, calcule os níveis de energia, faça o diagrama de níveis, faça o preenchimento eletrônico destes níveis, calcule e escreva todos os orbitais moleculares. 

Calcule para cada ligação (apresente sobre a ligação) as ordens de ligações e as populações eletrônicas. Faz sentido para ti? Analise seus dados sob a luz da química?

Esboce os orbitais moleculares HOMO e LUMO. 

Sua molécula é o 1,4Diazepina

File:1,4-diazepine.png - Wikimedia Commons

*(caso não seja adequado, use alguma biblioteca e compare com os resultados do QR)

Vou te guiar no processo de análise da 1,4-diazepina pelo método de Hückel:

    Numeração dos átomos:
    A 1,4-diazepina tem um anel de sete membros com dois nitrogênios nas posições 1 e 4. A numeração segue no sentido horário:
        N1, C2, C3, N4, C5, C6, C7.

    Determinante secular:
    Considerando apenas os orbitais pzpz​, a matriz de Hückel será construída com:
        αα (energia de Coulomb) nos elementos diagonais.
        ββ (energia de ressonância) para interações entre átomos vizinhos.

    Cálculo dos níveis de energia:
    A partir do determinante secular, resolveremos a equação característica para encontrar os valores de energia.

    Diagrama de níveis e preenchimento eletrônico:
    O diagrama mostrará os níveis de energia com preenchimento conforme o número de elétrons ππ da molécula.

    Orbitais moleculares:
    Calcularemos as combinações lineares dos orbitais atômicos pzpz​ para determinar os orbitais moleculares.

    Ordens de ligação e populações eletrônicas:
    Com os coeficientes dos orbitais, calcularemos a densidade eletrônica e as ordens de ligação.

    Orbitais HOMO e LUMO:
    Esboçaremos os orbitais de maior energia ocupada (HOMO) e menor energia desocupada (LUMO).

Se quiser, posso começar construindo o determinante secular ou avançar direto para os cálculos! -->

# Avaliação 7

A 1,4 - Diazepina é um dos principais elementos das benzodiazepinas e tienodiazepinas, componentes esses de algumas drogas farmacêuticas e recreativas.
Para a molecula de 1,4 Diazepina foi numerado seguindo a seguinte ordem,
  1 : N,
  2 : C,
  3 : C,
  4 : N,
  5 : C,
  6 : C,
  7 : C
Como a Hamiltoniana pode ser expressa em termos de soma de unidade de eletrons:
$$ \hat{H} = \sum_{i}^{} \hat{h_i} $$

isso segue que energia total é a soma das energias de cada eletron. No caso de hidrocarbonetos fechados, como a 1,4 Diazepina, cada orbital molecular é ocupado duplamente, o que da o dobro da soma de energia em cada orbital molecular:

$$ E = 2 \sum_{i=1}^{m} \epsilon_i $$

Onde o termo $\epsilon$ é a energia orbital, temos tambem que cada orbital molecular, $\Psi_i$, é descrito como a combinação linear dos orbitais atomicos, $\phi_\mu$:

$$ \Psi_i = \sum_{\mu=1}^{n} C_{\mu i} \phi_\mu $$

Onde $n$ é o numero de atomos do sistema conjugado.

## Determinante Secular da 1,4 Diazepina

Primeiro a formação da matriz de energia da 1,4 Diazepina é dada por:

$$ \begin{vmatrix}
H_{11} - ES_{11} & H_{12} - ES_{12} & H_{13} - ES_{13} & H_{14} - ES_{14} & H_{15} - ES_{15} & H_{16} - ES_{16} & H_{17} - ES_{17} \\
H_{21} - ES_{21} & H_{22} - ES_{22} & H_{23} - ES_{23} & H_{24} - ES_{24} & H_{25} - ES_{25} & H_{26} - ES_{26} & H_{27} - ES_{27} \\
H_{31} - ES_{31} & H_{32} - ES_{32} & H_{33} - ES_{33} & H_{34} - ES_{34} & H_{35} - ES_{35} & H_{36} - ES_{36} & H_{37} - ES_{37} \\
H_{41} - ES_{41} & H_{42} - ES_{42} & H_{43} - ES_{43} & H_{44} - ES_{44} & H_{45} - ES_{45} & H_{46} - ES_{46} & H_{47} - ES_{47} \\
H_{51} - ES_{51} & H_{52} - ES_{52} & H_{53} - ES_{53} & H_{54} - ES_{54} & H_{55} - ES_{55} & H_{56} - ES_{56} & H_{57} - ES_{57} \\
H_{61} - ES_{61} & H_{62} - ES_{62} & H_{63} - ES_{63} & H_{64} - ES_{64} & H_{65} - ES_{65} & H_{66} - ES_{66} & H_{67} - ES_{67} \\
H_{71} - ES_{71} & H_{72} - ES_{72} & H_{73} - ES_{73} & H_{74} - ES_{74} & H_{75} - ES_{75} & H_{76} - ES_{76} & H_{77} - ES_{77} \\
\end{vmatrix} = 0 $$

Aplicando o Método das Ligações fortes temos as aproximações:

$$ H_{ii} = \alpha_i $$

Se os atomos $i$ e $j$ são vizinhos, então:
$$ H_{ij} = \beta_{ij} $$

Se os atomos $i$ e $j$ não são vizinhos, então:
$$ H_{ij} = 0 $$

$$ S_{ij} = \delta_{ij} $$

Sendo $\delta_{ij}$ o delta de Kronecker, que é 1 se $i = j$ e 0 caso contrário.

O Método de Huckel é uma simplificação do Método das Ligações Fortes, onde:

$$ \alpha_i = cte $$

e para dois orbitais atômicos $p_z$ que formem uma ligaçao $\pi$:

$$ \beta_{ij} = cte $$

Onde temos,

$$ \alpha_{x} = \alpha_{0} + h_{x}\beta_{0} = \alpha_{c} + h_{x}\beta_{cc} $$
$$ \beta_{xx'} = \beta_{0} k_{xx'} = \beta_{cc} k_{xx'} $$

Para as ligações pertencentes a 1,4 Diazepina temos:

$$ \alpha_{0} = 0.0 $$
$$ \beta_{0} = -2.5 eV $$
$$ \alpha_{N} = 0.0 + 1.5 \times -2.5 = -3.75 eV $$
$$ \alpha_{C} = 0.0 + 0.0 \times -2.5 = 0.0 eV $$
$$ \beta_{CC} = 1.0 \times -2.5 = -2.5 eV $$
$$ \beta_{NC} = 1.0 \times -2.5 = -2.5 eV $$

A matrix apos a aplicação do método das ligações fortes é:

$$ \begin{vmatrix}
\alpha_{N} & \beta_{NC} & 0 & 0 & 0 & 0 & \beta_{NC} \\
\beta_{NC} & \alpha_{C} & \beta_{CC} & 0 & 0 & 0 & 0 \\
0 & \beta_{CC} & \alpha_{C} & \beta_{NC} & 0 & 0 & 0 \\
0 & 0 & \beta_{NC} & \alpha_{N} & \beta_{NC} & 0 & 0 \\
0 & 0 & 0 & \beta_{NC} & \alpha_{C} & \beta_{CC} & 0 \\
0 & 0 & 0 & 0 & \beta_{CC} & \alpha_{C} & \beta_{NC} \\
\beta_{NC} & 0 & 0 & 0 & 0 & \beta_{NC} & \alpha_{N} \\
\end{vmatrix} 
-
\begin{vmatrix}
E & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & E & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & E & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & E & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & E & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & E & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & E \\
\end{vmatrix} = 0 $$


Para a 1,4 Diazepina, temos a seguinte matriz de Huckel:

$$\begin{vmatrix}
-3.750 & -2.500 & 0.000 & 0.000 & 0.000 & 0.000 & -2.500 & \\
-2.500 & 0.000 & -2.500 & 0.000 & 0.000 & 0.000 & 0.000 & \\
0.000 & -2.500 & 0.000 & -2.500 & 0.000 & 0.000 & 0.000 & \\
0.000 & 0.000 & -2.500 & -3.750 & -2.500 & 0.000 & 0.000 & \\
0.000 & 0.000 & 0.000 & -2.500 & 0.000 & -2.500 & 0.000 & \\
0.000 & 0.000 & 0.000 & 0.000 & -2.500 & 0.000 & -2.500 & \\
-2.500 & 0.000 & 0.000 & 0.000 & 0.000 & -2.500 & 0.000 & \\
\end{vmatrix} = E\bold{I}$$

Para calcular os niveis de energia, precisamos encontrar os autovalores da matriz de Huckel. Para isso, usaremos o método QR.

Onde temos os autovalores, $E_i$:

$$\begin{vmatrix}
-6.601421 \\ -5.629135 \\ -3.214767 \\ -0.457460 \\ 0.768840 \\ 3.610295 \\ 4.023648 \\ \end{vmatrix}$$

Para a 1,4 Diazepina, temos 8 eletrons, 1 de cada Carbono, 2 do Nitrogenio n1, e 1 do Nitrogenio n4, $\pi$ para preencher os orbitais moleculares. O preenchimento eletronico é feito de acordo com o principio de Aufbau, de modo que os orbitais de menor energia são preenchidos primeiro, sendo dois eletrons por orbital, devido ao principio de Pauli, e com spins opostos, de acordo com o principio de Hund.

Então, temos que os niveis de energia ocupados são:
- $E_1 = -6.601 eV$
- $E_2 = -5.629 eV$
- $E_3 = - 3.215 eV$
- $E_4 = -0.457 eV$

e os autovetores:

$$\begin{vmatrix}
-0.536919 & 0.622134 & 0.106998 & 0.406630 & 0.191329 & -0.276301 & 0.186874 & \\
-0.327276 & 0.191329 & 0.374239 & -0.497702 & 0.276301 & 0.622134 & -0.071614 & \\
-0.327276 & -0.191329 & 0.374239 & -0.497702 & -0.276301 & -0.622134 & -0.071614 & \\
-0.536919 & -0.622134 & 0.106998 & 0.406630 & -0.191329 & 0.276301 & 0.186874 & \\
-0.285116 & -0.276301 & -0.397147 & -0.037837 & 0.622134 & -0.191329 & -0.509462 & \\
-0.215951 & -0.000000 & -0.617691 & -0.413554 & 0.000000 & 0.000000 & 0.633085 & \\
-0.285116 & 0.276301 & -0.397147 & -0.037837 & -0.622134 & 0.191329 & -0.509462 & \\
\end{vmatrix}$$

Iremos agora calcular os orbitais moleculares, as ordens de ligação e as populações eletrônicas.

## Orbitais Moleculares

Os orbitais moleculares são obtidos a partir dos autovetores, que são combinações lineares dos orbitais atômicos. Para a 1,4 Diazepina, temos 7 orbitais moleculares, $\Psi_i$, que são combinações lineares dos orbitais atômicos, $\phi_\mu$.

$$
\psi_{1}  = (-0.537) \cdot \phi_{1} + (-0.327) \cdot \phi_{2} + (-0.327) \cdot \phi_{3} + (-0.537) \cdot \phi_{4} + (-0.285) \cdot \phi_{5} + (-0.216) \cdot \phi_{6} + (-0.285) \cdot \phi_{7}\\
\psi_{2}  = (0.622) \cdot \phi_{1} + (0.191) \cdot \phi_{2} + (-0.191) \cdot \phi_{3} + (-0.622) \cdot \phi_{4} + (-0.276) \cdot \phi_{5} + (-0.000) \cdot \phi_{6} + (0.276) \cdot \phi_{7}\\
\psi_{3}  = (0.107) \cdot \phi_{1} + (0.374) \cdot \phi_{2} + (0.374) \cdot \phi_{3} + (0.107) \cdot \phi_{4} + (-0.397) \cdot \phi_{5} + (-0.618) \cdot \phi_{6} + (-0.397) \cdot \phi_{7}\\
\psi_{4}  = (0.407) \cdot \phi_{1} + (-0.498) \cdot \phi_{2} + (-0.498) \cdot \phi_{3} + (0.407) \cdot \phi_{4} + (-0.038) \cdot \phi_{5} + (-0.414) \cdot \phi_{6} + (-0.038) \cdot \phi_{7}\\
\psi_{5}  = (0.191) \cdot \phi_{1} + (0.276) \cdot \phi_{2} + (-0.276) \cdot \phi_{3} + (-0.191) \cdot \phi_{4} + (0.622) \cdot \phi_{5} + (0.000) \cdot \phi_{6} + (-0.622) \cdot \phi_{7}\\
\psi_{6}  = (-0.276) \cdot \phi_{1} + (0.622) \cdot \phi_{2} + (-0.622) \cdot \phi_{3} + (0.276) \cdot \phi_{4} + (-0.191) \cdot \phi_{5} + (0.000) \cdot \phi_{6} + (0.191) \cdot \phi_{7}\\
\psi_{7}  = (0.187) \cdot \phi_{1} + (-0.072) \cdot \phi_{2} + (-0.072) \cdot \phi_{3} + (0.187) \cdot \phi_{4} + (-0.509) \cdot \phi_{5} + (0.633) \cdot \phi_{6} + (-0.509) \cdot \phi_{7}\\
$$

Dados os orbitais moleculares, podemos calcular as ordens de ligação e as populações eletrônicas.

Para as populações eletrônicas, temos que a densidade eletrônica é dada por:

$$
PET_{0} = 1.704259 \\ PET_{1} = 1.062956 \\ PET_{2} = 1.062956 \\ PET_{3} = 1.704259 \\ PET_{4} = 0.633581 \\ PET_{5} = 1.198407 \\ PET_{6} = 0.633581 \\ $$

E para as ordens de ligação temos portanto, 

$$P_{1,2} = 0.264830\\P_{2,3} = 0.916530\\P_{3,4} = 0.264830\\P_{4,5} = 0.534203\\P_{5,6} = 0.645065\\P_{6,7} = 0.645065\\P_{7,1} = 0.534203\\$$

### **Análise Química dos Resultados**

1. **Níveis de Energia e Preenchimento Eletrônico:**  
   - Os quatro primeiros níveis de energia ($E_1$ a $E_4$) estão ocupados, totalizando 8 elétrons π.  
   - O **HOMO** (Highest Occupied Molecular Orbital) é o orbital de energia $E_4 = -0.457\, eV$.  
   - O **LUMO** (Lowest Unoccupied Molecular Orbital) é o orbital de energia $E_5 = 0.769\, eV$.  
   - A diferença de energia entre HOMO e LUMO ($\Delta E \approx 1.226\, eV$) indica moderada estabilidade e possível reatividade química.

2. **Ordens de Ligação:**  
   - As ligações C2–C3 ($P_{2,3} = 0.917$) e C5–C6/C6–C7 ($P_{5,6} = P_{6,7} = 0.645$) têm **altas ordens de ligação**, sugerindo regiões de alta densidade eletrônica e estabilidade.  
   - As ligações N1–C2 e C3–N4 ($P_{1,2} = P_{3,4} = 0.265$) são mais fracas, sugerindo **menor densidade eletrônica**, o que pode ser regiões de reatividade.

3. **Populações Eletrônicas:**  
   - **Nitrogênios (N1 e N4):** Têm populações eletrônicas maiores ($\sim1.70$), indicando acumulação de carga negativa.  
   - **Carbonos C5 e C7:** Menores populações eletrônicas ($\sim0.63$), sugerindo maior eletropositividade, o que pode favorecer ataques nucleofílicos.

## Codigo usado para analise:
```C++
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

``` 