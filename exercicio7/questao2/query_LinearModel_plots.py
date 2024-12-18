from mp_api.client import MPRester
from itertools import product
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import os
from scipy import stats

# Substitua com sua chave API válida do Materials Project
API_KEY = "bMpDCd5C2QI8845pHRTxb4JMk496t9Ru"

# Defina os elementos comuns para A e B em perovskitas
A_elements = ["Ca", "Sr", "Ba", "K", "Na", "Pb"]
B_elements = ["Ti", "Zr", "Nb", "Ta", "Mn", "Fe"]
X_element = "O"  # Oxigênio como o elemento comum X

def query_materials_filtered():
    """
    Consulta o Materials Project para obter dados de band_gap e estrutura cristalina.
    Exclui materiais com band_gap igual a 0.
    """
    embeddings = []
    band_gaps = []

    with MPRester(API_KEY) as mpr:
        try:
            # Gera todas as combinações químicas possíveis no formato "A-B-O"
            chemsys_filter = [f"{a}-{b}-{X_element}" for a, b in product(A_elements, B_elements)]

            # Realiza a consulta
            perovskites = mpr.materials.summary.search(
                chemsys=chemsys_filter,
                fields=["material_id", "formula_pretty", "band_gap", "structure"],  # Campos a serem recuperados
            )

            # Processa cada material encontrado
            for material in perovskites:
                if material.band_gap is not None and material.band_gap > 1.0 and material.structure is not None:
                    # Gera um embedding fixo baseado na estrutura cristalina
                    embedding = create_embedding(material.structure)
                    embeddings.append(embedding)
                    band_gaps.append(material.band_gap)

            print(f"Total de materiais processados (band_gap > 1.0): {len(embeddings)}")
            return np.array(embeddings), np.array(band_gaps)

        except Exception as e:
            print(f"Erro durante a consulta: {e}")
            return np.array(embeddings), np.array(band_gaps)

def create_embedding(structure):
    """
    Gera um embedding fixo para uma estrutura cristalina.
    Combina estatísticas das posições atômicas e os tipos de átomos.
    """
    species = [site.specie.Z for site in structure]  # Número atômico
    coords = np.array([site.frac_coords for site in structure])  # Coordenadas fracionárias

    # Estatísticas das coordenadas
    mean_coords = np.mean(coords, axis=0)
    std_coords = np.std(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    min_coords = np.min(coords, axis=0)

    # Estatísticas dos tipos atômicos
    mean_species = np.mean(species)
    std_species = np.std(species)

    # Concatena todas as features em um vetor fixo
    embedding = np.concatenate([mean_coords, std_coords, max_coords, min_coords, [mean_species, std_species]])
    return embedding

# Incluindo as funções de visualização

def plot_histogram_band_gaps(band_gaps, filename="histogram_band_gaps.png"):
    """
    Gera e salva um histograma dos valores de band_gap.
    """
    plt.figure(figsize=(8, 5))
    plt.hist(band_gaps, bins=20, edgecolor="k", alpha=0.7)
    plt.title("Histograma dos valores de band_gap")
    plt.xlabel("band_gap (eV)")
    plt.ylabel("Frequência")
    plt.grid(True)
    plt.savefig(filename)  # Salva o gráfico
    plt.close()  # Fecha o gráfico para evitar sobreposição

def plot_embedding_histograms(embeddings, output_dir="embedding_histograms"):
    """
    Gera e salva histogramas das features do embedding.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)  # Cria o diretório para salvar os gráficos

    num_features = embeddings.shape[1]
    for i in range(num_features):
        plt.figure(figsize=(5, 4))
        plt.hist(embeddings[:, i], bins=20, edgecolor="k", alpha=0.7)
        plt.title(f"Feature {i + 1}")
        plt.xlabel("Valor")
        plt.ylabel("Frequência")
        filename = f"{output_dir}/feature_{i + 1}.png"
        plt.savefig(filename)  # Salva o gráfico
        plt.close()  # Fecha o gráfico

def plot_prediction_vs_actual(y_test, y_pred, filename="prediction_vs_actual.png"):
    """
    Gera e salva um gráfico de predição vs valores reais.
    """
    plt.figure(figsize=(8, 8))
    plt.scatter(y_test, y_pred, alpha=0.7, edgecolors="k")
    plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], color="r", linestyle="--", linewidth=2)
    plt.title("Predição vs Valor Real")
    plt.xlabel("Band_gap Real (eV)")
    plt.ylabel("Band_gap Predito (eV)")
    plt.grid(True)
    plt.savefig(filename)  # Salva o gráfico
    plt.close()  # Fecha o gráfico

def train_and_evaluate(embeddings, band_gaps):
    """
    Treina e avalia um modelo linear para prever o band_gap a partir dos embeddings.
    """
    # Divisão dos dados
    X_train, X_test, y_train, y_test = train_test_split(embeddings, band_gaps, test_size=0.2, random_state=42)

    # Modelo linear
    model = LinearRegression()
    model.fit(X_train, y_train)

    # Avaliação
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    print(f"Mean Squared Error: {mse:.4f}")

    # Gera os gráficos
    os.makedirs("plots", exist_ok=True)
    plot_histogram_band_gaps(band_gaps, filename="plots/histogram_band_gaps.png")
    plot_prediction_vs_actual(y_test, y_pred, filename="plots/prediction_vs_actual.png")

    return model

# Executa o pipeline completo
embeddings, band_gaps = query_materials_filtered()
if len(embeddings) > 0:
    model = train_and_evaluate(embeddings, band_gaps)
