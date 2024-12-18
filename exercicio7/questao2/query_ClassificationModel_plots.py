from mp_api.client import MPRester
from itertools import product
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import os


# Substitua com sua chave API válida do Materials Project
API_KEY = "bMpDCd5C2QI8845pHRTxb4JMk496t9Ru"

# Defina os elementos comuns para A e B em perovskitas
A_elements = ["Ca", "Sr", "Ba", "K", "Na", "Pb"]
B_elements = ["Ti", "Zr", "Nb", "Ta", "Mn", "Fe"]
X_element = "O"  # Oxigênio como o elemento comum X

def query_materials_for_classification():
    """
    Consulta o Materials Project para obter dados de band_gap e estrutura cristalina.
    Classifica os materiais em dois grupos: Metálicos (band_gap == 0) e Semicondutores (0.5 <= band_gap <= 3).
    """
    embeddings = []
    labels = []  # 0 para Metálicos, 1 para Semicondutores

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
                if material.band_gap is not None and material.structure is not None:
                    # Classificação por band_gap
                    if material.band_gap == 0:
                        label = 0  # Metálico
                    elif 0.5 <= material.band_gap <= 3:
                        label = 1  # Semicondutor
                    else:
                        continue  # Ignora materiais fora dos critérios

                    # Gera um embedding fixo baseado na estrutura cristalina
                    embedding = create_embedding(material.structure)
                    embeddings.append(embedding)
                    labels.append(label)

            print(f"Total de materiais processados: {len(embeddings)}")
            return np.array(embeddings), np.array(labels)

        except Exception as e:
            print(f"Erro durante a consulta: {e}")
            return np.array(embeddings), np.array(labels)

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

def train_and_evaluate_classifier(embeddings, labels):
    """
    Treina e avalia um modelo de classificação para prever o tipo de material.
    """
    # Divisão dos dados
    X_train, X_test, y_train, y_test = train_test_split(embeddings, labels, test_size=0.2, random_state=42)

    # Modelo de classificação
    model = RandomForestClassifier(random_state=42)
    model.fit(X_train, y_train)

    # Avaliação
    y_pred = model.predict(X_test)
    print("Relatório de Classificação:")
    print(classification_report(y_test, y_pred))
    print("Matriz de Confusão:")
    print(confusion_matrix(y_test, y_pred))

    return model

# Executa o pipeline completo
embeddings, labels = query_materials_for_classification()
if len(embeddings) > 0:
    model = train_and_evaluate_classifier(embeddings, labels)


import matplotlib.pyplot as plt
import os

def plot_histogram_labels(labels, filename="histogram_labels.png"):
    """
    Gera e salva um histograma dos rótulos (labels) de classificação.
    """
    plt.figure(figsize=(8, 5))
    plt.hist(labels, bins=3, edgecolor="k", alpha=0.7)
    plt.title("Histograma dos Rótulos de Classificação")
    plt.xlabel("Rótulo (0: Metálico, 1: Semicondutor)")
    plt.ylabel("Frequência")
    plt.grid(True)
    plt.savefig(filename)  # Salva o gráfico
    plt.close()  # Fecha o gráfico

def plot_confusion_matrix(cm, classes, filename="confusion_matrix.png"):
    """
    Gera e salva a matriz de confusão.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
    plt.title("Matriz de Confusão")
    plt.colorbar()
    tick_marks = range(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    thresh = cm.max() / 2.0
    for i, j in np.ndindex(cm.shape):
        plt.text(j, i, f"{cm[i, j]}", horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel("Rótulo Verdadeiro")
    plt.xlabel("Rótulo Predito")
    plt.tight_layout()
    plt.savefig(filename)  # Salva o gráfico
    plt.close()  # Fecha o gráfico

# Ajuste na função de treinamento para salvar os gráficos
def train_and_evaluate_classifier_with_plots(embeddings, labels):
    """
    Treina e avalia um modelo de classificação para prever o tipo de material.
    Gera e salva gráficos de avaliação.
    """
    # Divisão dos dados
    X_train, X_test, y_train, y_test = train_test_split(embeddings, labels, test_size=0.2, random_state=42)

    # Modelo de classificação
    model = RandomForestClassifier(random_state=42)
    model.fit(X_train, y_train)

    # Avaliação
    y_pred = model.predict(X_test)
    print("Relatório de Classificação:")
    print(classification_report(y_test, y_pred))

    # Gera a matriz de confusão
    cm = confusion_matrix(y_test, y_pred)
    print("Matriz de Confusão:")
    print(cm)

    # Salva os gráficos
    os.makedirs("plots", exist_ok=True)
    plot_histogram_labels(labels, filename="plots/histogram_labels.png")
    plot_confusion_matrix(cm, classes=["Metálico", "Semicondutor"], filename="plots/confusion_matrix.png")

    return model

# Executa o pipeline completo com os novos gráficos
embeddings, labels = query_materials_for_classification()
if len(embeddings) > 0:
    model = train_and_evaluate_classifier_with_plots(embeddings, labels)

