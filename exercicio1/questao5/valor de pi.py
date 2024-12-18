import random
import matplotlib.pyplot as plt

def estimate_pi(num_iterations):
    num_points_inside_circle = 0
    total_num_points = 0
    pi_estimates = []

    for i in range(num_iterations):
        # Gerar coordenadas aleatórias dentro do quadrado unitário
        x = random.uniform(0, 1)
        y = random.uniform(0, 1)

        # Calcular a distância em relação ao centro do quadrado unitário
        distance = x**2 + y**2

        # Verificar se a distância está dentro do círculo unitário
        if distance <= 1:
            num_points_inside_circle += 1

        total_num_points += 1

        # Calcular a estimativa de π até o momento
        pi_estimate = 4 * (num_points_inside_circle / total_num_points)
        pi_estimates.append(pi_estimate)

        # Plotar o valor estimado de π até o momento
        plt.plot(i+1, pi_estimate, 'b.')

    return pi_estimate

# Definir o número de iterações para o cálculo
num_iterations = 1000

# Calcular a estimativa de π e plotar os valores
pi_estimate = estimate_pi(num_iterations)

# Adicionar a linha pontilhada com o valor real de π
plt.axhline(y=3.14159, color='r', linestyle='--', label='Valor real de π')

# Configurar o gráfico
plt.xlabel('Número de Tentativas')
plt.ylabel('Estimativa de π')
plt.title('Estimativa de π em Função do Número de Tentativas')
plt.legend()
plt.grid(True)

# Mostrar o gráfico
plt.show()

