<!-- Considere a equação de movimento para um pêndulo simples:

d2θdt2=gLsinθ


onde θ
é o deslocamento angular. A aceleração angular nesse caso e α=−(g/L)sinθ

. Analiticamente, você deve ter resolvido esse problema utilizando a aproximação para oscilações com pequenas amplitudes angulares:

d2θdt2=gLsinθ≈gLθ

 

que possuem solução do tipo:

θ(t)=θ0cos(2πT+ϕ)

sendo o período (nessa aproximação) dado por:

T=2πLg−−√

Resolveremos esse sistema sem fazer a aproximação de ângulos pequenos. Para grandes amplitudes de oscilação, o período sofre correções que dependem da amplitiude θ0

:

T=2πLg−−√(1+116θ20+…)


onde θ0

é o deslocamento angular inicial do pêndulo. A energia total do sistema é:

E=12mL2ω2−mgLcosθ

(a) Faça o fluxograma do programa Euler deste link. Modifique o programa (Euler) de tal forma que este também calcule a solução analítica para θ(t)

(b) Considere θ0=10o,g=10,L=1,dt=0.05,tmax=20
. Compare o resultado analítico com o obtido numericamente para θ(t), utilizando os métodos de Euler e Euler-Cromer. Plote os três resultados no mesmo gráfico θ(t)×t

. Discuta o resultado obtido para a energia, comparando os dois métodos numéricos.

(c) Utilizando os mesmo parâmetros do item anterior, determine numericamente um valor máximo de dt tal que a precisão obtida pelo método de Euler para θ(t=5)

(quando comparado à solução exata) seja de 1%.

(d) Utilize o dt obtido no exercício anterior, utilize o método Euler-Cromer considerando θ0=90o,g=10,L=1
para obter os gráficos θ0×t,ω(t)×teE(t)×t para 0

. Determine o valor do período nesse caso. Compare o resultado com as soluções analíticas correspondentes ao período. 

(e) Obtenha o gráfico T×θ0
, para 10o<θ0<120o utilizando o método de Euler-Cromer. Que conclusão você pode obter sobre a solução analítica a partir deste gráfico? -->

# Avaliação 8 - Questão 2

## a)
<!-- Faça o fluxograma do programa Euler deste link. Modifique o programa (Euler) de tal forma que este também calcule a solução analítica para θ(t) -->

```Mermaid
graph TD
    Start([Início]) --> DefineParams["Inicialização dos parâmetros e condições iniciais"]
    DefineParams --> CheckTime{t_i <= t_max?}
    
    CheckTime -- Sim --> UpdateTime[Atualizar t_f = t_i + dt]
    UpdateTime --> CalcValues["Calcular valores novos valores de theta, omega e energia"]
    CalcValues --> WriteFiles[Escrever nos arquivos]
    WriteFiles --> UpdateVars[Atualizar t_i, theta, omega]
    UpdateVars --> CheckTime

    CheckTime -- Não --> CloseFiles[Fechar arquivos]
    CloseFiles --> End([Encerrar])

```

