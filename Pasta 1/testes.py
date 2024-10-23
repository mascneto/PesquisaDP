import scipy.io
from scipy.signal import lfilter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dados = scipy.io.loadmat(r"C:\Users\jmlnn\Downloads\barra_6kv_1_1.mat")

# Exibir as variáveis contidas no arquivo .mat
# print(dados.keys())

# Acessar uma variável específica dentro do arquivo .mat
Ch2 = dados['Ch2']

# Criando o gráfico
plt.plot(dados.get('Ch1')*1000)
plt.plot(dados.get('Ch2')*100)
plt.xlabel("Pontos adquiridos (x10e6)")
plt.ylabel("Amplitude (mV)")
#plt.plot(dados.get('Ch3'))
plt.grid(True)
#plt.show()

# Encontrar os máximos e mínimos em Ch2
maior = np.max(Ch2)
ind_maior = np.argmax(Ch2)

menor = np.min(Ch2)
ind_menor = np.argmin(Ch2)

# print(ind_maior, ind_menor)

# Calcular as variáveis c_o, n_g, zero, início e tamanho do vetor
c_o = ind_menor - ind_maior
n_g = c_o / 2 # corte origem ?
zero = ind_maior - n_g
inicio = float(zero + 45000)  # Ajuste dependendo da unidade de amostra
tamanho_vetor = float(ind_menor + n_g + 70000)

Ch2 = np.array(dados.get('Ch2')) #  Tratamento dos dados de Ch2 para que vire um array numpy
Ch2_1d = np.ravel(Ch2) # tratamento dos dados de Ch2 para que seja um array numpy unidimensional


from scipy.signal import lfilter

def corte_seno_mm(Ch_seno):
    """
    Função que determina os zeros de uma senoide adquirida,
    aplicando um filtro de média móvel.
    """

    # Parâmetros do filtro
    windowSize = 5000
    b = np.ones(windowSize) / windowSize  # Coeficientes do numerador
    a = 1  # Coeficiente do denominador

    # Aplica o filtro de média móvel
    y = lfilter(b, a, Ch_seno)

    # Encontra os índices de passagem por zero
    indices_passagem_zero_janela = np.where(np.diff(np.sign(y)) != 0)[0]

    # Calcula o índice de fim e início conforme a lógica
    fim = indices_passagem_zero_janela[2] + 2000
    in_idx = indices_passagem_zero_janela[1] - (indices_passagem_zero_janela[2] - indices_passagem_zero_janela[1]) - 2000

    return print(in_idx, fim)


corte_seno_mm(Ch2_1d)