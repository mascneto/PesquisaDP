import os
from pathlib import Path
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter

# Parâmetros
n_pastas = 1
n_arq = 4

#Acumular os pulsos nas devidas variáveis:
nn = 0
condi_sinal = []
tensao_barra = []
antena_antena = []
charge_ldic = []
antena_wave = []
condi_wave = []

# Diretório base contendo as pastas para análise
caminho = r"C:\Users\jmlnn\Documents\PesquisaDP\OneDrive-2024-11-12"
diretorio_base = Path(caminho) # C:\Users\jmlnn\Documents\PesquisaDP\PesquisaDP\OneDrive-2024-11-12

# Inicializar contadores
total_pastas = 0
total_arquivos = 0

# Verificar se o diretório base existe
if not diretorio_base.is_dir():
    raise FileNotFoundError(f"Diretório base não encontrado: {diretorio_base}")

# Listar todas as pastas dentro do diretório base
pastas = [item for item in diretorio_base.iterdir() if item.is_dir()]

# Garantir que temos 50 pastas para iterar (limitar a 50 pastas)
for folder in pastas[:n_pastas]:
    total_pastas += 1
    print(f"Analisando a pasta: {folder.name}")
    
   # Listar todos os arquivos .mat dentro da pasta
    arquivos = [file for file in folder.iterdir() if file.is_file() and file.suffix == '.mat']
    
    # Garantir que temos 50 arquivos para iterar em cada pasta
    for file in arquivos[:n_arq]:
        total_arquivos += 1
        print(f"Processando arquivo: {file.name}")
        mat_data = scipy.io.loadmat(file) # Carregar o arquivo .mat

# Resumo da análise
print(f"Total de pastas analisadas: {total_pastas}")
print(f"Total de arquivos .mat analisados: {total_arquivos}")

# Identificar os caminhos dos arquivos a serem processados
caminho_tamos = arquivos[0]
caminho_Ch_sen = arquivos[1]
dados_tamos = scipy.io.loadmat(caminho_tamos)
dados_Ch = scipy.io.loadmat(caminho_Ch_sen)
Ch_sen = dados_Ch['Ch2']
tamos = dados_tamos['dt']


def corte_seno_mm_tempo(tamos, Ch_sen):
        #Extraindo o seno do .mat:
        Ch_seno = dados_Ch['Ch2']
        Ch_seno = np.array(Ch_seno).ravel()

        #Aplicando o lfilter no sinal:
        #     
        # Parâmetros do filtro
        windowSize = 5000
        b = np.ones(windowSize) / windowSize  # Coeficientes do numerador
        a = 1  # Coeficiente do denominador

        # Aplica o filtro de média móvel
        y = lfilter(b, a, Ch_seno)

        #Reconstrução do vetor de tempo:
        Qamos = len(y) #Quantidade de amostras

        tfinal = tamos*Qamos

        vetor_t = np.linspace(0,tfinal,Qamos)

        indices_passagem_zero_janela = np.where(np.diff(np.sign(y)) != 0)[0] 

        # Calcula o índice de fim e início
        fim = indices_passagem_zero_janela[2] + 2000

        in_idx = indices_passagem_zero_janela[1] - (indices_passagem_zero_janela[2] - indices_passagem_zero_janela[1]) - 2000

        indices_passagem_zero_janela = np.where(np.diff(np.sign(y)) != 0)[0]

        vetor_t = vetor_t - vetor_t[indices_passagem_zero_janela[0]]

        return [vetor_t[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]],
                y[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]],
                indices_passagem_zero_janela[0],
                indices_passagem_zero_janela[2],
                in_idx,
                fim]


for arquivo in arquivos:
    # print(arquivo)
    dados = scipy.io.loadmat(arquivo)

    #Extraindo os Dados dos Canais dos osciloscopios:
    ch1 = np.array(dados['Ch1']).ravel()
    ch2 = np.array(dados['Ch2']).ravel()
    ch3 = np.array(dados['Ch3']).ravel()
    ch4 = np.array(dados['Ch4']).ravel()

    #Canal 3
    ampldetec = ch3[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]] # começo e fim dos pulsos do seno cortado
    condi_sinal = np.append(condi_sinal, ampldetec)

    #Canal 1
    antena_pura = ch1[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]] # pega um arquivo .mat por vez
    antena_antena = np.append(antena_antena, antena_pura) # soma todos os valores das antenas .mat

    #Canal 4
    charge = ch4[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]]
    charge_ldic = np.append(charge_ldic, charge)

antena = condi_sinal # no cod. original, Almir renomeia esta variável
quadrado = abs(antena)
pulso29 = np.zeros_like(quadrado)

#plt.plot(quadrado)
#plt.show()

# print(f"ampldetec: {len(ampldetec)}") # para 1 pasta e 4 arquivos, retorna 2082275
print(ampldetec.shape)

atenuacao = 0.1
print(condi_sinal)


