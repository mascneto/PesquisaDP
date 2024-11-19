import os
from pathlib import Path
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter, find_peaks

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


for arquivo in arquivos: # aquisições de um arquivo    
    # Load .mat file
    mat_data = scipy.io.loadmat(arquivo)
    
    # Assuming Ch1, Ch2, Ch3, Ch4 are in the loaded mat file
    # You might need to adjust these based on exact variable names in your .mat file
    Ch1 = mat_data['Ch1']
    Ch2 = mat_data['Ch2']
    Ch3 = mat_data['Ch3']
    Ch4 = mat_data['Ch4']
    
    # Extract data from inicio to end of vector
    ampldetec = Ch3[corte_seno_mm_tempo(tamos,Ch2)[2]:corte_seno_mm_tempo(tamos,Ch2)[3]]
    condi_sinal.append(ampldetec)
    
    antena_pura = Ch1[corte_seno_mm_tempo(tamos,Ch2)[2]:corte_seno_mm_tempo(tamos,Ch2)[3]]
    antena_antena.append(antena_pura)
    
    senoide_16 = Ch2[corte_seno_mm_tempo(tamos,Ch2)[2]:corte_seno_mm_tempo(tamos,Ch2)[3]]
    tensao_barra.append(senoide_16)
    
    charge = Ch4[corte_seno_mm_tempo(tamos,Ch2)[2]:corte_seno_mm_tempo(tamos,Ch2)[3]]
    charge_ldic.append(charge)

# Convert lists to numpy arrays
condi_sinal = np.sum(np.column_stack(condi_sinal), axis=1)
antena_antena = np.sum(np.column_stack(antena_antena), axis=1)
tensao_barra = np.sum(np.column_stack(tensao_barra), axis=1)
charge_ldic = np.sum(np.column_stack(charge_ldic), axis=1)

# Assume the following variables are available
tempo = np.arange(1, len(senoide_16)+1)
senoide = tensao_barra
antena = condi_sinal

# Calculate statistics
# media = np.mean(antena, axis=0)
# moda = stats.mode(antena, axis=0)[0]

# Process the data
quadrado = np.abs(antena)
maximo_teste = np.max(quadrado, axis=0)
maxi_maxi = np.max(maximo_teste)
pulso29 = np.zeros_like(antena)

limiar = 30 # or 2.56e-6

for i in range(len(tempo)):
    if quadrado[i] > (limiar/100) * maxi_maxi:
        pulso29[i] = quadrado[i]
    else:
        pass

pulsoreal29 = np.where(pulso29 != 0, antena, 0)

amplitudes = []
fases = []
pks, locs = find_peaks(np.abs(pulsoreal29), distance=3440)
fases.append(locs)
amplitudes.append(pks)

# Plot the data
plt.figure(1)
plt.plot(tempo, senoide*0.1, label='senoide')
plt.plot(tempo, antena, label='antena')
plt.plot(tempo, pulsoreal29, label='pulsoreal29')
plt.legend()
plt.show()