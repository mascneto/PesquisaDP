import os
import scipy.io
from scipy.signal import lfilter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import zipfile

# caminho = 'C:\\Users\yuutr\\Desktop\\Pasta Pessoal - Boaventura Filho\\PIBIC_Almir\\Dados\\Barra'
# caminho = "C:\Users\jmlnn\Downloads\barra_6kv_1_1.mat"

# Exemplo de uso (nome do arquivo .zip sem o caminho completo)
nome_arquivo_zip = "OneDrive-2024-11-12.zip"

CAMINHO_ARQUIVO_ZIP = os.path.dirname(os.path.abspath(__file__)) 

# print(f"caminho_arquivo_zip: {CAMINHO_ARQUIVO_ZIP}")
# Obtendo o caminho do arquivo .zip da variável de ambiente
caminho_zip = Path(os.getenv('CAMINHO_ARQUIVO_ZIP', 'C:/Users/jmlnn/Documents/PesquisaDP/PesquisaDP/Pasta 1/OneDrive-2024-11-12.zip'))

def renumerar_arquivo(caminho,aq):
    with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
        lista_arquivos = zip_ref.namelist()
        print(lista_arquivos)
    for cont in range(len(lista_arquivos)):
        arquivo = f'{aq}_{cont}.mat'
        arq_cond = f'{aq}_condi_{cont}.mat'
        arq_charge = f'{aq}_charge_{cont}.mat'

        caminho_arquivo_antigo = os.path.join(caminho, lista_arquivos[cont]) 
        caminho_arquivo_novo = os.path.join(caminho, arquivo)
        os.rename(caminho_arquivo_antigo, caminho_arquivo_novo)
        
    print(os.listdir(caminho))


# Abrir o arquivo ZIP
with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
    # Listar os arquivos dentro do ZIP
    lista_arquivos = zip_ref.namelist()
    # print("Arquivos no ZIP:", lista_arquivos)

    caminho_tamos = lista_arquivos[0]
    caminho_Ch_sen = lista_arquivos[1]

    #print(f'caminho tamos:{caminho_tamos}')

    # Verificar se os arquivos foram encontrados
    if caminho_tamos and caminho_Ch_sen:
        # Abrir e carregar os arquivos .mat diretamente do ZIP
        with zip_ref.open(caminho_tamos) as file_tamos, zip_ref.open(caminho_Ch_sen) as file_Ch_sen:
            dados_tamos = scipy.io.loadmat(file_tamos)
            dados_Ch = scipy.io.loadmat(file_Ch_sen)

        #print("Dados Tamos:", dados_tamos)
        #print("Dados Ch:", dados_Ch)
    else:
        pass


Ch_sen = dados_Ch['Ch2']
tamos = dados_tamos['dt']


def corte_seno_mm_tempo(tamos, Ch_sen):
        #Extraindo o seno do .mat:
        Ch_seno = dados_Ch['Ch2']
        Ch_seno = np.array(Ch_seno)
        Ch_seno = np.ravel(Ch_seno)
        #Aplicando o lfilter no sinal:
            
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

        vetor_t = vetor_t - vetor_t[indices_passagem_zero_janela[0]]
        vetor_t = np.ravel(vetor_t)

        return [vetor_t[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]], \
                y[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]], \
                indices_passagem_zero_janela[0], \
                indices_passagem_zero_janela[2]]


plt.plot(corte_seno_mm_tempo(tamos,Ch_sen)[0],corte_seno_mm_tempo(tamos,Ch_sen)[1])
#plt.show()

#Acumular os pulsos nas devidas variáveis:
condi_sinal = []
tensao_barra = []
antena_antena = []
charge_ldic = []
antena_wave = []
condi_wave = []

# Abrir o arquivo ZIP
with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
    # Listar os arquivos dentro do ZIP
    lista_arquivos = zip_ref.namelist()
    # print("Arquivos no ZIP:", lista_arquivos)

    # Loop para processar os arquivos
    for cont in range(0, len(lista_arquivos)):  # Começar de 1, porque 0 pode ser o arquivo ZIP ou arquivo que você não quer processar
        caminho_arquivo = lista_arquivos[cont]  # Caminho do arquivo dentro do ZIP (não é necessário usar os.path.join aqui)

        # Abrir o arquivo dentro do ZIP
        with zip_ref.open(caminho_arquivo) as file:
            # Carregar os dados usando scipy.io.loadmat
            Dados = scipy.io.loadmat(file)
            # print(caminho_arquivo)
            # Extraindo os Dados dos Canais dos Osciloscópio
            ch1 = np.ravel(np.array(Dados['Ch1']))
            ch2 = np.ravel(np.array(Dados['Ch2']))
            ch3 = np.ravel(np.array(Dados['Ch3']))
            ch4 = np.ravel(np.array(Dados['Ch4']))

        #Canal 3
        ampldetec = ch3[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]] # começo e fim dos pulsos do seno cortado
        condi_sinal = [condi_sinal, ampldetec] # condi_sinal é a lista vazia

        #Canal 1
        antena_pura = ch1[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]] # pega um arquivo .mat por vez
        antena_antena = np.append(antena_antena, antena_pura) # soma todos os valores das antenas .mat

        #Canal 4
        charge = ch4[corte_seno_mm_tempo(tamos,ch2)[2]:corte_seno_mm_tempo(tamos,ch2)[3]]
        charge_ldic = [charge_ldic, charge]

aten = 0.1 # atenuação

max_antena = np.max(charge_ldic)
min_antena = np.min(charge_ldic)

"""
% maximo_antena = max(charge_ldic);
% minimo_antena = min(charge_ldic);
% maximo_maximo_antena = max(maximo_antena);
% minimo_minimo_antena = min(minimo_antena);
% media_ruido = mean(charge_ldic,1);
"""

print(f"antena_pura: {antena_pura}")
print(f"antena_antena: {antena_antena}")