import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io

# criar variavel de caminho do arquivo (esperando almir enviar dados)

# Caminho para os arquivos
caminho = 'E:/Almir Carlos/barra/'  # Altere conforme o seu sistema operacional


for cont in range(49):
    # Definir o nome dos arquivos
    arquivo = f'barra_6kv_{cont}'
    arq_cond = f'barra_6kv_condi_{cont}'
    arq_charge = f'barra_6kv_charge_{cont}'

    # f-strings para concatenar strings e variáveis
    medicao = f"barra_6kv_{cont}/"  # Adiciona o num. da contagem ao dado da barra
    medicao_path = os.path.join('barra_6kv', f"{cont}") # lidar com o caminho do diretório
    
    print("Caminho com f-string:", medicao)
    print("Caminho usando os.path.join:", medicao_path)

# Variaveis
nn=0
condi_sinal = []
tensao_barra = []
antena_antena = []
charge_ldic = []
antena_wave = []
condi_wave = []

