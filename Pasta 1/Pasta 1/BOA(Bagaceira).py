import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.io import loadmat

# criar variavel de caminho do arquivo (esperando almir enviar dados)

# Caminho para os arquivos
caminho = 'C:\\Users\yuutr\\Desktop\\Pasta Pessoal - Boaventura Filho\\PIBIC_Almir\\Dados\\Barra'  # Altere conforme o seu sistema operacional


#Função para renumerar os arquivos
#Descrição das variáveis:
#caminho - Path para a pastas na qual voce quer modificar e enumerar seus dados
#aq - nome que vc quer no arquivo
def renumerar_arquivo(caminho,aq):
    lista_arquivos = os.listdir(caminho)
    for cont in range(len(lista_arquivos)):
        arquivo = f'{aq}_{cont}.mat'
        arq_cond = f'{aq}_condi_{cont}.mat'
        arq_charge = f'{aq}_charge_{cont}.mat'

        caminho_arquivo_antigo = os.path.join(caminho, lista_arquivos[cont]) 
        caminho_arquivo_novo = os.path.join(caminho, arquivo)
        os.rename(caminho_arquivo_antigo, caminho_arquivo_novo)
        
    return(os.listdir(caminho))

renumerar_arquivo(caminho,'barra_6kv')
