import os
import scipy.io
from scipy.signal import lfilter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

caminho = 'C:\\Users\yuutr\\Desktop\\Pasta Pessoal - Boaventura Filho\\PIBIC_Almir\\Dados\\Barra'

# Exibir as variáveis contidas no arquivo .mat
# print(dados.keys())

#Função para renumerar os arquivos
#Descrição das variáveis:
#caminho - Path para a pastas na qual voce quer modificar e enumerar seus dados
#aq - nome que vc quer no arquivo
def renumerar_arquivo(caminho,aq):
    lista_arquivos = os.listdir(caminho)
    print(lista_arquivos)
    for cont in range(len(lista_arquivos)):
        arquivo = f'{aq}_{cont}.mat'
        arq_cond = f'{aq}_condi_{cont}.mat'
        arq_charge = f'{aq}_charge_{cont}.mat'

        caminho_arquivo_antigo = os.path.join(caminho, lista_arquivos[cont]) 
        caminho_arquivo_novo = os.path.join(caminho, arquivo)
        os.rename(caminho_arquivo_antigo, caminho_arquivo_novo)
        
    print(os.listdir(caminho))


# Acessar uma variável específica dentro do arquivo .mat
        

#Ch2 = dados['Ch2']

# Criando o gráfico
#plt.plot(dados.get('Ch1')*1000)
#plt.plot(dados.get('Ch2')*100)
#plt.xlabel("Pontos adquiridos (x10e6)")
#plt.ylabel("Amplitude (mV)")
#plt.plot(dados.get('Ch3'))
#plt.grid(True)
#plt.show()

# Encontrar os máximos e mínimos em Ch2
#maior = np.max(Ch2)
#ind_maior = np.argmax(Ch2)

#menor = np.min(Ch2)
#ind_menor = np.argmin(Ch2)

# print(ind_maior, ind_menor)

# Calcular as variáveis c_o, n_g, zero, início e tamanho do vetor
#c_o = ind_menor - ind_maior
#n_g = c_o / 2 # corte origem ?
#zero = ind_maior - n_g
#inicio = float(zero + 45000)  # Ajuste dependendo da unidade de amostra
#tamanho_vetor = float(ind_menor + n_g + 70000)

#Ch2 = np.array(dados.get('Ch2')) #  Tratamento dos dados de Ch2 para que vire um array numpy
#Ch2_1d = np.ravel(Ch2) # tratamento dos dados de Ch2 para que seja um array numpy unidimensional


from scipy.signal import lfilter

#Processamento para extrair os inputs do cortesen
lista_arquivos = os.listdir(caminho)
print(lista_arquivos)
caminho_tamos = os.path.join(caminho, lista_arquivos[0])
caminho_Ch_sen = os.path.join(caminho,lista_arquivos[1])
dados_tamos = scipy.io.loadmat(caminho_tamos)
dados_Ch = scipy.io.loadmat(caminho_Ch_sen)

Ch_sen = dados_Ch['Ch2']
tamos = dados_tamos['dt']


def corte_seno_mm(Ch_seno):
    Ch_seno = dados_Ch['Ch2']
    Ch_seno = np.array(Ch_seno)
    Ch_seno = np.ravel(Ch_seno)

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

    #Reconstrução do vetor de tempo:
    Qamos = len(y)
    tfinal = tamos*Qamos
            
            
    vetor_t = np.linspace(0,tfinal,Qamos)

    vetor_t = vetor_t - vetor_t[in_idx]
    vetor_t = np.ravel(vetor_t)

    return [in_idx, fim, vetor_t[in_idx:fim], y[in_idx:fim]]


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

        return [vetor_t[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]], y[indices_passagem_zero_janela[0]:indices_passagem_zero_janela[2]],indices_passagem_zero_janela[0], indices_passagem_zero_janela[2]]

plt.plot(corte_seno_mm_tempo(tamos,Ch_sen)[0],corte_seno_mm_tempo(tamos,Ch_sen)[1])
plt.show()

#Acumular os pulsos nas devidas variáveis:
condi_sinal = pd.DataFrame()
tensao_barra = pd.DataFrame()
antena_antena = pd.DataFrame()
charge_ldic = pd.DataFrame()
antena_wave = pd.DataFrame()
condi_wave = pd.DataFrame()

for cont in range(1,len(lista_arquivos)):
    
    caminho_arquivo = os.path.join(caminho, lista_arquivos[cont]) 
    Dados = scipy.io.loadmat(caminho_arquivo)

    #Extraindo os Dados dos Canais dos osciloscopios:
    DfCH1 = pd.DataFrame(Dados['Ch1'])
    DfCH2 = pd.DataFrame(Dados['Ch2'])
    DfCH3 = pd.DataFrame(Dados['Ch3'])
    DfCH4 = pd.DataFrame(Dados['Ch4'])

    #Canal 3
    ampldetec = DfCH3[corte_seno_mm_tempo(tamos,DfCH2)[2]:corte_seno_mm_tempo(tamos,DfCH2)[3]]
    condi_sinal = [condi_sinal,ampldetec]

    #Canal 1
    antena_pura = DfCH1[corte_seno_mm_tempo(tamos,DfCH2)[2]:corte_seno_mm_tempo(tamos,DfCH2)[3]]
    antena_antena = [antena_antena, antena_pura]

    #Canal 4
    charge = DfCH4[corte_seno_mm_tempo(tamos,DfCH2)[2]:corte_seno_mm_tempo(tamos,DfCH2)[3]]
    charge_ldic = [charge_ldic, charge]



