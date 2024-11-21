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
tempo = np.arange(1, len(senoide_16)+1)

def process_signal(
    tensao_barra, 
    condi_sinal, 
    tempo, 
    limiar=30, 
    peak_distance=3440, 
    plot=True
):
    """
    Process signal with peak detection and analysis.

    Parameters:
    -----------
    tensao_barra : array
        Input voltage bar signal
    condi_sinal : array
        Input condition signal 
    tempo : array
        Time array
    limiar : float, optional
        Threshold percentage (default 30)
    peak_distance : int, optional
        Minimum distance between peaks (default 3440)
    plot : bool, optional
        Whether to generate plots (default False)

    Returns:
    --------
    dict : Dictionary containing processed signal data
    """
    # Assign signals
    senoide = tensao_barra
    antena = condi_sinal

    # Ensure numpy arrays
    tempo = np.array(tempo)
    antena = np.array(antena)
    senoide = np.array(senoide)

    # Process the data
    quadrado = np.abs(antena)
    maximo_teste = np.max(quadrado, axis=0)
    maxi_maxi = np.max(maximo_teste)

    # Thresholding
    pulso29 = np.zeros_like(antena)
    pulso29[quadrado > (limiar/100) * maxi_maxi] = quadrado[quadrado > (limiar/100) * maxi_maxi]

    # Create pulsoreal29
    pulsoreal29 = np.where(pulso29 != 0, antena, 0)

    # Find peaks
    pks, _ = scipy.signal.find_peaks(np.abs(pulsoreal29), distance=peak_distance)

    # Get peak values
    locs = np.array(pks, dtype=np.int64)
    peak_amplitudes = pulsoreal29[locs]
    peak_times = tempo[locs]

    # Horizontal concatenation of amplitudes and phases
    novo = np.column_stack((peak_amplitudes, peak_times))

    # Create x-axis steps
    passo_plot = 360/len(antena)
    eixo_x = np.arange(passo_plot, 360 + passo_plot, passo_plot)

    # First column of senoide (if 2D)
    senoide_16 = senoide[:, 0] if len(senoide.shape) > 1 else senoide

    # 360-degree axis
    eixo_360 = eixo_x

    # Auxiliary conversion vector
    tam360 = eixo_360.shape[0]
    auxiliar_conversao = (novo[:, 1] / tam360) * 360

    # 2D amplitude vector
    amplitudes_2d = peak_amplitudes.reshape(-1, 1)
    novo_graus = np.hstack([amplitudes_2d, auxiliar_conversao.reshape(-1, 1)])

    # Absolute values
    absoluto = np.abs(novo_graus)

    # Optional plotting
    if plot:
        plt.figure(figsize=(12, 6))
        plt.plot(tempo, pulsoreal29, 'b-', label='Signal', alpha=0.7)
        plt.plot(peak_times, peak_amplitudes, 'ro', label='Peaks', markersize=8)
        plt.plot(tempo, senoide*0.1, label='senoide')
        plt.plot(tempo, antena, label='antena')
        plt.plot(tempo, pulsoreal29, label='pulsoreal29')
        plt.title('Signal with Detected Peaks')
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()

    # Return results in a dictionary
    return {
        'pulso29': pulso29,
        'pulsoreal29': pulsoreal29,
        'peak_amplitudes': peak_amplitudes,
        'peak_times': peak_times,
        'novo': novo,
        'novo_graus': novo_graus,
        'absoluto': absoluto,
        'senoide_16': senoide_16,
        'eixo_360': eixo_360
    }

# Example usage
# result = process_signal(tensao_barra, condi_sinal, tempo, plot=True)

result = process_signal(
    tensao_barra, 
    condi_sinal, 
    tempo, 
    limiar=30,  # optional
    peak_distance=3440,  # optional
    plot=True  # optional, shows plot if True
)

# Access results
peak_amplitudes = result['peak_amplitudes']
peak_times = result['peak_times']