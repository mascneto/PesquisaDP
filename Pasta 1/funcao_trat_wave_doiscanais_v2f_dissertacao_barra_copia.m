clc
close all
clear all

cont = 1
tic

%%
%%%%% Bloco para renomear e numerar arquivos
for cont=1:50 % numero de pastas barra_1 a barra_50
Caminho = 'E:\Almir Carlos\barra\'
% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
Medicao = 'barra_6kv';
Medicao = strcat(Medicao,'_',num2str(cont),'\')
% Medicao = strcat(Medicao,'_',cont_med);
Arquivo = 'barra_6kv'; 
cont2 = cont+0;
Arq_cond = strcat(Arquivo,'_condi_',num2str(cont2))
Arq_charge = strcat(Arquivo,'_charge_',num2str(cont2))
Arquivo = strcat(Arquivo,'_',num2str(cont))
s = strcat(Caminho,Medicao,Arquivo,'_1')
load (s)


% function [absoluto] = tratamento_basico(s)
% pause(3)


format long

nn=0;
condi_sinal = [];
tensao_barra = [];
antena_antena = [];
charge_ldic = [];
antena_wave = [];
condi_wave = [];

atenuacao = 0.001; %se barra se outros 0.001

% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\Isolador\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Medicao = 'isolador_30kv_10\';
% Medicao = strcat(Medicao,'_',cont_med);
% Arquivo = 'isolador_30kv_10';
% s = strcat(Caminho,Medicao,Arquivo,'_1');
load(s)
% plot(Ch2)
% [maior,ind_maior] = max(Ch2);
% [menor,ind_menor] = min(Ch2);
% c_o=ind_menor-ind_maior;
% n_g = c_o/2;
% zero = ind_maior-n_g;
% inicio = zero+45000;
% tamanho_vetor = ind_menor+n_g+70000;
[inicio, tamanho_vetor] = corte_seno_mm(Ch2);
limiar = 30;

%    figure
%    plot(Ch1*1000) 
% %    title('Antena pura')
%    hold on
%    plot(Ch2*100)
%    
%     set(gca,'FontName','Times New Roman','FontSize',10)
%     yl = get(gca,'YTickLabel');
%     new_yl = strrep(yl(:),'.',',');
%     set(gca,'YTickLabel',new_yl)
%     xl = get(gca,'XTickLabel');
%     new_xl = strrep(xl(:),'.',',');
%     set(gca,'XTickLabel',new_xl)
%     xlabel('Pontos adquiridos(x10^6)')
%     ylabel('Amplitude (mV)')

for nn = 1:1:10
%     if nn>30 & nn<40
%         nn = nn+10;
%     end
    
%     load TP_26_DETETOR_20.mat
%     Arquivo=strcat(Caminho,Medicaok,'_',num2str(i));
%     load(Arquivo,Sensor(s,:),'dt','T')
    s = strcat(Caminho,Medicao,Arquivo,'_',num2str(nn));
%       s = strcat(Caminho,Medicao,'_',num2str(nn));
      nn

      
        load(s)
%         tamanho_vetor = length(Ch3);
        ampldetec = Ch3(inicio:tamanho_vetor,1);
        condi_sinal = [condi_sinal,ampldetec];
        
        antena_pura = Ch1(inicio:tamanho_vetor,1);
        antena_antena = [antena_antena, antena_pura];
        
        senoide_16 = Ch2(inicio:tamanho_vetor,1);
        tensao_barra = [tensao_barra, senoide_16];
        
        charge = Ch4(inicio:tamanho_vetor,1);
        charge_ldic = [charge_ldic, charge];
        
     antena_trat = wdenoise(antena_pura,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    antena_wave = [antena_trat,antena_wave];
    
     condi_trat = wdenoise(ampldetec,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    condi_wave = [condi_trat,condi_wave];
        
end
%
% plot(x,y,'Color',[1 0.5 0])
aten = 0.1; % atenuação
tam_plot = length(ampldetec);
% condi_sinal antena_antena charge_ldic
   figure
   plot(abs(condi_wave)*1000)
%    title('Antena pura')
   hold on
   plot(senoide_16*1000,'Color',[1 0.3 0])
   grid on
   xlim([0 tam_plot])
    set(gca,'FontName','Times New Roman','FontSize',10)
    yl = get(gca,'YTickLabel');
    new_yl = strrep(yl(:),'.',',');
    set(gca,'YTickLabel',new_yl)
    xl = get(gca,'XTickLabel');
    new_xl = strrep(xl(:),'.',',');
    set(gca,'XTickLabel',new_xl)
    xlabel('Pontos adquiridos (x10^6)')
    ylabel('Amplitude (mV)')

       figure
   plot(abs(antena_wave)*1000)
%    title('Antena pura')
   hold on
   plot(senoide_16*10,'Color',[1 0.3 0])
   grid on
   xlim([0 tam_plot])
    set(gca,'FontName','Times New Roman','FontSize',10)
    yl = get(gca,'YTickLabel');
    new_yl = strrep(yl(:),'.',',');
    set(gca,'YTickLabel',new_yl)
    xl = get(gca,'XTickLabel');
    new_xl = strrep(xl(:),'.',',');
    set(gca,'XTickLabel',new_xl)
    xlabel('Pontos adquiridos (x10^6)')
    ylabel('Amplitude (mV)')
    
           figure
   plot(charge_ldic*1000)
%    title('Antena pura')
   hold on
   plot(senoide_16*1000,'Color',[1 0.3 0])
   grid on
   xlim([0 tam_plot])
    set(gca,'FontName','Times New Roman','FontSize',10)
    yl = get(gca,'YTickLabel');
    new_yl = strrep(yl(:),'.',',');
    set(gca,'YTickLabel',new_yl)
    xl = get(gca,'XTickLabel');
    new_xl = strrep(xl(:),'.',',');
    set(gca,'XTickLabel',new_xl)
    xlabel('Pontos adquiridos (x10^6)')
    ylabel('Amplitude (mV)')

    
     figure
%    plot(charge_ldic*1000)
   plot(condi_sinal*1000)
%    title('Antena pura')
   hold on
   plot(senoide_16*1000)
   
    set(gca,'FontName','Times New Roman','FontSize',10)
    yl = get(gca,'YTickLabel');
    new_yl = strrep(yl(:),'.',',');
    set(gca,'YTickLabel',new_yl)
    xl = get(gca,'XTickLabel');
    new_xl = strrep(xl(:),'.',',');
    set(gca,'XTickLabel',new_xl)
    xlabel('Pontos adquiridos (x10^6)')
    ylabel('Amplitude (mV)')

%% Somente para determinaÃ§Ã£o de ruÃ­do

% maximo_antena = max(charge_ldic);
% minimo_antena = min(charge_ldic);
% maximo_maximo_antena = max(maximo_antena);
% minimo_minimo_antena = min(minimo_antena);
% media_ruido = mean(charge_ldic,1);

%%
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = condi_sinal; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada); # de onde vem essa derivada?
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
   
amplitudes = [];
fases = [];
    for j = 1:1:size(antena,2)
    [pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
    fases = [fases;locs];
    amplitudes = [amplitudes;pks];     

    end
    
novo = horzcat(amplitudes,fases);
   
   passo_plot = 360/(size(antena,1));
   eixo_x = [passo_plot:passo_plot:360];
   senoide_16 = senoide(:,1);
   eixo_360 = eixo_x';
%    
tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_cond)
save(nome_arquivo,'absoluto')



clearvars -except cont cont2 charge_ldic Arq_charge limiar senoide_16 tensao_barra

    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = charge_ldic; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
%    
amplitudes = [];
fases = [];
for j = 1:1:size(antena,2)
[pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
fases = [fases;locs];
amplitudes = [amplitudes;pks];     

end
novo = horzcat(amplitudes,fases);
 
passo_plot = 360/(size(antena,1));
eixo_x = [passo_plot:passo_plot:360];
senoide_16 = senoide(:,1);
eixo_360 = eixo_x';

% [linha,coluna,dimensao] = ind2sub(size(novo,1), find(fases==eixo_360))

tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_charge)
save(nome_arquivo,'absoluto')

clearvars -except cont
end
clear all
%%
for cont=1:50

Caminho = 'E:\Almir Carlos\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
Medicao = 'barra_6kv';
Medicao = strcat(Medicao,'_',num2str(cont),'\')
% Medicao = strcat(Medicao,'_',cont_med);
Arquivo = 'barra_6kv'; 
cont2 = cont+50;
Arq_cond = strcat(Arquivo,'_condi_',num2str(cont2))
Arq_charge = strcat(Arquivo,'_charge_',num2str(cont2))
Arquivo = strcat(Arquivo,'_',num2str(cont))
s = strcat(Caminho,Medicao,Arquivo,'_1')
load (s)


% function [absoluto] = tratamento_basico(s)
% pause(3)


format long

nn=0;
condi_sinal = [];
tensao_barra = [];
antena_antena = [];
charge_ldic = [];
antena_wave = [];
condi_wave = [];

atenuacao = 0.001; %se barra se outros 0.001

% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\Isolador\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Medicao = 'isolador_30kv_10\';
% Medicao = strcat(Medicao,'_',cont_med);
% Arquivo = 'isolador_30kv_10';
% s = strcat(Caminho,Medicao,Arquivo,'_1');
load(s)
% plot(Ch2)
% [maior,ind_maior] = max(Ch2);
% [menor,ind_menor] = min(Ch2);
% c_o=ind_menor-ind_maior;
% n_g = c_o/2;
% zero = ind_maior-n_g;
% inicio = zero;
% tamanho_vetor = ind_menor+n_g;
[inicio, tamanho_vetor] = corte_seno_mm(Ch2);
limiar = 30;
for nn = 11:1:20
    
%     load TP_26_DETETOR_20.mat
%     Arquivo=strcat(Caminho,Medicaok,'_',num2str(i));
%     load(Arquivo,Sensor(s,:),'dt','T')
    s = strcat(Caminho,Medicao,Arquivo,'_',num2str(nn));
%       s = strcat(Caminho,Medicao,'_',num2str(nn));
      nn

      
        load(s)
%         tamanho_vetor = length(Ch3);
        ampldetec = Ch3(inicio:tamanho_vetor,1);
        condi_sinal = [condi_sinal,ampldetec];
        
        antena_pura = Ch1(inicio:tamanho_vetor,1);
        antena_antena = [antena_antena, antena_pura];
        
        senoide_16 = Ch2(inicio:tamanho_vetor,1);
        tensao_barra = [tensao_barra, senoide_16];
        
        charge = Ch4(inicio:tamanho_vetor,1);
        charge_ldic = [charge_ldic, charge];
        
     antena_trat = wdenoise(antena_pura,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    antena_wave = [antena_trat,antena_wave];
    
     condi_trat = wdenoise(ampldetec,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    condi_wave = [condi_trat,condi_wave];
        
end
aten = 0.1;

%% Somente para determinaÃ§Ã£o de ruÃ­do

maximo_antena = max(condi_sinal);
minimo_antena = min(condi_sinal);
maximo_maximo_antena = max(maximo_antena);
minimo_minimo_antena = min(minimo_antena);
media_ruido = mean(condi_sinal,1);

%%
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = condi_sinal; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
   
amplitudes = [];
fases = [];
    for j = 1:1:size(antena,2)
    [pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
    fases = [fases;locs];
    amplitudes = [amplitudes;pks];     

    end
    
novo = horzcat(amplitudes,fases);
   
   passo_plot = 360/(size(antena,1));
   eixo_x = [passo_plot:passo_plot:360];
   senoide_16 = senoide(:,1);
   eixo_360 = eixo_x';
%    
tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_cond)
save(nome_arquivo,'absoluto')

% clearvars -except cont cont2 antena_wave Arq_charge limiar senoide_16 tensao_barra
clearvars -except cont cont2 charge_ldic Arq_charge limiar senoide_16 tensao_barra
    
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = charge_ldic; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
%    
amplitudes = [];
fases = [];
for j = 1:1:size(antena,2)
[pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
fases = [fases;locs];
amplitudes = [amplitudes;pks];     

end
novo = horzcat(amplitudes,fases);
 
passo_plot = 360/(size(antena,1));
eixo_x = [passo_plot:passo_plot:360];
senoide_16 = senoide(:,1);
eixo_360 = eixo_x';

% [linha,coluna,dimensao] = ind2sub(size(novo,1), find(fases==eixo_360))

tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_charge)
save(nome_arquivo,'absoluto')

clearvars -except cont
end

for cont=1:50

Caminho = 'E:\Almir Carlos\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
Medicao = 'barra_6kv';
Medicao = strcat(Medicao,'_',num2str(cont),'\')
% Medicao = strcat(Medicao,'_',cont_med);
Arquivo = 'barra_6kv'; 
cont2 = cont+100;
Arq_cond = strcat(Arquivo,'_condi_',num2str(cont2))
Arq_charge = strcat(Arquivo,'_charge_',num2str(cont2))
Arquivo = strcat(Arquivo,'_',num2str(cont))
s = strcat(Caminho,Medicao,Arquivo,'_1')
load (s)


% function [absoluto] = tratamento_basico(s)
% pause(3)


format long

nn=0;
condi_sinal = [];
tensao_barra = [];
antena_antena = [];
charge_ldic = [];
antena_wave = [];
condi_wave = [];

atenuacao = 0.001; %se barra se outros 0.001

% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\Isolador\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Medicao = 'isolador_30kv_10\';
% Medicao = strcat(Medicao,'_',cont_med);
% Arquivo = 'isolador_30kv_10';
% s = strcat(Caminho,Medicao,Arquivo,'_1');
load(s)
% plot(Ch2)
% [maior,ind_maior] = max(Ch2);
% [menor,ind_menor] = min(Ch2);
% c_o=ind_menor-ind_maior;
% n_g = c_o/2;
% zero = ind_maior-n_g;
% inicio = zero;
% tamanho_vetor = ind_menor+n_g;
[inicio, tamanho_vetor] = corte_seno_mm(Ch2);
limiar = 30;
for nn = 21:1:30
    
%     load TP_26_DETETOR_20.mat
%     Arquivo=strcat(Caminho,Medicaok,'_',num2str(i));
%     load(Arquivo,Sensor(s,:),'dt','T')
    s = strcat(Caminho,Medicao,Arquivo,'_',num2str(nn));
%       s = strcat(Caminho,Medicao,'_',num2str(nn));
      nn

      
        load(s)
%         tamanho_vetor = length(Ch3);
        ampldetec = Ch3(inicio:tamanho_vetor,1);
        condi_sinal = [condi_sinal,ampldetec];
        
        antena_pura = Ch1(inicio:tamanho_vetor,1);
        antena_antena = [antena_antena, antena_pura];
        
        senoide_16 = Ch2(inicio:tamanho_vetor,1);
        tensao_barra = [tensao_barra, senoide_16];
        
        charge = Ch4(inicio:tamanho_vetor,1);
        charge_ldic = [charge_ldic, charge];
        
     antena_trat = wdenoise(antena_pura,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    antena_wave = [antena_trat,antena_wave];
    
     condi_trat = wdenoise(ampldetec,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    condi_wave = [condi_trat,condi_wave];
        
end
aten = 0.1;

%% Somente para determinaÃ§Ã£o de ruÃ­do

maximo_antena = max(condi_sinal);
minimo_antena = min(condi_sinal);
maximo_maximo_antena = max(maximo_antena);
minimo_minimo_antena = min(minimo_antena);
media_ruido = mean(condi_sinal,1);

%%
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = condi_sinal; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
   
amplitudes = [];
fases = [];
    for j = 1:1:size(antena,2)
    [pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
    fases = [fases;locs];
    amplitudes = [amplitudes;pks];     

    end
    
novo = horzcat(amplitudes,fases);
   
   passo_plot = 360/(size(antena,1));
   eixo_x = [passo_plot:passo_plot:360];
   senoide_16 = senoide(:,1);
   eixo_360 = eixo_x';
%    
tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_cond)
save(nome_arquivo,'absoluto')



% clearvars -except cont cont2 antena_wave Arq_charge limiar senoide_16 tensao_barra
clearvars -except cont cont2 charge_ldic Arq_charge limiar senoide_16 tensao_barra
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = charge_ldic; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
%    
amplitudes = [];
fases = [];
for j = 1:1:size(antena,2)
[pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
fases = [fases;locs];
amplitudes = [amplitudes;pks];     

end
novo = horzcat(amplitudes,fases);
 
passo_plot = 360/(size(antena,1));
eixo_x = [passo_plot:passo_plot:360];
senoide_16 = senoide(:,1);
eixo_360 = eixo_x';

% [linha,coluna,dimensao] = ind2sub(size(novo,1), find(fases==eixo_360))

tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_charge)
save(nome_arquivo,'absoluto')

clearvars -except cont
end

for cont=1:50

Caminho = 'E:\Almir Carlos\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
Medicao = 'barra_6kv';
Medicao = strcat(Medicao,'_',num2str(cont),'\')
% Medicao = strcat(Medicao,'_',cont_med);
Arquivo = 'barra_6kv'; 
cont2 = cont+150;
Arq_cond = strcat(Arquivo,'_condi_',num2str(cont2))
Arq_charge = strcat(Arquivo,'_charge_',num2str(cont2))
Arquivo = strcat(Arquivo,'_',num2str(cont))
s = strcat(Caminho,Medicao,Arquivo,'_1')
load (s)


% function [absoluto] = tratamento_basico(s)
% pause(3)


format long

nn=0;
condi_sinal = [];
tensao_barra = [];
antena_antena = [];
charge_ldic = [];
antena_wave = [];
condi_wave = [];

atenuacao = 0.001; %se barra se outros 0.001

% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\Isolador\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Medicao = 'isolador_30kv_10\';
% Medicao = strcat(Medicao,'_',cont_med);
% Arquivo = 'isolador_30kv_10';
% s = strcat(Caminho,Medicao,Arquivo,'_1');
load(s)
% plot(Ch2)
% [maior,ind_maior] = max(Ch2);
% [menor,ind_menor] = min(Ch2);
% c_o=ind_menor-ind_maior;
% n_g = c_o/2;
% zero = ind_maior-n_g;
% inicio = zero;
% tamanho_vetor = ind_menor+n_g;
[inicio, tamanho_vetor] = corte_seno_mm(Ch2);
limiar = 30;
for nn = 31:1:40
    
%     load TP_26_DETETOR_20.mat
%     Arquivo=strcat(Caminho,Medicaok,'_',num2str(i));
%     load(Arquivo,Sensor(s,:),'dt','T')
    s = strcat(Caminho,Medicao,Arquivo,'_',num2str(nn));
%       s = strcat(Caminho,Medicao,'_',num2str(nn));
      nn

      
        load(s)
%         tamanho_vetor = length(Ch3);
        ampldetec = Ch3(inicio:tamanho_vetor,1);
        condi_sinal = [condi_sinal,ampldetec];
        
        antena_pura = Ch1(inicio:tamanho_vetor,1);
        antena_antena = [antena_antena, antena_pura];
        
        senoide_16 = Ch2(inicio:tamanho_vetor,1);
        tensao_barra = [tensao_barra, senoide_16];

        charge = Ch4(inicio:tamanho_vetor,1);
        charge_ldic = [charge_ldic, charge];
        
        
     antena_trat = wdenoise(antena_pura,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    antena_wave = [antena_trat,antena_wave];
    
     condi_trat = wdenoise(ampldetec,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    condi_wave = [condi_trat,condi_wave];
        
end
aten = 0.1;

%% Somente para determinaÃ§Ã£o de ruÃ­do

maximo_antena = max(condi_sinal);
minimo_antena = min(condi_sinal);
maximo_maximo_antena = max(maximo_antena);
minimo_minimo_antena = min(minimo_antena);
media_ruido = mean(condi_sinal,1);

%%
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = condi_sinal; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
   
amplitudes = [];
fases = [];
    for j = 1:1:size(antena,2)
    [pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
    fases = [fases;locs];
    amplitudes = [amplitudes;pks];     

    end
    
novo = horzcat(amplitudes,fases);
   
   passo_plot = 360/(size(antena,1));
   eixo_x = [passo_plot:passo_plot:360];
   senoide_16 = senoide(:,1);
   eixo_360 = eixo_x';
%    
tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_cond)
save(nome_arquivo,'absoluto')



% clearvars -except cont cont2 antena_wave Arq_charge limiar senoide_16 tensao_barra
clearvars -except cont cont2 charge_ldic Arq_charge limiar senoide_16 tensao_barra
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = charge_ldic; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
%    
amplitudes = [];
fases = [];
for j = 1:1:size(antena,2)
[pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
fases = [fases;locs];
amplitudes = [amplitudes;pks];     

end
novo = horzcat(amplitudes,fases);
 
passo_plot = 360/(size(antena,1));
eixo_x = [passo_plot:passo_plot:360];
senoide_16 = senoide(:,1);
eixo_360 = eixo_x';

% [linha,coluna,dimensao] = ind2sub(size(novo,1), find(fases==eixo_360))

tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_charge)
save(nome_arquivo,'absoluto')

clearvars -except cont
end
%%
for cont=1:50

Caminho = 'E:\Almir Carlos\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\barra\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
Medicao = 'barra_6kv';
Medicao = strcat(Medicao,'_',num2str(cont),'\')
% Medicao = strcat(Medicao,'_',cont_med);
Arquivo = 'barra_6kv'; 
cont2 = cont+200;
Arq_cond = strcat(Arquivo,'_condi_',num2str(cont2))
Arq_charge = strcat(Arquivo,'_charge_',num2str(cont2))
Arquivo = strcat(Arquivo,'_',num2str(cont))
s = strcat(Caminho,Medicao,Arquivo,'_1')
load (s)


% function [absoluto] = tratamento_basico(s)
% pause(3)


format long

nn=0;
condi_sinal = [];
tensao_barra = [];
antena_antena = [];
charge_ldic = [];
antena_wave = [];
condi_wave = [];

atenuacao = 0.001; %se barra se outros 0.001

% Caminho = 'C:\Users\almir\OneDrive - ee.ufcg.edu.br\Medicoes Mestrado\Isolador\'%'G:\Meu Drive\Mestrado UFCG\Semestre 4\Aquisições mestrado\Ensaios14032023\';
% Medicao = 'isolador_30kv_10\';
% Medicao = strcat(Medicao,'_',cont_med);
% Arquivo = 'isolador_30kv_10';
% s = strcat(Caminho,Medicao,Arquivo,'_1');
load(s)
% plot(Ch2)
% [maior,ind_maior] = max(Ch2);
% [menor,ind_menor] = min(Ch2);
% c_o=ind_menor-ind_maior;
% n_g = c_o/2;
% zero = ind_maior-n_g;
% inicio = zero;
% tamanho_vetor = ind_menor+n_g;
[inicio, tamanho_vetor] = corte_seno_mm(Ch2);
limiar = 30;
for nn = 41:1:50
    
%     load TP_26_DETETOR_20.mat
%     Arquivo=strcat(Caminho,Medicaok,'_',num2str(i));
%     load(Arquivo,Sensor(s,:),'dt','T')
    s = strcat(Caminho,Medicao,Arquivo,'_',num2str(nn));
%       s = strcat(Caminho,Medicao,'_',num2str(nn));
      nn

      
        load(s)
%         tamanho_vetor = length(Ch3);
        ampldetec = Ch3(inicio:tamanho_vetor,1);
        condi_sinal = [condi_sinal,ampldetec];
        
        antena_pura = Ch1(inicio:tamanho_vetor,1);
        antena_antena = [antena_antena, antena_pura];
        
        senoide_16 = Ch2(inicio:tamanho_vetor,1);
        tensao_barra = [tensao_barra, senoide_16];
        
        charge = Ch4(inicio:tamanho_vetor,1);
        charge_ldic = [charge_ldic, charge];
        
     antena_trat = wdenoise(antena_pura,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    antena_wave = [antena_trat,antena_wave];
    
     condi_trat = wdenoise(ampldetec,4, ...
    'Wavelet', 'dB2', ...
    'DenoisingMethod', 'UniversalThreshold', ...
    'ThresholdRule', 'Soft', ...
    'NoiseEstimate', 'LevelIndependent'); 
    condi_wave = [condi_trat,condi_wave];
        
end
aten = 0.1;

%% Somente para determinaÃ§Ã£o de ruÃ­do

maximo_antena = max(condi_sinal);
minimo_antena = min(condi_sinal);
maximo_maximo_antena = max(maximo_antena);
minimo_minimo_antena = min(minimo_antena);
media_ruido = mean(condi_sinal,1);

%%
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = condi_sinal; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
   
amplitudes = [];
fases = [];
    for j = 1:1:size(antena,2)
    [pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
    fases = [fases;locs];
    amplitudes = [amplitudes;pks];     

    end
    
novo = horzcat(amplitudes,fases);
   
   passo_plot = 360/(size(antena,1));
   eixo_x = [passo_plot:passo_plot:360];
   senoide_16 = senoide(:,1);
   eixo_360 = eixo_x';
%    
tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_cond)
save(nome_arquivo,'absoluto')



% clearvars -except cont cont2 antena_wave Arq_charge limiar senoide_16 tensao_barra
clearvars -except cont cont2 charge_ldic Arq_charge limiar senoide_16 tensao_barra
    tempo = [1:1:length(senoide_16)]';
    senoide = tensao_barra;
    antena = charge_ldic; %antena_barra;
    
    media = mean(antena);
    moda = mode(antena);
    
%%
%    quadrado = abs(derivada);
   quadrado = abs(antena);
%    quadrado = derivada.^2;
   maximo_teste = max(quadrado);
   maxi_maxi = max(maximo_teste);
   pulso29 = ones(size(antena));
%    figure
%    plot(quadrado)
%    hold on
%    plot(senoide_16*atenuacao)

for j =1:1:size(antena,2)    
   for i=1:1:length(tempo)
       if quadrado(i,j) > (limiar/100)*maxi_maxi%63 %2.56e-6
           pulso29(i,j) = quadrado(i,j);

       else
              pulso29(i,j)=0;
           i=i+1;
       end
   end
end
% 
    pulsoreal29 = ones(size(antena));
   for j =1:1:size(antena,2)  
      for i=1:1:length(tempo)
       if pulso29(i,j)==0
           pulsoreal29(i,j)=0;
       else
           pulsoreal29(i,j)=antena(i,j);
       end
      end
   end
%    
amplitudes = [];
fases = [];
for j = 1:1:size(antena,2)
[pks,locs] = findpeaks(abs(pulsoreal29(:,j)),tempo,'MinPeakDistance',3440);
fases = [fases;locs];
amplitudes = [amplitudes;pks];     

end
novo = horzcat(amplitudes,fases);
 
passo_plot = 360/(size(antena,1));
eixo_x = [passo_plot:passo_plot:360];
senoide_16 = senoide(:,1);
eixo_360 = eixo_x';

% [linha,coluna,dimensao] = ind2sub(size(novo,1), find(fases==eixo_360))

tam360 = size(eixo_360,1);
auxiliar_conversao = [];
for cont=1:1:size(novo,1)
    auxiliar_conversao(cont,1) = (novo(cont,2)/tam360)*360;    
end
novo_graus = horzcat(amplitudes,auxiliar_conversao);

absoluto = novo_graus;
if size(absoluto,1)>5000
    continue
end

%% 
absoluto = abs(novo_graus);
nome_arquivo = sprintf('%s.mat',Arq_charge)
save(nome_arquivo,'absoluto')

clearvars -except cont
end
toc