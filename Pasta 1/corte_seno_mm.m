function [in, fim] = corte_seno_mm (Ch_seno)
%Função que determina os zeros da senoide adquirida...
...passando um filtro de média móvel.
    ...este filtro já existe nativamente no matlab
    ... a janela eu ajustei empiricamente e coloquei ajustes
    ...no início e fim dos vetores.
    
windowSize = 5000; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,Ch_seno);
% figure
% plot(t,Ch_seno)
% hold on
% plot(t,y)
% legend('Input Data','Filtered Data')
% figure
% plot(y)
indices_passagem_zero_janela= find(diff(sign(y)) ~= 0);
% disp(indices_passagem_zero);
fim = indices_passagem_zero_janela(3,1)+2000;
in = indices_passagem_zero_janela(2,1) -(indices_passagem_zero_janela(3,1)-indices_passagem_zero_janela(2,1))-2000;
end

