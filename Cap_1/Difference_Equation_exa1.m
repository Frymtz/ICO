% --- Simulação de uma G(z) de Segunda Ordem via Equação a Diferenças ---

clear;
clc;
close all;

% --- 1. Definição dos Coeficientes da Equação a Diferenças ---
% G(z) = (b0*z^-1 + b1*z^-2) / (1 + a1*z^-1 + a2*z^-2)
% y(k) = -a1*y(k-1) - a2*y(k-2) + b0*u(k-1) + b1*u(k-2)

% Coeficientes do denominador (associados a y):
% Para 1 + a1_den*z^-1 + a2_den*z^-2 = 0
a1_den = -1.5; % Coeficiente de y(k-1) na forma G(z)
a2_den =  0.7;  % Coeficiente de y(k-2) na forma G(z)

% Coeficientes do numerador (associados a u):
% Para b0_num*z^-1 + b1_num*z^-2
b0_num =  0.1;  % Coeficiente de u(k-1)
b1_num =  0.05; % Coeficiente de u(k-2)

% Coeficientes para a equação a diferenças na forma y(k) = ...
% y(k) = c_y1*y(k-1) + c_y2*y(k-2) + c_u1*u(k-1) + c_u2*u(k-2)
c_y1 = -a1_den; % = 1.5
c_y2 = -a2_den; % = -0.7
c_u1 = b0_num;  % = 0.1
c_u2 = b1_num;  % = 0.05

% --- 2. Definição do Sinal de Entrada u(k) ---
Nsim = 100; % Número de amostras para a simulação
k_vec = (0:Nsim-1)'; % Vetor de tempo discreto (índices de amostra)

% Sinal de entrada: Degrau unitário começando em k=5
u = zeros(Nsim, 1);
start_step_time = 5;
u(k_vec >= start_step_time) = 1;

% Outras opções de entrada (descomente uma para testar):
% u = ones(Nsim, 1); % Degrau unitário em k=0 (se k_vec começasse em 0)
% Impulso unitário em k=start_step_time
% u_impulse = zeros(Nsim, 1);
% u_impulse(k_vec == start_step_time) = 1;
% u = u_impulse;

% --- 3. Inicialização dos Vetores e Condições Iniciais ---
y = zeros(Nsim, 1); % Vetor para armazenar a saída do sistema

% Condições iniciais para os valores passados (assumindo sistema em repouso)
y_km1 = 0; % y(k-1) inicial
y_km2 = 0; % y(k-2) inicial
u_km1 = 0; % u(k-1) inicial
u_km2 = 0; % u(k-2) inicial

% --- 4. Implementação do Loop da Equação a Diferenças ---
for k = 1:Nsim % O índice k aqui corresponde ao instante de tempo k_vec(k)
    
    % Entrada atual u_k (para esta simulação, não é usada diretamente na
    % equação se não houver termo b_(-1)*u(k), mas útil para atualizar u_km1)
    u_k_atual = u(k); 
    
    % Calcular a saída y(k) usando a equação a diferenças
    % y(k) = c_y1*y(k-1) + c_y2*y(k-2) + c_u1*u(k-1) + c_u2*u(k-2)
    y_k_atual = c_y1*y_km1 + c_y2*y_km2 + c_u1*u_km1 + c_u2*u_km2;
    
    % Armazenar a saída calculada
    y(k) = y_k_atual;
    
    % Atualizar os valores passados para a próxima iteração
    % A ordem é importante: primeiro os mais antigos
    y_km2 = y_km1;
    y_km1 = y_k_atual;
    
    u_km2 = u_km1;
    u_km1 = u_k_atual; % u(k) atual se torna u(k-1) na próxima iteração
end

% --- 5. Verificação com a Função 'filter' do MATLAB (Opcional) ---
% Para a função filter:
% num_filter = [b_const, b_z_minus_1, b_z_minus_2, ...]
% den_filter = [a_const, a_z_minus_1, a_z_minus_2, ...]
% onde a_const é geralmente 1.
% A equação é a(1)y(k) + a(2)y(k-1) + ... = b(1)u(k) + b(2)u(k-1) + ...
% Nossa equação: y(k) - c_y1*y(k-1) - c_y2*y(k-2) = c_u1*u(k-1) + c_u2*u(k-2)
% y(k) + a1_den*y(k-1) + a2_den*y(k-2) = b0_num*u(k-1) + b1_num*u(k-2)
% Então:
den_filter = [1, a1_den, a2_den];    % [1, -1.5, 0.7]
num_filter = [0, b0_num, b1_num];    % [0,  0.1, 0.05] (o '0' é para o termo u(k))

y_filter = filter(num_filter, den_filter, u);

% --- 6. Plotar os Resultados ---
figure;
% Plot da entrada u(k)
subplot(2,1,1);
stem(k_vec, u, 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
title('Sinal de Entrada u(k)');
xlabel('Amostra (k)');
ylabel('Amplitude');
grid on;
ylim([-0.2, 1.2]);

% Plot da saída y(k) da simulação manual e da função filter
subplot(2,1,2);
stem(k_vec, y, 'r', 'Marker', 'o', 'DisplayName', 'Simulação Manual y(k)');
hold on;
plot(k_vec, y_filter, 'g--', 'LineWidth', 2, 'DisplayName', 'Saída de filter() y_{filter}(k)');
hold off;
title('Sinal de Saída y(k)');
xlabel('Amostra (k)');
ylabel('Amplitude');
legend show;
grid on;

sgtitle('Simulação de Sistema Discreto de 2ª Ordem');

% --- Análise Adicional (Opcional) ---
% Visualizar os pólos e zeros
% Gz = tf(num_filter, den_filter, -1); % -1 indica tempo de amostragem não especificado (discreto)
% disp('Função de Transferência G(z):');
% disp(Gz);
% pzmap(Gz);
% title('Diagrama de Pólos e Zeros de G(z)');