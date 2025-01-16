close all
clear all

Gs = tf(5,[1,0.6,0.05]); %funçao continuo
T = 0.5/5; % tempo de amostragem


%funçao discreta
G = c2d(Gs,T,'zoh');

%% Controlador P

%plots
figure(1)
rlocus(G)
grid on;

figure(2)
step(G)
grid on;

CP = 0.0549;

%malha fechada controlador P
G_MF_P = feedback(CP*G,1);

figure(3)
step(CP*G)
grid on;

figure(4)
step(G_MF_P)
grid on;

stepinfo(G_MF_P)

%erros regime permanente
syms z;
% Função de transferência simbólica
Gz_sym = poly2sym(G.Numerator{1}, z) / poly2sym(G.Denominator{1}, z);

% 1. Entrada degrau: R(z) = z / (z - 1)
R_step = z / (z - 1);
step_G = R_step / (1 + CP*Gz_sym);
Ess_step = limit(step_G, z, 1);
disp(['Erro de regime permanente para degrau: ', char(Ess_step)]);

% 2. Entrada rampa: R(z) = z / (z - 1)^2
R_ramp = z / (z - 1)^2;
Ess_ramp = limit(R_ramp / (1 + Gz_sym), z, 1);
disp(['Erro de regime permanente para rampa: ', char(Ess_ramp)]);

% 3. Entrada parábola: R(z) = z / (z - 1)^3
R_parabola = z / (z - 1)^3;
Ess_parabola = limit(R_parabola / (1 + Gz_sym), z, 1);
disp(['Erro de regime permanente para parábola: ', char(Ess_parabola)]);


%% Controle PI--------------------------------
syms K Ki

%valor para a igualdade de K/KI
igualdade = 0.990049833749677;

%soluçao da equaçao pra achar K/KI
equacao = (K - Ki*0.1/2)/(K + Ki*0.1/2) == igualdade;
K_ratio = solve(equacao, K); % Primeiro resolver para K
K_Ki = double(K_ratio/Ki);

PI_inic = tf([1, -igualdade],[1,-1],T)

figure(5)
rlocus(PI_inic*G) %ganho 0.045
grid on;

K_linha_inic = 0.0455;

%% rodar ate aqui e colocar o valor de K_Ki no lugar de 10
%o ganho proporcional e proximo de 1 entao considera igual 1.
equacao1 = K + Ki*T/2 == 0.0455;

%valor de K_Ki.
equacao2 = K/Ki == 10; 

% Resolver o sistema de equações
solucao = solve([equacao1, equacao2], [K, Ki]);

%valores de K e Ki para o K'
Kp_linha = double(solucao.K);
Ki_linha = double(solucao.Ki);

%ganho proporcional K' que acha na MA do I*G
K_linha = K_linha_inic;
%valor do zero do controlador
K_zero = (Kp_linha - Ki_linha*T/2)/(Kp_linha + Ki_linha*T/2);

%controlador PI
PI_disc = K_linha*tf([1,-K_zero],[1,-1],T)

figure(5)
rlocus(PI_disc*G) %ganho 0.045
grid on;

figure(6)
step(PI_disc*G)
grid on;

%malha fechada
G_MF_PI = feedback(PI_disc*G,1);

figure(7)
step(G_MF_PI)
grid on;

stepinfo(G_MF_PI)

%residuos
num_prod = conv(PI_disc.Numerator{1}, G.Numerator{1});
den_prod = conv(PI_disc.Denominator{1}, G.Denominator{1});
[r,p,k] = residue(num_prod, den_prod);



