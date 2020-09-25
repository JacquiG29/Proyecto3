% =========================================================================
% MT3006 - LABORATORIO 4
% -------------------------------------------------------------------------
% Emplee el código base que se presenta a continuación para desarrollar un
% EKF que permita estimar el estado de un robot móvil diferencial equipado
% con encoders incrementales en sus ruedas y un GPS, según lo indicado en
% la guía del laboratorio. Modifique sólo en los lugares del
% código en donde se le indica.
% =========================================================================

%% Parámetros de la simulación
dt = 0.01; % período de muestreo
t0 = 0; % tiempo inicial
tf = 20; % tiempo final
N = (tf-t0)/dt; % número de iteraciones

%% Inicialización y condiciones iniciales
% Condiciones iniciales 
x0 = 0; % ***** PUEDE MODIFICAR *****
y0 = 0; % ***** PUEDE MODIFICAR *****
theta0 = pi/2; % ***** PUEDE MODIFICAR *****

xi0 = [x0; y0; theta0; 0; 0];
mu0 = [0; 0];
xi = xi0; % vector de estado 
mu = mu0; % vector de entradas

% Arrays para almacenar las trayectorias de las variables de estado,
% entradas y trayectoria del punto de interés en el plano de imagen
XI = zeros(numel(xi),N+1);
MU = zeros(numel(mu),N+1);
% Inicialización de arrays
XI(:,1) = xi0;
MU(:,1) = mu0;

% Arrays auxiliares para dashboard de sensores
GPS = zeros(2, N+1); 
ENC_R = zeros(1,N+1); 
ENC_L = zeros(1,N+1); 

%% Parámetros del Pioneer3-DX
robot_r = (195/1000)/2; % radio de las ruedas
robot_ell = 381/1000; % distancia entre ruedas
robot_N = 256; % ticks por revolución de los encoders

%% Definición de mapas, matrices y jacobianos para el EKF
% Modelo de odometría
Func = @(x,drho,dtheta) [x(1)+drho*cos(x(3)); x(2)+drho*sin(x(3)); x(3)+dtheta]; 

% Jacobianos del modelo de odometría
dFunc_dx = @(x,drho,dtheta) [1,0, -drho*sin(x(3)); 0,1, drho*cos(x(3)); 0,0,1];
dFunc_dw = @(x,drho,dtheta) [cos(x(3)),0; sin(x(3)),0; 0,1];

% Matriz del GPS
Cd = [1,0,0; 0,1,0];

%% Varianzas de proceso y medición estimadas
% Varianza del ruido de proceso (encoders)
varcase_sel = 3;
% Caso 1: sobre-estimando la varianza
if(varcase_sel == 1)
    sw1 = 0.1; % desplazamiento lineal, en metros 
    sw2 = 1*(pi/180); % desplazamiento angular, en radianes 
% Caso 2: varianza cercana a la real
elseif(varcase_sel == 2)
    sw1 = 0.0012; % desplazamiento lineal, en metros;
    sw2 = 0.0062; % desplazamiento angular, en radianes 
% Caso 3: valor intermedio
elseif(varcase_sel == 3)
    sw1 = 0.5*(0.1-0.0012) + 0.0012; % desplazamiento lineal, en metros;
    sw2 = 0.5*(1*(pi/180) - 0.0062) + 0.0062; % desplazamiento angular, en radianes 
end
    
% Varianza del ruido de observación (GPS)
sv1 = 0.4; % coordenada x, en metros
sv2 = 0.4; % coordenada y, en metros

%% Definición de variables para la odometría
% Medición previa del número de ticks por rueda
ticksR = 0; 
ticksL = 0;
% Longitud de los arcos de las ruedas y del centro
DL = 0;
DR = 0;
% Incrementos de posición lineal y angular
drho = 0;   
dtheta = 0;

% *************************************************************************
% DEFINICIONES E INICIALIZACIÓN PARA EL FILTRO DE KALMAN AQUÍ
% *************************************************************************
% Matrices de covarianza del ruido de proceso y observación
Qw = diag([sw1^2, sw2^2]);
Qv = diag([sv1^2, sv2^2]);

% Estimados a-PRIORi y a-POSTeriori del vector de estado
xhat0 = [x0, y0, theta0]';
sys_order = numel(xhat0);
xhat_prior = zeros(sys_order,1); 
xhat_post = xhat0;

% Inicialización de la matriz de covarianza del error de estimación
sigma_e = 0.01; % desviación pequeña porque conocemos con alta certeza x0
%Si aumento este valor,va a converger más rápido
%porque le estoy diciendo que estoy lejos de la estimación
%Hay que saber tunear este filtro de Kalman
P_prior = zeros(sys_order, sys_order);
P_post = sigma_e^2*eye(sys_order); 

% Array para almacenar la salida del filtro de Kalman
XHAT = zeros(numel(xhat_post),N+1);
XHAT(:,1) = xhat0; 
% *************************************************************************

%% Solución recursiva del sistema dinámico
for n = 0:N
    % *********************************************************************
    % COLOCAR EL CONTROLADOR AQUÍ
    % *********************************************************************
    % Velocidades de las ruedas del robot
    phiR = 5;
    phiL = 2;
    
    % Vector de entrada del sistema
    mu = [phiR; phiL];
    % *********************************************************************
    
    % Se adelanta un paso en la simulación del robot diferencial y se
    % obtienen las mediciones otorgadas por los sensores a bordo
    [gps_position, encoder_rticks, encoder_lticks, xi] = ...
        differential_drive(xi, mu, dt);
    
    % *********************************************************************
    % ESTIMAR EL ESTADO DEL ROBOT AQUÍ
    % *********************************************************************
    % Se calculan los incrementos de posición lineal y angular mediante
    % odometría
    DL = 2*pi*robot_r*((encoder_lticks - ticksL)/robot_N);
    DR = 2*pi*robot_r*((encoder_rticks - ticksR)/robot_N);
    ticksL = encoder_lticks; 
    ticksR = encoder_rticks;
    drho = (DR + DL)/2;
    dtheta = (DR - DL)/robot_ell;
    
    % Se generan los jacobianos a emplear en el EKF
    Ad = dFunc_dx(xhat_post, drho, dtheta);
    Fd = dFunc_dw(xhat_post, drho, dtheta);
    
    % Predicción
    xhat_prior = Func(xhat_post, drho, dtheta);
    P_prior = Ad*P_post*Ad' + Fd*Qw*Fd';
    
    % Corrección
    y = gps_position;
    L_k = P_prior*Cd'*(Qv + Cd*P_prior*Cd')^(-1);
    xhat_post = xhat_prior + L_k*(y - Cd*xhat_prior);
    P_post = P_prior - L_k*Cd*P_prior;
    
    % Se guarda la salida del filtro de Kalman
    XHAT(:,n+1) = xhat_post;
    % *********************************************************************
    
    % Se guardan las trayectorias del estado y las entradas
    XI(:,n+1) = xi;
    MU(:,n+1) = mu;
    
    % Se recolectan las mediciones de los sensores para desplegarlas en un
    % dashboard
    GPS(:,n+1) = gps_position;
    ENC_R(:,n+1) = encoder_rticks;
    ENC_L(:,n+1) = encoder_lticks;
end

% Trayectoria real del estado del robot
% ***** MODIFICAR PARA GRAFICAR ENCIMA EL ESTADO ESTIMADO *****
figure;
t = t0:dt:tf;
plot(t, XI(1:3,:)', 'LineWidth', 1);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, XHAT', '--', 'LineWidth', 1);
hold off;
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{x}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('$x(t)$', '$y(t)$', '$\theta(t)$', 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;

% Factor de escala para el "mundo" del robot
sf = 2; %***** PUEDE MODIFICAR *****

% Gráfica adicional que compara la salida del EKF con las mediciones del
% GPS
figure;
plot(t, GPS', 'LineWidth', 1);
hold on;
plot(t, XHAT(1:2,:)', 'LineWidth', 1);
hold off;
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathrm{pos}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
title('Comparación salida EKF vs GPS');
grid minor;

%% Animación y generación de figuras (NO modificar)
figure;
subplot(2,2,1);
plot(t, GPS', 'LineWidth', 1);
title('GPS');
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 14);
l = legend('$x(t)$', '$y(t)$', 'Location', 'best', 'Orientation', 'horizontal');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;
subplot(2,2,2);
stairs(t, ENC_R', 'LineWidth', 1);
title('Encoder rueda derecha');
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 14);
grid minor;
subplot(2,2,3);
stairs(t, ENC_L', 'LineWidth', 1);
title('Encoder rueda izquierda');
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 14);
grid minor;

figure;
xlim(sf*[-1, 1]+[-0.5, 0.5]);
ylim(sf*[-1, 1]+[-0.5, 0.5]);
grid minor;
hold on;

q = XI(:,1);
x = q(1); y = q(2); theta = q(3);

trajplot = plot(x, y, '--k', 'LineWidth', 1);

BV = [-0.1, 0, 0.1; 0, 0.3, 0];
IV = [cos(theta-pi/2), -sin(theta-pi/2); sin(theta-pi/2), cos(theta-pi/2)] * BV;
bodyplot = fill(IV(1,:) + x, IV(2,:) + y, [0.5,0.5,0.5]);

xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 16);
axis square;
hold off;

for n = 2:N+1
    q = XI(:,n);
    x = q(1); y = q(2); theta = q(3);
    
    trajplot.XData = [trajplot.XData, x];
    trajplot.YData = [trajplot.YData, y];
    
    BV = [-0.1, 0, 0.1; 0, 0.3, 0];
    IV = [cos(theta-pi/2), -sin(theta-pi/2); sin(theta-pi/2), cos(theta-pi/2)] * BV;
    bodyplot.XData = IV(1,:) + x;
    bodyplot.YData = IV(2,:) + y;
    
    pause(dt);
end