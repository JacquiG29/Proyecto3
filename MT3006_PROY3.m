% =========================================================================
% MT3006 - PROYECTO 3
% -------------------------------------------------------------------------
% Jacqueline Guarcax,16142
% Alejandro Windevoxhel,16047
% =========================================================================

%% Parámetros de la simulación
dt = 0.01; % período de muestreo
t0 = 0; % tiempo inicial
tf = 10; % tiempo final ***** PUEDE MODIFICAR *****
N = (tf-t0)/dt; % número de iteraciones

%% Generación de landmarks 
dist_obs=2;%forma en que se distribuiran obstaculos
num_obs=10;
obs=zeros(2,num_obs);

if dist_obs ==1
    obs = [-2.14,  2.90, -3.80, -2.13, 4.97, -0.38,  3.59, -1.39, 2.18, -1.97;
        -0.09, -2.64,  1.81,  2.99, 1.30, -4.56, -0.75,  0.08, 2.80,  4.46];
elseif dist_obs ==2    
    radio=4;
    for i=1:num_obs
        obs(:,i)=radio*[cos(2*pi*(i/num_obs));sin(2*pi*(i/num_obs))];
    end
end

%% Inicialización y condiciones iniciales
% Condiciones iniciales 
x0 = 4; % ***** PUEDE MODIFICAR *****
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

% Inicialización del mapa (array de landmarks encontrados)
map = zeros(2,10);

%% Generación de la trayectoria del robot
waypoints = [-3,-3; 2,-4; 5,-2; 4,1; 1,2; -2,4; -3,3; -3,1; -4,-1; -3,-3]';
tintrp = t0:dt:round(tf/2);
refx = interp1(linspace(t0,round(tf/2),10), waypoints(1,:), tintrp, 'spline');
refy = interp1(linspace(t0,round(tf/2),10), waypoints(2,:), tintrp, 'spline');
ref = [refx; refy];

%% Parámetros del robot y del controlador
robot_r = (195/1000)/2; % radio de las ruedas
robot_ell = 381/1000; % distancia entre ruedas
robot_N = 256; % ticks por revolución de los encoders

% Acercamiento exponencial
v0 = 4*2; % (ajustar según rendimiento del EKF)
alpha = 0.5;

% PID orientación
kpO = 4*2; % (ajustar según rendimiento del EKF)
kiO = 0.0001; 
kdO = 0;
EO = 0;
eO_1 = 0;

%% Funciones y jacobianos empleados en el mapeo
% Modelo de medición
G = @(xi, pi) [norm(pi-xi(1:2)); atan2(pi(2)-xi(2),pi(1)-xi(1)) - xi(3)];

% Función que estima la posición de los obstáculos observados con base en
% la pose actual del robot
g = @(xi, yi) [xi(1) + yi(1)*cos(xi(3)+yi(2));
             xi(2) + yi(1)*sin(xi(3)+yi(2))];

% Matriz de observación asociada al i-ésimo obstáculo encontrado
C_i = @(xi, pi) [(pi(1)-xi(1))/norm(pi-xi(1:2)), (pi(2)-xi(2))/norm(pi-xi(1:2));
                 -(pi(2)-xi(2))/norm(pi-xi(1:2))^2, (pi(1)-xi(1))/norm(pi-xi(1:2))^2];
             
C_eta = @(xi, pi) [-(pi(1)-xi(1))/norm(pi-xi(1:2)), -(pi(2)-xi(2))/norm(pi-xi(1:2)),0;
                   (pi(2)-xi(2))/norm(pi-xi(1:2))^2, -(pi(1)-xi(1))/norm(pi-xi(1:2))^2,-1];
% Jacobiano de g con respecto de las mediciones
G_eta = @(xi, yi) [1,0,-yi(1)*sin(xi(3)+yi(2));
                   0,1,yi(1)*cos(xi(3)+yi(2))];             
             
% Jacobiano de g con respecto de las mediciones
Gy = @(xi, yi) [cos(xi(3)+yi(2)), -yi(1)*sin(xi(3)+yi(2));
                sin(xi(3)+yi(2)), yi(1)*cos(xi(3)+yi(2))];

% Jacobiano de inserción para el aumento de la matriz de covarianza del
% error de estimación
ins_jacob = @(n, xi, yi) [eye(n), zeros(n,2);G_eta(xi,yi),zeros(2,n-3), Gy(xi,yi)]; 
 
% Modelo de odometría
Func = @(x,drho,dtheta) [x(1)+drho*cos(x(3)); x(2)+drho*sin(x(3)); x(3)+dtheta]; 

% Jacobianos del modelo de odometría
dFunc_dx = @(x,drho,dtheta) [1,0, -drho*sin(x(3)); 0,1, drho*cos(x(3)); 0,0,1];
dFunc_dw = @(x,drho,dtheta) [cos(x(3)),0; sin(x(3)),0; 0,1];

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
%% Parámetros e inicialización del filtro de Kalman
% Varianza del ruido de proceso (encoders)
sw1 = 0.0012; % desplazamiento lineal, en metros;
sw2 = 0.0062; % desplazamiento angular, en radianes
% Varianzas del ruido de observación (estimadas experimentalmente)
s_kappa = 0.0998;
s_beta = 0.0176;
% Varianzas del ruido de proceso (estimadas experimentalmente)
s_x = 0.0342;
s_y = 0.0996;
% Matrices de covarianza
Qw = diag([sw1^2, sw2^2]);
QvG = diag([s_kappa^2, s_beta^2]);
Qvg = diag([s_x^2, s_y^2]);
% Vector de estado a-priori y a-posteriori
x_hat_prior = [x0, y0, theta0]';
x_hat_post = [0,0,0]';
% Matriz de covarianza a-priori y a-posteriori
a=numel(x_hat_prior)
P_prior = [zeros(a,a)];
P_post = [zeros(a,a)];
% Número de obstáculos encontrados
M = 0;
% Tolerancia empleada para la correspondencia de obstáculos previamente 
% encontrados (ajustar según rendimiento del EKF)
tol = 4*sqrt(s_x^2 + s_y^2);
controlador=0;
%% Solución recursiva del sistema dinámico
for n = 0:N
    % Se adelanta un paso en la simulación del robot diferencial y se
    % obtienen las mediciones otorgadas por los sensores a bordo
    [gps_position, encoder_rticks, encoder_lticks, xi] = ...
        differential_drive(xi, mu, dt);
    [range, bearing] = distance_sensor(xi, obs);
    
    % Posición y orientación real del robot
    x = xi(1); y = xi(2); theta = xi(3);
    
    % Se calculan los incrementos de posición lineal y angular mediante
    % odometría
    DL = 2*pi*robot_r*((encoder_lticks - ticksL)/robot_N);
    DR = 2*pi*robot_r*((encoder_rticks - ticksR)/robot_N);
    ticksL = encoder_lticks;
    ticksR = encoder_rticks;
    %calcular odometria
    drho = (DR + DL)/2;
    dtheta = (DR - DL)/robot_ell;
    
    % *********************************************************************
    % COLOCAR EL CONTROLADOR AQUÍ
    % *********************************************************************
    % Controlador PID con acercamiento exponencial
    goal = ref(:, mod(n, length(ref))+1);
    e = goal - [x; y];
    thetag = atan2(e(2), e(1));

    eP = norm(e);
    eO = thetag - theta;
    eO = atan2(sin(eO), cos(eO));

    % Control de velocidad lineal
    kP = v0 * (1-exp(-alpha*eP^2)) / eP;
    v = kP*eP;

    % Control de velocidad angular
    eO_D = eO - eO_1;
    EO = EO + eO;
    w = kpO*eO + kiO*EO + kdO*eO_D;
    eO_1 = eO;

    % Velocidades de las ruedas del robot
    if controlador==1
        phiR = (2*v + 2*w*robot_ell) / (2*robot_r);
        phiL = (2*v - 2*w*robot_ell) / (2*robot_r);
    else
        phiR = 10;
        phiL = 8;
    end
    
    % Vector de entrada del sistema
    mu = [phiR; phiL];
    % *********************************************************************
    
    % *********************************************************************
    % IMPLEMENTAR EL EKF PARA LA CONSTRUCCIÓN DEL MAPA AQUÍ
    % *********************************************************************
    % Predicción
    % ---------------------------------------------------------------------
    x_Eta=Func(x_hat_prior(1:3), drho, dtheta);
    x_hat_prior(1:3)=x_Eta;
    x_hat_prior = x_hat_post;
    %P_prior = P_post;
    
    % Se generan los jacobianos a emplear en el EKF
    A = dFunc_dx(x_hat_prior(1:3), drho, dtheta);
    F = dFunc_dw(x_hat_prior(1:3), drho, dtheta);
    
    P_ee_prior=P_prior(1:3,1:3);
    P_ep_prior=P_prior(1:3,4:end);
    P_ep_prior_T=P_prior(4:end,1:3);
    P_pp_prior=P_prior(4:end,4:end);
    
    P_prior=[A*P_ee_prior*A'+F*Qw*F'   A*P_ep_prior;
             P_ep_prior_T*A'             P_pp_prior];
    % Corrección
    % ---------------------------------------------------------------------
    % Se pre-asigna el vector de innovaciones y la matriz de observación
    % con base en el número de obstáculos actuales
    z = zeros(2*M,1);
    C = [zeros(2*M,3),zeros(2*M,2*M)];
    
    % Se estima la posición de todos los obstáculos medidos con el sensor
    for i = 1:length(range)
        yi = [range(i); bearing(i)]; % medición
        pi_hat = g(xi,yi); % posición estimada
        obs_found = 0; % monitoreo de coincidencia de obstáculos previos
        % Se compara el obstáculo estimado con los estimados de todos los
        % obstáculos previamente encontrados, buscando hacer una
        % correspondencia
        for j = 1:M
            if((z(2*j-1) == 0) && (z(2*j) == 0))
                % Si hay una correspondencia se calcula la innovación y la
                % matriz de observación asociadas al obstáculo previo
                if(norm(pi_hat - x_hat_prior(2*j-1:2*j)) < tol)
                    z(2*j-1:2*j) = yi - G(xi, x_hat_post(2*j-1:2*j)); 
                    C(2*j-1:2*j,1:3)= C_eta(xi, x_hat_post(2*j-1:2*j));
                    C(2*j-1:2*j,2*j+2:2*j+3) = C_i(xi, x_hat_post(2*j-1:2*j));
                    obs_found = 1; % hay coincidencia
                end
            end
        end
        
        % En caso de no haber coincidencia con ninguno de los obstáculos
        % previos, se incluye el nuevo obstáculo y se aumenta el estado y
        % la matriz de covarianza del error
        if(~obs_found)
            x_hat_prior = [x_hat_prior; pi_hat];
            n_x=2*M+3;
            Yy = ins_jacob(n_x, xi, yi);
            P_prior = Yy*[P_prior, zeros(n_x,2); zeros(2,n_x), Qvg]*Yy';
            M = M + 1; % se incrementa el número de obstáculos encontrados
            % Se asegura que tanto la innovación del nuevo obstáculo como
            % su matriz de observación sean cero, para no tomarse en cuenta
            % dentro de la corrección sino hasta la siguiente observación
            % del obstáculo
            z = [z; 0;0];
            C_temp = C;
            C = [zeros(2*M,3),zeros(2*M,2*M)]; 
            C(1:2*(M-1), 1:3+2*(M-1)) = C_temp;
        end
    end
    
    % Se calcula la ganancia de Kalman y se actualizan los estimados
    % a-posteriori
    L = P_prior*C'*( C*P_prior*C' + kron(eye(M),QvG) )^(-1);
    x_hat_post = x_hat_prior + L*z;
    P_post = P_prior - L*C*P_prior;
    
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
%%
% Se construye el mapa con las posiciones estimadas de los obstáculos,
% estos se grafican como rombos verdes en la visualización de la
% simulación
map = reshape(x_hat_post(4:end), [2, numel(x_hat_post(4:end))/2]);

% Trayectoria real del estado del robot
figure;
t = t0:dt:tf;
plot(t, XI(1:3,:)', 'LineWidth', 1);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{x}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('$x(t)$', '$y(t)$', '$\theta(t)$', 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;

% Factor de escala para el "mundo" del robot
sf = 5; 


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

%plot(refx,refy,'--g');

q = XI(:,1);
x = q(1); y = q(2); theta = q(3);

obstacles = scatter(obs(1,:), obs(2,:), 'r', 'filled');
hold on;
scatter(map(1,:), map(2,:), 'g', 'd');
sensor_plot = gobjects(1,10);
for i = 1:length(obs)
    sensor_plot(i) = plot([xi(1),xi(1)], [xi(2),xi(2)], '--b', 'LineWidth', 1);
end
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
    
    for i = 1:length(obs)
        if(norm(obs(:,i) - q(1:2)) <= 2)
            sensor_plot(i).XData = [x, obs(1,i)];
            sensor_plot(i).YData = [y, obs(2,i)];
        else
            sensor_plot(i).XData = [x, x];
            sensor_plot(i).YData = [y, y];
        end
    end
    
    trajplot.XData = [trajplot.XData, x];
    trajplot.YData = [trajplot.YData, y];
    
    BV = [-0.1, 0, 0.1; 0, 0.3, 0];
    IV = [cos(theta-pi/2), -sin(theta-pi/2); sin(theta-pi/2), cos(theta-pi/2)] * BV;
    bodyplot.XData = IV(1,:) + x;
    bodyplot.YData = IV(2,:) + y;
    
    pause(dt);
end