% =========================================================================
% MT3006 - LABORATORIO 5
% -------------------------------------------------------------------------
% Emplee el c�digo base que se presenta a continuaci�n para desarrollar un
% EKF que permita la construcci�n del mapa con los obst�culos/landmarks
% dados, seg�n lo indicado en la gu�a del laboratorio. Modifique s�lo en
% los lugares del c�digo en donde se le indica.
% =========================================================================

%% Par�metros de la simulaci�n
dt = 0.01; % per�odo de muestreo
t0 = 0; % tiempo inicial
tf = 20; % tiempo final ***** PUEDE MODIFICAR *****
N = (tf-t0)/dt; % n�mero de iteraciones

%% Generaci�n de landmarks 
obs = [-2.14,  2.90, -3.80, -2.13, 4.97, -0.38,  3.59, -1.39, 2.18, -1.97;
       -0.09, -2.64,  1.81,  2.99, 1.30, -4.56, -0.75,  0.08, 2.80,  4.46];

%% Inicializaci�n y condiciones iniciales
% Condiciones iniciales 
x0 = -2; % ***** PUEDE MODIFICAR *****
y0 = -2; % ***** PUEDE MODIFICAR *****
theta0 = pi/2; % ***** PUEDE MODIFICAR *****

xi0 = [x0; y0; theta0; 0; 0];
mu0 = [0; 0];
xi = xi0; % vector de estado 
mu = mu0; % vector de entradas

% Arrays para almacenar las trayectorias de las variables de estado,
% entradas y trayectoria del punto de inter�s en el plano de imagen
XI = zeros(numel(xi),N+1);
MU = zeros(numel(mu),N+1);
% Inicializaci�n de arrays
XI(:,1) = xi0;
MU(:,1) = mu0;

% Arrays auxiliares para dashboard de sensores
GPS = zeros(2, N+1); 
ENC_R = zeros(1,N+1); 
ENC_L = zeros(1,N+1); 

% Inicializaci�n del mapa (array de landmarks encontrados)
map = zeros(2,10);

%% Generaci�n de la trayectoria del robot
waypoints = [-3,-3; 2,-4; 5,-2; 4,1; 1,2; -2,4; -3,3; -3,1; -4,-1; -3,-3]';
tintrp = t0:dt:round(tf/2);
refx = interp1(linspace(t0,round(tf/2),10), waypoints(1,:), tintrp, 'spline');
refy = interp1(linspace(t0,round(tf/2),10), waypoints(2,:), tintrp, 'spline');
ref = [refx; refy];

%% Par�metros del robot y del controlador
robot_r = (195/1000)/2; % radio de las ruedas
robot_ell = 381/1000; % distancia entre ruedas
% Acercamiento exponencial
v0 = 4*2; % (ajustar seg�n rendimiento del EKF)
alpha = 0.5;

% PID orientaci�n
kpO = 4*2; % (ajustar seg�n rendimiento del EKF)
kiO = 0.0001; 
kdO = 0;
EO = 0;
eO_1 = 0;

%% Funciones y jacobianos empleados en el mapeo
% Modelo de medici�n
G = @(xi, pi) [norm(pi-xi(1:2)); atan2(pi(2)-xi(2),pi(1)-xi(1)) - xi(3)];

% Funci�n que estima la posici�n de los obst�culos observados con base en
% la pose actual del robot
g = @(xi, yi) [xi(1) + yi(1)*cos(xi(3)+yi(2));
             xi(2) + yi(1)*sin(xi(3)+yi(2))];

% Matriz de observaci�n asociada al i-�simo obst�culo encontrado
C_i = @(xi, pi) [(pi(1)-xi(1))/norm(pi-xi(1:2)), (pi(2)-xi(2))/norm(pi-xi(1:2));
                 -(pi(2)-xi(2))/norm(pi-xi(1:2))^2, (pi(1)-xi(1))/norm(pi-xi(1:2))^2];

% Jacobiano de g con respecto de las mediciones
Gy = @(xi, yi) [cos(xi(3)+yi(2)), -yi(1)*sin(xi(3)+yi(2));
                sin(xi(3)+yi(2)), yi(1)*cos(xi(3)+yi(2))];

% Jacobiano de inserci�n para el aumento de la matriz de covarianza del
% error de estimaci�n
ins_jacob = @(n, xi, yi) [eye(n), zeros(n,2); zeros(2,n), Gy(xi,yi)]; 
         
%% Par�metros e inicializaci�n del filtro de Kalman
% Varianzas del ruido de observaci�n (estimadas experimentalmente)
s_kappa = 0.0998;
s_beta = 0.0176;
% Varianzas del ruido de proceso (estimadas experimentalmente)
s_x = 0.0342;
s_y = 0.0996;
% Matrices de covarianza
QvG = diag([s_kappa^2, s_beta^2]);
Qvg = diag([s_x^2, s_y^2]);
% Vector de estado a-priori y a-posteriori
x_hat_prior = [];
x_hat_post = [];
% Matriz de covarianza a-priori y a-posteriori
P_prior = [];
P_post = [];
% N�mero de obst�culos encontrados
M = 0;
% Tolerancia empleada para la correspondencia de obst�culos previamente 
% encontrados (ajustar seg�n rendimiento del EKF)
tol = 4*sqrt(s_x^2 + s_y^2);

%% Soluci�n recursiva del sistema din�mico
for n = 0:N
    % Se adelanta un paso en la simulaci�n del robot diferencial y se
    % obtienen las mediciones otorgadas por los sensores a bordo
    [gps_position, encoder_rticks, encoder_lticks, xi] = ...
        differential_drive(xi, mu, dt);
    [range, bearing] = distance_sensor(xi, obs);
    
    % Posici�n y orientaci�n real del robot
    x = xi(1); y = xi(2); theta = xi(3);
    
    % *********************************************************************
    % COLOCAR EL CONTROLADOR AQU�
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
    phiR = (2*v + 2*w*robot_ell) / (2*robot_r);
    phiL = (2*v - 2*w*robot_ell) / (2*robot_r);
%     phiR = 6;
%     phiL = 4;
    
    % Vector de entrada del sistema
    mu = [phiR; phiL];
    % *********************************************************************
    
    % *********************************************************************
    % IMPLEMENTAR EL EKF PARA LA CONSTRUCCI�N DEL MAPA AQU�
    % *********************************************************************
    % Predicci�n
    % ---------------------------------------------------------------------
    x_hat_prior = x_hat_post;
    P_prior = P_post;
    
    % Correcci�n
    % ---------------------------------------------------------------------
    % Se pre-asigna el vector de innovaciones y la matriz de observaci�n
    % con base en el n�mero de obst�culos actuales
    z = zeros(2*M,1);
    C = zeros(2*M,2*M);
    
    % Se estima la posici�n de todos los obst�culos medidos con el sensor
    for i = 1:length(range)
        yi = [range(i); bearing(i)]; % medici�n
        pi_hat = g(xi,yi); % posici�n estimada
        obs_found = 0; % monitoreo de coincidencia de obst�culos previos
        % Se compara el obst�culo estimado con los estimados de todos los
        % obst�culos previamente encontrados, buscando hacer una
        % correspondencia
        for j = 1:M
            if((z(2*j-1) == 0) && (z(2*j) == 0))
                % Si hay una correspondencia se calcula la innovaci�n y la
                % matriz de observaci�n asociadas al obst�culo previo
                if(norm(pi_hat - x_hat_prior(2*j-1:2*j)) < tol)
                    z(2*j-1:2*j) = yi - G(xi, x_hat_post(2*j-1:2*j)); 
                    C(2*j-1:2*j,2*j-1:2*j) = C_i(xi, x_hat_post(2*j-1:2*j));
                    obs_found = 1; % hay coincidencia
                end
            end
        end
        
        % En caso de no haber coincidencia con ninguno de los obst�culos
        % previos, se incluye el nuevo obst�culo y se aumenta el estado y
        % la matriz de covarianza del error
        if(~obs_found)
            x_hat_prior = [x_hat_prior; pi_hat];
            Yy = ins_jacob(2*M, xi, yi);
            P_prior = Yy*[P_prior, zeros(2*M,2); zeros(2,2*M), Qvg]*Yy';
            M = M + 1; % se incrementa el n�mero de obst�culos encontrados
            % Se asegura que tanto la innovaci�n del nuevo obst�culo como
            % su matriz de observaci�n sean cero, para no tomarse en cuenta
            % dentro de la correcci�n sino hasta la siguiente observaci�n
            % del obst�culo
            z = [z; 0;0];
            C_temp = C;
            C = zeros(2*M, 2*M); 
            C(1:2*(M-1), 1:2*(M-1)) = C_temp;
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

% Se construye el mapa con las posiciones estimadas de los obst�culos,
% estos se grafican como rombos verdes en la visualizaci�n de la
% simulaci�n
map = reshape(x_hat_post, [2, numel(x_hat_post)/2]);

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


%% Animaci�n y generaci�n de figuras (NO modificar)
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