%% Execute file for Kepler's Equation.

% Kepler's equation determines the relation of the time and angular
% displacement within an orbit. This execute file provides a Kepler
% Analsis of a given TLE.

% Prepare MATLAB Workspace.

clear
clc
format compact
format long

% Initialize computation time.

tic;

% TLE Data.

e = 0.0002714;                  % Eccentricity
M = 2.5;                        % Mean Anomaly, radians
n = 1.0027118*((2*pi)/86400);   % Mean Motion, rad/s

% Constants.

mu = 3.986004415e5;             % Gravitational parameter, km^3/s^2

if -pi < M && M < 0 || M > pi

    E_0 = M - e;

else

    E_0 = M - e;

end

% Find semimajor axis.

a = (mu/(n^2))^(1/3);           % Semimajor Axis, km

[E_n] = Keplers_Equation(E_0, e, M);

% Find True Anomaly.

nu = 2*atan(sqrt((1 + e)/(1 - e))*tan(E_n/2))*(180/pi);

if nu < 0
    
    nu = nu + 360;
    
else
    
    nu = nu;
    
end

E_star = E_n;

[count, epsilon_t] = Kepler_Equation_Error(E_star, E_0, e, M);

%% Plot Newton-Raphson Error.

figure(1)
plot(count(1,:), epsilon_t(1,:)), axis ([0 5 0 6e-12]), grid on,
title('Newton-Raphson Error'), xlabel('Iteration'), ylabel('Error')

%% Mean Anomaly vs. True Anomaly Plot.

E = -pi:pi/60:pi;

M = @(E) E + e*sin(E);
nu = 2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))*(180/pi);

E_deg = E*(180/pi);
M_deg = (E + e*sin(E))*(180/pi);
figure(2)
plot(E_deg,M_deg), grid on, hold on,
title('TJS - 9, Mean Anomaly vs. True Anomaly') 
legend('e = 0.0004636'),xlabel('True Anomaly, \nu (degrees,^o)','FontWeight','bold'),
ylabel('Mean Anomaly, M (degrees, ^o)','FontWeight','bold'),
axis([-180 180 -180 180]), 
xticks([-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180]),
xticklabels({'-180^o','-150^o','-120^o','-90^o','-60^o','-30^o','0^o',...
    '30^o','60^o','90^o','120^o','150^o','180^o'})
yticks([-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180]),
yticklabels({'-180^o','-150^o','-120^o','-90^o','-60^o','-30^o','0^o',...
    '30^o','60^o','90^o','120^o','150^o','180^o'})

%% Mean Anomaly vs. Eccentric Anomaly Plot.

M = @(E) E + e*sin(E);

E_deg = E*(180/pi);
M_deg = (E + e*sin(E))*(180/pi);

figure(3)
plot(E_deg,M_deg), grid on, 
title('TJS - 9, Mean Anomaly vs. Eccentric Anomaly') 
legend('e = 0.0004636'),xlabel('Eccentric Anomaly, E (degrees,^o)','FontWeight','bold'),
ylabel('Mean Anomaly, M (degrees, ^o)','FontWeight','bold'), hold off

t_nit = toc;
fprintf('Calculation time = %4.6f [s] \n',t_nit);
