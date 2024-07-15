clear;
addpath('../FDMscheme/');

% Usage example
a = 1.0;                   % Advection speed
delta_t = 0.01;            % Time step size
delta_x = 0.01;            % Space step size
x_end = 1.0;               % End point of the spatial domain
t_end = 1.0;               % End point of the time domain
scheme = 'BW';             % Choose scheme ('LF', 'LW', 'CD', 'BW', 'FW')

% Boundary and initial conditions
fai = @(t) sin(2 * pi * t); % Boundary condition at x=0
g = @(x) 0;   % Initial condition at t=0

% Solve the PDE
u = FDM3points(a, delta_t, delta_x, x_end, t_end, scheme, fai, g)';

% Plot the results
x = linspace(0, x_end, floor(x_end / delta_x) + 1);
t = linspace(0, t_end, floor(t_end / delta_t) + 1);
[T, X] = meshgrid(t, x);

% Plot the surface
figure;
surf(X, T, u);
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Solution of the PDE using 3-points scheme');