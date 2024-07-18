clear;
addpath('../FDMscheme/');

% Advection speed
a = 1.0;

% Time step size
delta_t = 0.01;

% Space step size
delta_x = 0.01;

% Spatial domain
x_start = 0.0;
x_end = 1.0;

% Time domain
t_start = 0.0;
t_end = 1.0;

% Boundary and initial conditions
fai = @(t) sin(2 * pi * t); % Boundary condition at x=0
g = @(x) 0;   % Initial condition at t=0

% Solve the PDE
u = ImplicitCentralDifference(a, delta_t, delta_x, x_start, x_end, t_start, t_end, fai, g);

% Plot the result
x = linspace(x_start, x_end, floor((x_end - x_start) / delta_x) + 1);
t = linspace(t_start, t_end, floor((t_end - t_start) / delta_t) + 1);
[T, X] = meshgrid(t, x);

% Transpose u to align with T and X
u = u';

% Plot the surface
figure;
surf(X, T, u);
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Solution of the advection equation using the implicit central difference method');
