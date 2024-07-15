clear;
addpath('../FDMscheme/');

% Advection speed
a = 1.0;

% Time step size
delta_t = 0.01;

% Space step size
delta_x = 0.01;

% End point of the spatial domain
x_end = 1.0;

% End point of the time domain
t_end = 1.0;

% Boundary and initial conditions
fai = @(t) sin(2 * pi * t); % Boundary condition at x=0
g = @(x) 0;   % Initial condition at t=0

% Solve the PDE
u = LeapFrog(a, delta_t, delta_x, x_end, t_end, fai, g);

% Plot the result
x = linspace(0, x_end, floor(x_end / delta_x) + 1);
t = linspace(0, t_end, floor(t_end / delta_t) + 1);
[T, X] = meshgrid(t, x);

% Transpose u to align with T and X
u = u';

% Plot the surface
figure;
surf(X, T, u);
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Solution of the advection equation using the Leap-Frog method');
