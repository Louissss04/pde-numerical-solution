function u = NonlinearHyperbolicSolver(f, delta_t, delta_x, x_start, x_end, t_start, t_end, scheme, phi, g)
    % NonlinearHyperbolicSolver - Solves the PDE u_t + f(u)_x = 0 using different numerical schemes
    %
    % Syntax:
    %   u = NonlinearHyperbolicSolver(f, delta_t, delta_x, x_start, x_end, t_start, t_end, scheme, phi, g)
    %
    % Inputs:
    %   f         - Function handle for the nonlinear function f(u)
    %   delta_t   - Time step size
    %   delta_x   - Space step size
    %   x_start   - Start point of the spatial domain
    %   x_end     - End point of the spatial domain
    %   t_start   - Start point of the time domain
    %   t_end     - End point of the time domain
    %   scheme    - Scheme to use: 'BW' (Backward/Upwind), 'CD' (Central Difference), 
    %               'LF' (Lax-Friedrichs), 'LW' (Lax-Wendroff), 'LWorigin' (Lax-Wendroff Origin), 
    %               'CIR' (Courant-Isaacson-Rees), 'MCIR' (Murman-Courant-Isaacson-Rees)
    %   phi       - Function handle for boundary condition u(t,0)
    %   g         - Function handle for initial condition u(0,x)
    %
    % Outputs:
    %   u         - Solution matrix where each row is the state vector at a time step
    %
    %   Written by Qi Sun, July 2024.

    % Check if phi(t_start) equals g(x_start)
    if phi(t_start) - g(x_start) >= 1e-8
        error('Boundary condition phi(t_start) must equal initial condition g(x_start).');
    end

    % Calculate number of steps
    num_t_points = floor((t_end - t_start) / delta_t) + 1;
    num_x_points = floor((x_end - x_start) / delta_x) + 1;

    % Calculate lambda
    lambda = delta_t / delta_x;

    % Initialize the spatial domain and initial condition
    x = linspace(x_start, x_end, num_x_points);
    t = linspace(t_start, t_end, num_t_points);
    u = zeros(num_t_points, num_x_points);
    u(1, :) = g(x);

    % Get f_prime
    syms U
    f_sym = f(U);
    f_prime_sym = diff(f_sym, U);

    % Check if f_prime is a constant
    if isAlways(diff(f_prime_sym, U) == 0)
        f_prime_constant = double(f_prime_sym);
        f_prime = @(u) f_prime_constant;
    else
        f_prime = matlabFunction(f_prime_sym);
    end

    % Define the numerical flux function based on the scheme
    switch scheme
        case 'BW'  % Backward/Upwind
            flux = @(u, v) (0.5 * sign(double(f_prime(u))) + 0.5) .* (f(u) - f(v)) + f(v);
        case 'CD'  % Central Difference
            flux = @(u, v) 0.5 * (f(u) + f(v));
        case 'LF'  % Lax-Friedrichs
            flux = @(u, v) 0.5 * (f(u) + f(v)) + 0.5 * lambda * (u - v);
        case 'LW'  % Lax-Wendroff
            flux = @(u, v) 0.5 * (f(u) + f(v)) - 0.5 * lambda * double(f_prime((u + v) / 2)) .* (f(v) - f(u));
        case 'LWorigin'  % Lax-Wendroff Origin
            flux = @(u, v) 0.5 * (f(u) + f(v)) - 0.5 * lambda * (double(f_prime(u)).^2 + double(f_prime(v)).^2) / 2 .* (v - u);
        case 'CIR'  % Courant-Isaacson-Rees
            flux = @(u, v) 0.5 * (f(u) + f(v)) - 0.5 * abs(double(f_prime((u + v) / 2))) .* (v - u);
        case 'MCIR'  % Murman-Courant-Isaacson-Rees
            flux = @(u, v) 0.5 * (f(u) + f(v)) - 0.5 * abs(f(u) - f(v)) .* sign(v - u);
        otherwise
            error('Invalid scheme. Choose ''BW'', ''CD'', ''LF'', ''LW'', ''LWorigin'', ''CIR'', or ''MCIR''.');
    end

    % Time evolution
    for n = 1:num_t_points-1
        u_next = zeros(size(u(n, :)));

        % Apply the numerical scheme for internal points
        for j = 2:num_x_points-1
            F_j_plus_half = flux(u(n, j), u(n, j+1));
            F_j_minus_half = flux(u(n, j-1), u(n, j));
            u_next(j) = u(n, j) - lambda * (F_j_plus_half - F_j_minus_half);
        end

        % Apply boundary condition at x=0
        u_next(1) = phi(t(n+1));

        % Apply the numerical scheme for the boundary at x=x_end
        F_j_minus_half = flux(u(n, num_x_points-1), u(n, num_x_points));
        u_next(num_x_points) = u(n, num_x_points) - lambda * (flux(u(n, num_x_points), u(n, num_x_points)) - F_j_minus_half);

        % Update the solution
        u(n+1, :) = u_next;
    end
end
