function u = Boxscheme(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
    % Boxscheme - Solves the PDE u_t + a u_x = 0 using the Box scheme method
    %
    % Syntax:
    %   u = Boxscheme(a, delta_t, delta_x, x_start, x_end, t_start, t_end, phi, g)
    %
    % Inputs:
    %   a         - Advection speed
    %   delta_t   - Time step size
    %   delta_x   - Space step size
    %   x_start   - Start point of the spatial domain
    %   x_end     - End point of the spatial domain
    %   t_start   - Start point of the time domain
    %   t_end     - End point of the time domain
    %   fai       - Function handle for boundary condition u(t,0)
    %   g         - Function handle for initial condition u(0,x)
    %
    % Outputs:
    %   u         - Solution matrix where each row is the state vector at a time step
    %
    %   Written by Qi Sun, July 2024.

    % Check if fai(t_start) equals g(x_start)
    if phi(t_start) - g(x_start) >= 1e-8
        error('Boundary condition fai(t_start) must equal initial condition g(x_start).');
    end

    % Calculate number of steps
    num_t_points = floor((t_end - t_start) / delta_t) + 1;
    num_x_points = floor((x_end - x_start) / delta_x) + 1;

    % Calculate lambda and nu
    lambda = delta_t / delta_x;
    nu = a * lambda;

    % Initialize the spatial domain and initial condition
    x = linspace(x_start, x_end, num_x_points);
    t = linspace(t_start, t_end, num_t_points);
    u = zeros(num_t_points, num_x_points);
    u(1, :) = g(x);

    % Time evolution using Box scheme
    for n = 1:num_t_points-1
        % Temporary array to store intermediate step values
        u_next = zeros(1, num_x_points);

        % Apply boundary condition at x=0
        u_next(1) = phi(t(n+1));

        % Apply the Box scheme for internal points
        for j = 2:num_x_points
            u_next(j) = (u(n, j-1) + ((1 - nu)/(1 + nu)) * (u(n, j) - u_next(j-1)));
        end

        % Update the solution
        u(n+1, :) = u_next;
    end
end
