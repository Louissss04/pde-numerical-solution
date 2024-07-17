function u = LeapFrog(a, delta_t, delta_x, x_start, x_end, t_start, t_end, fai, g)
    % LeapFrog - Solves the PDE u_t + a u_x = 0 using the Leap-Frog method
    %
    % Syntax:
    %   u = LeapFrog(a, delta_t, delta_x, x_start, x_end, t_start, t_end, fai, g)
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
    if fai(t_start) ~= g(x_start)
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

    % Startup step: Use Backward/Upwind method to calculate u(2,:)
    for i = 2:num_x_points-1
        u(2, i) = u(1, i) - nu * (u(1, i) - u(1, i-1));
    end
    u(2, 1) = fai(t(2));  % Apply boundary condition at x=0
    u(2, end) = u(1, end) - nu * (u(1, end) - u(1, end-1));  % Apply boundary condition at x=x_end

    % Time evolution: Use Leap-Frog method
    for n = 2:num_t_points-1
        for i = 2:num_x_points-1
            u(n+1, i) = u(n-1, i) - nu * (u(n, i+1) - u(n, i-1));
        end
        u(n+1, 1) = fai(t(n+1));  % Apply boundary condition at x=0
        u(n+1, end) = u(n, end) - nu * (u(n, end) - u(n, end-1));  % Apply BW method at x=x_end
    end
end
