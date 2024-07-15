function u = LeapFrog(a, delta_t, delta_x, x_end, t_end, fai, g)
    % LeapFrog - Solves the PDE u_t + a u_x = 0 using the Leap-Frog method
    %
    % Syntax:
    %   u = LeapFrog(a, delta_t, delta_x, x_end, t_end, fai, g)
    %
    % Inputs:
    %   a         - Advection speed
    %   delta_t   - Time step size
    %   delta_x   - Space step size
    %   x_end     - End point of the spatial domain
    %   t_end     - End point of the time domain
    %   fai       - Function handle for boundary condition u(t,0)
    %   g         - Function handle for initial condition u(0,x)
    %
    % Outputs:
    %   u         - Solution matrix where each row is the state vector at a time step
    %
    %   Write by Qi Sun, July 2024.

    % Check if fai(0) equals g(0)
    if fai(0) ~= g(0)
        error('Boundary condition fai(0) must equal initial condition g(0).');
    end

    % Calculate number of steps
    num_t_points = floor(t_end / delta_t) + 1;
    num_x_points = floor(x_end / delta_x) + 1;

    % Calculate lambda and nu
    lambda = delta_t / delta_x;
    nu = a * lambda;

    % Initialize the spatial domain and initial condition
    x = linspace(0, x_end, num_x_points);
    u = zeros(num_t_points, num_x_points);
    u(1, :) = g(x);

    % Startup step: Use Backward/Upwind method to calculate u(2,:)
    for i = 2:num_x_points-1
        u(2, i) = u(1, i) - nu * (u(1, i) - u(1, i-1));
    end
    u(2, 1) = fai(delta_t);  % Apply boundary condition at x=0
    u(2, end) = u(1, end) - nu * (u(1, end) - u(1, end-1));  % Apply boundary condition at x=x_end

    % Time evolution: Use Leap-Frog method
    for n = 2:num_t_points-1
        for i = 2:num_x_points-1
            u(n+1, i) = u(n-1, i) - nu * (u(n, i+1) - u(n, i-1));
        end
        u(n+1, 1) = fai(n * delta_t);  % Apply boundary condition at x=0
        u(n+1, end) = u(n, end) - nu * (u(n, end) - u(n, end-1));  % Apply BW method at x=x_end
    end
end
