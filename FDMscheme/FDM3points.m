function u = FDM3points(a, delta_t, delta_x, x_end, t_end, scheme, fai, g)
    % FDM3points - Solves the PDE u_t + a u_x = 0 (a>0) using different
    %              3-points finite difference schemes
    % Syntax: u = FDM3points(a, delta_t, delta_x, x_end, t_end, scheme, fai, g)
    %
    % Inputs:
    %   a         - Advection speed
    %   delta_t   - Time step size
    %   delta_x   - Space step size
    %   x_end     - End point of the spatial domain
    %   t_end     - End point of the time domain
    %   scheme    - Scheme to use: 'LF' (Lax-Friedrichs), 'LW' (Lax-Wendroff), 
    %               'CD' (Central Difference), 'BW' (Backward/Upwind), 'FW' (Forward/Downwind)
    %   fai       - Function handle for boundary condition u(t,0)
    %   g         - Function handle for initial condition u(0,x)
    %
    % Outputs:
    %   u         - Solution matrix where each row is the state vector at a time step
    %
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

    % Determine Q based on the scheme
    switch scheme
        case 'LF'  % Lax-Friedrichs
            Q = 1;
        case 'LW'  % Lax-Wendroff
            Q = nu^2;
        case 'CD'  % Central Difference
            Q = 0;
        case 'BW'  % Backward/Upwind
            Q = nu;
        case 'FW'  % Forward/Downwind
            Q = -nu;
        otherwise
            error('Invalid scheme. Choose ''LF'', ''LW'', ''CD'', ''BW'', or ''FW''.');
    end

    % Initialize the spatial domain and initial condition
    x = linspace(0, x_end, num_x_points);
    u = zeros(num_t_points, num_x_points);
    u(1, :) = g(x);

    % Time evolution
    for n = 1:num_t_points-1
        u_next = zeros(size(u(n, :)));

        % Apply the finite difference scheme for internal points
        for i = 2:num_x_points-1
            u_next(i) = u(n, i) - (nu / 2) * (u(n, i+1) - u(n, i-1)) + ...
                        (Q / 2) * (u(n, i+1) - 2 * u(n, i) + u(n, i-1));
        end

        % Apply boundary condition at x=0
        u_next(1) = fai((n-1) * delta_t);

        % Apply upwind scheme for boundary at x=x_end
        u_next(end) = u(n, end) - (nu * (u(n, end) - u(n, end-1)));

        % Update the solution
        u(n+1, :) = u_next;
    end
end

