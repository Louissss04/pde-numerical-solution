function u_next = FDM3points(u, a, delta_t, delta_x, scheme)
    % solve_pde - Solves the PDE u_t + a u_x = 0 using different finite difference schemes
    %
    % Syntax: u_next = solve_pde(u, a, delta_t, delta_x, scheme)
    %
    % Inputs:
    %   u        - Current state vector
    %   a        - Advection speed
    %   delta_t  - Time step size
    %   delta_x  - Space step size
    %   scheme   - Scheme to use: 'LF' (Lax-Friedrichs), 'LW' (Lax-Wendroff), 
    %              'CD' (Central Difference), 'BW' (Backward/Upwind), 'FW' (Forward/Downwind)
    %
    % Outputs:
    %   u_next   - Next state vector
    
    % Important constants
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
    
    % Get the number of spatial points
    n = length(u);
    
    % Initialize the next state vector
    u_next = zeros(size(u));
    
    % Apply the finite difference scheme
    for i = 2:n-1
        u_next(i) = u(i) - (nu / 2) * (u(i+1) - u(i-1)) + ...
                    (Q / 2) * (u(i+1) - 2 * u(i) + u(i-1));
    end
    
    % Handle boundary conditions (assuming periodic boundaries for simplicity)
    u_next(1) = u(1) - (nu / 2) * (u(2) - u(n)) + ...
                (Q / 2) * (u(2) - 2 * u(1) + u(n));
    u_next(n) = u(n) - (nu / 2) * (u(1) - u(n-1)) + ...
                (Q / 2) * (u(1) - 2 * u(n) + u(n-1));
end
