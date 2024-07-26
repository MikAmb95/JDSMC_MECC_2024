function  [x,solError,iterCount] = accelProjGradSolver_SHMPC(H,f,x0,xl,xu,MaxIter,eHat)
% Solves the quadratic programming problem:
% min 0.5*x'*H*x + f'*x   subject to:  xl <= x <= xu
% using the projected gradient method...

% * May eventually be updated to solve for general linear constraints, but
% current this code solves for box constraints only

% Check to make sure provided initial condition matches the size of the QP
% matrices
if length(x0) ~= size(H,2)
    error('Mismatch in the size of x0. Check provided x0')
end


x = x0; % Initialize x
n = length(x);
iterCount = 0;

% Find the Lipshitz constant L and set step-size
eigVec = eig(H);
L = max(eigVec);
m = min(eigVec);
kappa = L/m;
t = 1/L;

% Initialize Nesterov quantity
theta_prev = 0;
theta = 1;
z = x*0;

% Initialize cost gap
eigFactor = 1/sqrt(m);

% Calculate unconstrained optimizer and cost
COST = 0.5*x'*H*x + f'*x;
x_unc = -H\f;
J_unc = 0.5*x_unc'*H*x_unc + f'*x_unc;
solError = 2*eHat; % initialize this way so we at least do 1 iteration
% solError = eigFactor*(COST - J_unc);

while iterCount < MaxIter && solError > eHat
    % Store previous x value to check the convergence tolerance
    xPrev = x;
    
    % Compute intermediate Nesterov quantities
    gamma = theta_prev^2*L;
    y = x + (theta*gamma)/(gamma + m*theta) * (z - x);
    
    % Current gradient direction and step
    g = H*x + f; 
    xStar = y - t*g;
    
    % Project any components exceeding the constraints back onto the box
    % constraint
    for i = 1:n
        if xStar(i) > xu(i)
            x(i) = xu(i);
        elseif xStar(i) < xl(i)
            x(i) = xl(i);
        else
            x(i) = xStar(i);
        end
    end
    
    % Compute z and theta for the next iteration
    z = xPrev + (1/theta)*(x - xPrev);
    zeta = theta^2 - 1/kappa;
    theta_prev = theta;
    theta = (zeta/2)*( sqrt(1 + 4*theta^2/zeta^2) - 1 );
    
    % Calculate cost and solution gap
    COST = 0.5*x'*H*x + f'*x;
    solError = eigFactor*(COST - J_unc);
    
    % Track iteration number
    iterCount = iterCount + 1;
end


end

