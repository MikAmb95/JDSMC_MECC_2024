function  [x,iterCount,lastRes,JCost] = accelProjGradSolver(H,f,x0,xl,xu,MaxIter,xTol,printFlag)
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
iterCount = 1;

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

outerLoopRunCond = 1;
while outerLoopRunCond
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
    

    %% Check convergence criteria
    
    % Check break conditions
%     xConv = max(abs(x - xPrev));
    lastRes = norm(x - xPrev,Inf);
    xConv = lastRes;
        
    if iterCount >= MaxIter
        if printFlag
            fprintf('-------------------------------------------------- \n')
            fprintf(' *** Maximum iteration count of %0.0i exceeded. APG Algorithm Complete *** \n',MaxIter)
            fprintf(' *** Max update on completition is dx = %0.2e *** \n',xConv)
            fprintf('-------------------------------------------------- \n')
        end
        outerLoopRunCond = 0;
    elseif xConv < xTol
        if printFlag
            fprintf('-------------------------------------------------- \n')
            fprintf(' *** max(dx) < tolx = %0.2e. APG Algorithm Complete *** \n',xTol)
            fprintf(' *** Number of iteration for completition = %0.0i *** \n',iterCount)
            fprintf('-------------------------------------------------- \n')
        end
        outerLoopRunCond = 0;
    end
    
    % Track iteration number
    JCost = 1/2*x'*H*x + f'*x;  % REMEMBER THAT WE NEGLECT A TERM IN THE COST, THAT'S WHY ITS NEGATIVE

    if printFlag
        fprintf('* APG Iternation: %0.0i, Cost: %0.02e, Residual: %0.02e \n',iterCount,JCost,xConv)
    end

    iterCount = iterCount + 1;
    
   

end


end

