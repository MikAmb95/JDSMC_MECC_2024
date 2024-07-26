function [X,U] = runShrinkingHorizonMPC_generalDist(x0,U_MPC,controlArgs)
% Calculates all of the shirnking horizon constants and stores the relevant
% info in consts. Simulatenously runs the simulation 


global Ad Bd 
n  = size(Ad,1);
m  = size(Bd,2);
    
% Get MPC parameters
Q = controlArgs.Q;
R = controlArgs.R;
P = controlArgs.P;
K = controlArgs.K;
N = controlArgs.N;
d = controlArgs.dMax;

% Create state and input matrices for outputs
X = zeros(n,N+1);
U = zeros(m,N);
X(:,1) = x0;
Ui  = controlArgs.U_Init;

% Initialize the online estimated parameters
VStar = zeros(N,1);

% -------------- Calculate SHMPC constants for i = 1 step --------------
% Run loop and find constants
for i=1:N
    x_k = X(:,i);
    
    % Apply solution and perturb
    u_k = U_MPC(i);
    %X(:,i+1) = Ad*x_k + Bd*u_k + sign(x_k(1))*[d;0];
    X(:,i+1) = Ad*x_k + Bd*u_k;
    
    % Store
    U(:,i) = u_k;
    
    
    
end


end

