function [X,U,iter] = runShrinkingHorizonMPC_fixedIter(x0,controlArgs,output)
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
MaxIter = controlArgs.MaxIter;
umax = controlArgs.umax;
umin = controlArgs.umin;
d = controlArgs.dMax;
lambda_R = max(eig(R));
gamma = controlArgs.gamma;
alpha = controlArgs.alpha;
lambda = controlArgs.lambda;
betaVec = controlArgs.betaVec;

% Create state and input matrices for outputs
X = zeros(n,N+1);
U = zeros(m,N);
X(:,1) = x0;
Ui  = controlArgs.U_Init*0;

% Initialize the online estimated parameters
VStar = zeros(N,1);

iter = [];

% -------------- Calculate SHMPC constants for i = 1 step --------------
% Run loop and find constants
for i=1:N
        N_k = N - i + 1;
    x_k = X(:,i);
    
    % Get n*m length vectors u_low and u_up
    u_low = zeros(m*N_k,1);
    u_up = zeros(m*N_k,1);
    for j = 1:N
        u_low(1 + m*(j-1) : m + m*(j-1)) = umin;
        u_up(1 + m*(j-1) : m + m*(j-1)) = umax;
    end

    % Generate the constant QP matrices (fixed horizon lengths)
    H = output.matrices(i).H;
    G = output.matrices(i).G;
    W = output.matrices(i).W;
    
    % Calculate true value function and bound
    H_QP = H;
    f_QP = G*x_k;
    
    % Calculate shifted solution
    if i > 1
        uBar_i = Ui(m+1:end);
    else
        uBar_i = Ui;
    end

        
    % Calculate approx. MPC solution
    [Ui,iter(i)] = accelProjGradSolver(H_QP,f_QP,uBar_i,u_low,u_up,output.ellBarVec(i),1e-6,0);
    u_k = Ui(1:m);

    % Apply solution
    X(:,i+1) = Ad*x_k + Bd*u_k + sign(x_k(1))*[d;0];
    %X(:,i+1) = Ad*x_k + Bd*u_k;
    % Store
    U(:,i) = u_k;
end


end

