function [output] = calcParamsSHMPC(controlArgs)
% Calculates all of the shirnking horizon constants and stores the relevant
% info in consts. Simulatenously runs the simulation 


global Ad Bd 
n  = size(Ad,1);
m  = size(Bd,2);
    
% Get MPC parameters
Q = controlArgs.Q;
R = controlArgs.R;
K = dlqr(Ad,Bd,Q,R);
P = controlArgs.P;
N = controlArgs.N;
umax = controlArgs.umax;
umin = controlArgs.umin;
gamma = controlArgs.gamma;
alpha = controlArgs.alpha;
lambda = controlArgs.lambda;
betaVec = controlArgs.betaVec;
dMax = controlArgs.dMax;
X0 = controlArgs.X0;

% Calculate some stuff that we'll use later
eigP = eig(P);
lambda_p = lambda - min(eig(sqrtm(inv(P))*Q*sqrtm(inv(P))));

% -------------- Calculate SHMPC constants for i = 1 step --------------
% Generate the first set of maatrices
i = 1;
k = i - 1;
N_k = N - k;
xmax = ones(n,1);
xmin = -ones(n,1); % not actually used, just need to pass something
[H,G,W,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_k,Ad,Bd,lambda*P,Q,R,xmax,xmin,umax,umin);
H = 0.5*(H+H');
% Store matrices and calculate some things we need
eigH = eig(H);
initStruct.H = H;
initStruct.G = G;
initStruct.W = W;
initStruct.A = ACon;
initStruct.F = FCon;
initStruct.L = LCon;
initStruct.AHat = S;
initStruct.BHat = M;
initStruct.eigH = eigH;
matrices = repmat(initStruct,N,1);


% Calculate wBar constant for i  = 1... Note that we need a for-loop to
wBar = zeros(N,1);
rho_i =  sqrt(max(eigH));
term1 = sqrt((N - k - 1)*gamma + lambda*alpha - betaVec(i+1));
term2 = max([sqrt((N - k - 1)*gamma + lambda*alpha - betaVec(i)), sqrt((lambda-lambda_p)*alpha)]);
wBar_i = (1/rho_i)*(term1 - term2);
wBar(i) = wBar_i;

% % Calculate ellStar for 1st timestep
% ellStar.coldStart = zeros(N,1); % initialize
% ellStar.warmStart = zeros(N,1);
% eta = zeros(N,1);
% eta_i = (condH - 1)/(condH + 1);
% eta(1) = eta_i;
% ellStar.coldStart(1) = (1/log(eta_i))*( log(sqrt(min(eigH))*sqrt(min(eigP))*epsilonBar_i) ...
%     - log(norm(sqrtm(invH)*G)*sqrt(VBar(i))));
% ellStar.warmStart(1) = ellStar.coldStart(1); % fill first as CS since no WS

% Run loop and find wBar
condHVec = zeros(N,1);
condHVec(1) = max(eigH)/min(eigH);
rhoVec = zeros(N,1);
rhoVec(1) = rho_i;
VBarVec = zeros(N,1);
VBarVec(1) = N*gamma + lambda*alpha - betaVec(1);
omegaVec = zeros(N,1);
sigmaVec = zeros(N,1);
for i=2:N
    k = i -1;
    N_k = N - k;
    
    % Generate the constant QP matrices (fixed horizon lengths)
    [H,G,W,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_k,Ad,Bd,lambda*P,Q,R,xmax,xmin,umax,umin);
    H = 0.5*(H+H');
    % Store matrices
%     invH = inv(H);
%     sqrtInvH = sqrtm(invH);
    eigH = eig(H);
    condH = max(eigH)/min(eigH);
    matrices(i).H = H;
    matrices(i).G = G;
    matrices(i).W = W;
    matrices(i).A = ACon;
    matrices(i).F = FCon;
    matrices(i).L = LCon;
    matrices(i).AHat = S;
    matrices(i).BHat = M;
    matrices(i).eigH = eigH;
    condHVec(i) = condH;
    
    % Calculate wBar constant for i  = 1... Note that we need a for-loop to
    if i < N
        rho_i =  sqrt(max(eigH));
        term1 = sqrt((N - k - 1)*gamma + lambda*alpha - betaVec(i+1));
        term2 = max([sqrt((N - k - 1)*gamma + lambda*alpha - betaVec(i)), sqrt((lambda-lambda_p)*alpha)]);
        wBar_i = (1/rho_i)*(term1 - term2);
        wBar(i) = wBar_i;
    else
        rho_i = sqrt(max(eigP));
        term1 = sqrt(lambda*alpha);
        term2 = max([sqrt(lambda*alpha - betaVec(i)), sqrt((lambda-lambda_p)*alpha)]);
        wBar_i = (1/rho_i)*(term1 - term2);
        wBar(i) = wBar_i;
    end
    rhoVec(i) = rho_i;
    VBarVec(i) = (N-k)*gamma + lambda*alpha - betaVec(i);
end

% Check the value of dBar against the minimum valaue of wBar
if min(wBar) < dMax
    warning('Provided dMax is too large. Setting to zero')
    dMax = 0;
end

% Calculate eBar 
nuVec = zeros(N,1);
eBarVec = zeros(N,1);
for i = 1:N
    nu_i = norm(sqrtm(matrices(i).W)*Bd);
    eBar_i = (rhoVec(i)/nu_i)*(wBar(i) - dMax);
    eBarVec(i) = eBar_i;
    nuVec(i) = nu_i;
end

% Calculate ellBar - for APGM
ellBarVec = zeros(N,1);
for i = 1:N
    % Unpack things
    eigH = matrices(i).eigH;
    kappa = max(eigH)/min(eigH);
    H = matrices(i).H;
    G = matrices(i).G;
    sqrtInvH = inv(sqrtm(H));

    % Calculate WS parameters
    if i == 1
        a_i = sqrt(min(eigH))*eBarVec(i)/(norm(sqrtInvH*G*X0));
        ellBarVec(i) = inverseEta_APGM(a_i,kappa);
    else
        omega_i = (1/sqrt(min(eigH)))*norm(sqrtInvH*G*Bd);
        sigma_i = (1/sqrt(min(eigH)))*norm(sqrtInvH*G);
        a_i = eBarVec(i)/((1+omega_i)*eBarVec(i-1) + sigma_i*dMax);
        ellBarVec(i) = inverseEta_APGM(a_i,kappa);
        omegaVec(i) = omega_i;
        sigmaVec(i) = sigma_i;
    end
end


% Store things
output.matrices = matrices;
output.VBarVec = VBarVec;
output.wBar = wBar;
output.condHVec = condHVec;
output.betaVec = betaVec;
output.rhoVec = rhoVec;
output.nuVec = nuVec;
output.eBarVec = eBarVec;
output.ellBarVec = ellBarVec;
output.omegaVec = omegaVec;
output.sigmaVec = sigmaVec;
end


function ell =  inverseEta_APGM(a,kappa)
   ell  = 2*(log(a) - log(sqrt(kappa)))/log(1 - (1/sqrt(kappa))) + 1;
end
