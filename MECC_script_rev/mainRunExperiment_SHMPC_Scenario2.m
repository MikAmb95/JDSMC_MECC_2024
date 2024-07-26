%Matlab simulation for paper submitted to MECC2024
%Author: Michele Ambrosino 

clc, close all,clear all

%folder with some utilities functions
addpath('extraFunctionV2_new') %
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))


%% * Set simulation parameters

%These are the three different configuration we used in the paper
Nn = [150 150 100];
Rr = [20 20*0.7 20];

for k = 1:3

% * Time vector
dt = 0.1;
N = Nn(k);
TEnd = N*dt;
t = linspace(0,TEnd,N+1);

% * Initial Conditions
X0 = [1;1]; 

% * State space matrices   
global A B Ad Bd

% State Vector (X): omega1
%                   omega2


% Continuous-time dynamics
n_spin = 1; % rad/s, base spin rate
J_ratio = 0.05; % J3/J 
lambda_const = n_spin*(1 - J_ratio); % rad/s^2, lambda = n*(J-J3)/J
Jt = 1; % Nm
J3 = J_ratio*Jt;
A = [0,lambda_const;
     -lambda_const, 0];
B = [0; 1/Jt];

% Get discrete state-space matrices
[Ad,Bd] = c2d(A,B,dt);
n = size(A,1);
m = size(B,2);


%% * Set control parameters and run the simulation

% Set horizon
T_MPC = dt;
N_MPC = N;

% * Set P, Q, R
Q = diag([1 1]); %  State error cost
R = Rr(k); % Control effort cost
[K,P] = dlqr(Ad,Bd,Q,R);

% * Set constraints 
umax = 0.2*ones(m,1);
umin = -umax;
xmax = ones(n,1)*1000; % just a placeholder, doesn't get used since we use PGM
xmin = -xmax;
MaxIter = 10000; % 7 and 6200 at 0.2 s, R = 0.1, 10 and 8000 at 0.15 s, 15 and 5500 at 0.1 s
xTol_MPC = 1e-8;  % ||x||_inf optimization tolerance for MPC optimizer
if mod(T_MPC,dt) ~= 0
    error('Mismatch in divisibility of T and dt')
end


%% Calculate alpha, d, and lambda

% -------------- Calculate constant d and alpha --------------
% Compute alpha according to Marco's formula - taking the minimum
% OmegaBar = { x | -Kx in U }
ACon = [-K; K];
bCon = [umax; -umax];
alpha = 1e12; % initialize
invP = inv(P);
for i = 1:length(bCon)
    b_i = bCon(i);
    a_i = ACon(i,:)';
    VNew = b_i^2/(a_i'*invP*a_i); % assumes x_eq = 0
    if VNew < alpha
        alpha = VNew;
    end
end
sqrtQ  = sqrtm(Q);
invsqrtQ = inv(sqrtQ);
eigVec = eig(invsqrtQ*P*invsqrtQ);
gamma = alpha/max(eigVec);

InitialNorm = X0'*P*X0;
alpha;



%% Calculate lambda and beta

% Get initial cost given lambda
lambda = 10;
[H,G,W,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_MPC,Ad,Bd,lambda*P,Q,R,xmax,xmin,umax,umin);
H = 0.5*(H+H');
% 
% Get n*m length vectors u_low and u_up
u_low = zeros(m*(N-1),1);
u_up = zeros(m*(N-1),1);
for j = 1:N
    u_low(1 + m*(j-1) : m + m*(j-1)) = umin;
    u_up(1 + m*(j-1) : m + m*(j-1)) = umax;
end

% Calculate true value function 
H_QP = H;
f_QP = G*X0;
[U_MPC,iterCount] = accelProjGradSolver(H_QP,f_QP,zeros(m*N,1),u_low,u_up,700,1e-12,0);
 
% Calculate cost
COST = U_MPC'*H*U_MPC + 2*U_MPC'*G*X0 + X0'*W*X0;
BOUND = N*gamma + lambda*alpha;
beta = BOUND - COST;


% Get initial controls and states
u_plot = reorderMPCVector(U_MPC,m,N);
XVec = S*X0 + M*U_MPC;
x_plot = reorderMPCVector(XVec,n,N);

if InitialNorm < alpha
    error('Initial state is in terminal region')
end
if beta < 0
    error('Negative beta')
end

%% Get shrinking horizon constants

betaVec = zeros(N,1);
betaVec(1) = beta;

% LINEAR BETA SEQUENCE
for jj = 2:N
    kk = jj - 1;
    betaVec(jj) = min([beta*(N - kk)/N, 0.99*(N-kk)*gamma + lambda*alpha]);
end




% * Add controller arguments to a structure
% controlArgs.P = Sf;
controlArgs.P = P;
controlArgs.Q = Q;
controlArgs.R = R; 
controlArgs.K = K;
controlArgs.N = N_MPC; % Horzion lengt
controlArgs.T = T_MPC; % Sampling period
controlArgs.MaxIter = MaxIter;
controlArgs.xTol = xTol_MPC;
controlArgs.umax = umax; % Control constraints
controlArgs.umin = umin;
controlArgs.alpha = alpha;
controlArgs.gamma = gamma;
controlArgs.lambda = lambda;
controlArgs.dMax = 0;
controlArgs.betaVec = betaVec;
controlArgs.X0 = X0;

% Run function for the first sequence
output1 = calcParamsSHMPC(controlArgs);
init1 = output1.ellBarVec(1);


%% Run the analysis with disturbances

% Use the linear one
controlArgs.betaVec = betaVec;

% Run offline analysis using 3 different disturbances
dBar  = (min(output1.wBar));
d1Bar = 9/10*(min(output1.wBar));
controlArgs.dMax = d1Bar;
output1_dist = calcParamsSHMPC(controlArgs);



%%  Run the online analysis
invsqrtH = inv(sqrtm(H));
eigH = eig(H);

% Online disturbance
dOnline = 0.0075; % 0.008; % 4*(min(output1.wBar))
controlArgs.dMax = dOnline;

controlArgs.V0 = COST;
controlArgs.U_Init = U_MPC;
kappa_i = cond(H);
eta_i = sqrt(kappa_i)* (1 - 1/sqrt(kappa_i))^(99/2);
controlArgs.InitError = eta_i*(1/sqrt(min(eigH)))*norm(invsqrtH*G*X0);


[X_ol,U_ol] = runShrinkingHorizonMPC_openLoop(X0,U_MPC,controlArgs);
COST_opt = U_ol*H*U_ol' + 2*U_ol*G*X0 + X0'*W*X0

iter = output1_dist.ellBarVec; 
output1_dist.ellBarVec = iter;
[X_fixed,U_fixed,iter_out] = runShrinkingHorizonMPC_fixedIter(X0,controlArgs,output1_dist);
COST_sub = U_fixed*H*U_fixed' + 2*U_fixed*G*X0 + X0'*W*X0
sum(iter_out)
%output1_dist.ellBarVec = iter_out;
[Gap,R_U,R_X,delta,L,psi,psi_bar,bound] = eval_R_v2(Ad,Bd,controlArgs,output1_dist,dBar,N);
Gap



res.xopt(k).X = X_ol;
res.x(k).X = X_fixed;
res.u(k).U = U_fixed;
res.uopt(k).U = U_ol;
res.cost(k) = COST_sub;
res.optcost(k) = COST_opt;
res.nIter(k).Iter = iter_out;
res.Gap(k) = Gap; 
res.alpha(k) = alpha; 
res.N(k) = N;
res.P(k).P = P;
end


%% Plotting

clear i

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesFontSize',23)
labelsize = 35;
legendsize = 22;
colorMatrix = [0 0 0;  135 135 135; 60 60 60; 55 126 184; 228 26 28; 77 175 74]./255;


figure(3)
% Plot the sublevel set
for i = 1:2
hold on
NGrid = 250;
x1 = linspace(-1,1,NGrid);
x2 = linspace(-1,1,NGrid);
[XX,YY]=meshgrid(x1,x2);
P = res.P(i).P;
Z= P(1,1)*XX.^2 + P(1,2).*XX.*YY + P(2,1).*XX.*YY+ P(2,2)*YY.^2;
[~,hContour] = contourf(XX,YY,-Z,-[res.alpha(i) res.alpha(i)]);
setColor(1,:) = [255 135 1]/255;
setColor(2,:) = [255 10 1]/255;
setColor(3,:) = [135 135 255]/255;
set(hContour,'facecolor',setColor(i,:),'edgecolor',setColor(i,:),'linewidth',1)
drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
hFills = hContour.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 0.4*255;   % default=255
end

end
hold on; box on; grid on
h1 = plot(res.x(1).X(1,:),res.x(1).X(2,:),'--o');
h2 = plot(res.x(2).X(1,:),res.x(2).X(2,:),'color',[1 0 0],'linestyle','-.');
h3 = plot(res.x(3).X(1,:),res.x(3).X(2,:),'color',[0 1 0],'linestyle','--');
%h4 = plot(res.xX_ol(1,:),X_ol(2,:),'color',colorMatrix(1,:));
xlabel('$\dot{\theta}_1 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)
ylabel('$\dot{\theta}_2 \; [{{rad}\over{s}}]$','interpreter','latex','fontsize',labelsize)

plot(res.x(3).X(1,end),res.x(3).X(2,end),'r>','LineWidth',15);